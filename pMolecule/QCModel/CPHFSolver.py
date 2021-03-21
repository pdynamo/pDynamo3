"""A class for solving the CPHF equations and obtaining the Z-matrix."""

#-----------------------------------------------------------------------------------------------------------------------------------
#
# . In the CPHF equations there are five groups of orbitals:
#
#       inActive double, active double, fractional, active virtual, inActive virtual.
#
#   Their interactions make up 15 blocks ( ( 5 * ( 5 + 1 ) ) / 2 ).
#
# . Variables:
#
#   Non-redundant - all orbital pairs, xy, which have different occupancies (eight blocks).
#   Redundant     - all orbital pairs involving at least one active orbital that have the same occupancy (five blocks).
#
#   The two inActive blocks (double/double and virtual/virtual) do not enter.
#
#-----------------------------------------------------------------------------------------------------------------------------------

from  math                      import fabs
from  pCore                     import AttributableObject          , \
                                       DataType
from  pScientific.Arrays        import Array                       , \
                                       StorageType
from  pScientific.LinearAlgebra import CGLinearEquationSolver      , \
                                       CGLinearEquationSolverState
from .CICPHF                    import CICPHF_ApplyCPHFMatrix      , \
                                       CICPHF_CalculateCPHFVectors , \
                                       CICPHF_Transform

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
_OccupancyTolerance = 1.0e-06

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class CPHFSolver ( AttributableObject ):
    """A class for solving the CPHF equations."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "aDiagonal"                 : None ,
                             "energies"                  : None ,
                             "fCore"                     : None ,
                             "indicesNR"                 : None ,
                             "indicesR"                  : None ,
                             "moTEI234"                  : None ,
                             "nActive"                   : 0    ,
                             "nCore"                     : 0    ,
                             "nOrbitals"                 : 0    ,
                             "numberDegenerateRedundant" : 0    ,
                             "numberNonRedundant"        : 0    ,
                             "numberRedundant"           : 0    ,
                             "report"                    : None ,
                             "occupancies"               : None ,
                             "occupancyTolerance"        : _OccupancyTolerance ,
                             "onePDM"                    : None ,
                             "onePDMMO"                  : None ,
                             "orbitals"                  : None ,
                             "preconditioner"            : None ,
                             "qNR"                       : None ,
                             "qR"                        : None ,
                             "rhs"                       : None ,
                             "solver"                    : None ,
                             "solverState"               : None ,
                             "target"                    : None ,
                             "twoElectronIntegrals"      : None ,
                             "twoPDM"                    : None ,
                             "warnings"                  : None ,
                             "work1"                     : None ,
                             "work2"                     : None ,
                             "zMatrix"                   : None } )

    # . The following two methods are needed by the CG solver. 
    def ApplyMatrix ( self, x, y ):
        """Apply the A matrix to x and put in y."""
        CICPHF_ApplyCPHFMatrix ( self.numberNonRedundant   ,
                                 self.indicesNR            ,
                                 self.numberNonRedundant   ,
                                 self.indicesNR            ,
                                 self.aDiagonal            ,
                                 x                         ,
                                 self.orbitals             ,
                                 self.twoElectronIntegrals ,
                                 self.work1                ,
                                 self.work2                ,
                                 y                         )

    def ApplyPreconditioner ( self, x, y ):
        """Apply the diagonal preconditioner to x and put in y."""
        x.CopyTo ( y )
        y.Multiply ( self.preconditioner )

    @classmethod
    def FromTarget ( selfClass, target ):
        """Constructor from target."""
        self           = selfClass ( )
        self.target    = target
        self.nActive   = target.qcModel.activeOrbitals
        self.nCore     = target.qcState.coreOrbitals
        self.nOrbitals = len ( target.qcState.orbitalBases )
        return self

    def GatherTargetAliases ( self ):
        """Gather the target aliases."""
        scratch                   = self.target.scratch
        node                      = scratch.ci
        self.energies             = scratch.orbitalsP.energies
        self.fCore                = node.fCore
        self.moTEI234             = node.moTEI234
        self.occupancies          = scratch.orbitalsP.occupancies
        self.onePDM               = node.onePDM
        self.onePDMMO             = node.onePDMMO
        self.orbitals             = scratch.orbitalsP.orbitals
        self.twoElectronIntegrals = scratch.twoElectronIntegrals
        self.twoPDM               = node.twoPDM

    def SetUp ( self ):
        """Set up the CPHF calculation."""
        # . Initialization.
        numberActive        = self.nActive  
        numberCore          = self.nCore    
        numberOrbitals      = self.nOrbitals
        self.warnings       = []
        # . Reset the target aliases.
        self.GatherTargetAliases ( )
        # . Determine some counters.
        fractionalOccupancy = 0.0
        numberActiveDouble  = 0
        numberActiveVirtual = 0
        numberFractional    = 0
        isOK                = True
        for o in self.occupancies[numberCore:numberCore+numberActive]:
            if   fabs ( o - 2.0e+00 ) < self.occupancyTolerance: numberActiveDouble  += 1
            elif fabs ( o           ) < self.occupancyTolerance: numberActiveVirtual += 1
            else:
                if ( numberFractional > 0 ) and ( fabs ( o - fractionalOccupancy ) > self.occupancyTolerance ): isOK = False
                fractionalOccupancy  = o
                numberFractional    += 1
        if not isOK: self.warnings.append ( "Unequal active orbital fractional occupancies in gradients." )
        numberInActiveVirtual = numberOrbitals - ( numberCore + numberActive )
        # . Number of variables. */
        numberDegenerateRedundant = 0
        numberNonRedundant        =  numberCore           * ( numberFractional + numberActiveVirtual + numberInActiveVirtual ) + \
                                     numberActiveDouble   * ( numberFractional + numberActiveVirtual + numberInActiveVirtual ) + \
                                     numberFractional     * (                    numberActiveVirtual + numberInActiveVirtual )
        numberRedundant           =  numberCore * numberActiveDouble + \
                                   ( numberActiveDouble   * ( numberActiveDouble  - 1 ) ) // 2 + \
                                   ( numberFractional     * ( numberFractional    - 1 ) ) // 2 + \
                                   ( numberActiveVirtual  * ( numberActiveVirtual - 1 ) ) // 2 + \
                                     numberActiveVirtual  * numberInActiveVirtual
        # . Allocate variable arrays. */
        if self.numberNonRedundant != numberNonRedundant:
            self.indicesNR      = Array.WithExtents ( numberNonRedundant , 2, dataType = DataType.Integer )
            self.aDiagonal      = Array.WithExtent  ( numberNonRedundant )
            self.qNR            = Array.WithExtent  ( numberNonRedundant )
            self.rhs            = Array.WithExtent  ( numberNonRedundant )
            self.preconditioner = Array.WithExtent  ( numberNonRedundant )
        if self.numberRedundant != numberRedundant:
            self.indicesR       = Array.WithExtents ( numberRedundant , 2, dataType = DataType.Integer )
            self.qR             = Array.WithExtent  ( numberRedundant )
        # . Constant arrays.
        if self.work1   is None: self.work1   = Array.WithExtent ( numberOrbitals, storageType = StorageType.Symmetric )
        if self.work2   is None: self.work2   = Array.WithExtent ( numberOrbitals, storageType = StorageType.Symmetric )
        if self.zMatrix is None: self.zMatrix = Array.WithExtent ( numberOrbitals, storageType = StorageType.Symmetric )
        # . Save the counters. */
        self.numberDegenerateRedundant = numberDegenerateRedundant
        self.numberNonRedundant        = numberNonRedundant
        self.numberRedundant           = numberRedundant

    def Solve ( self ):
        """Solve the CPHF equations for zNR and form the Z-matrix."""
        # . Calculate the CPHF vectors.
        ( self.numberDegenerateRedundant, self.numberRedundant ) = CICPHF_CalculateCPHFVectors ( self.nActive                   ,
                                                                                                 self.nCore                     ,
                                                                                                 self.nOrbitals                 ,
                                                                                                 self.twoElectronIntegrals      ,
                                                                                                 self.twoPDM                    ,
                                                                                                 self.energies                  ,
                                                                                                 self.occupancies               ,
                                                                                                 self.orbitals                  ,
                                                                                                 self.moTEI234                  ,
                                                                                                 self.fCore                     ,
                                                                                                 self.onePDM                    ,
                                                                                                 self.onePDMMO                  ,
                                                                                                 self.work1                     ,
                                                                                                 self.work2                     ,
                                                                                                 self.numberDegenerateRedundant ,
                                                                                                 self.numberNonRedundant        ,
                                                                                                 self.numberRedundant           ,
                                                                                                 self.indicesNR                 ,
                                                                                                 self.indicesR                  ,
                                                                                                 self.aDiagonal                 ,
                                                                                                 self.qNR                       ,
                                                                                                 self.qR                        ,
                                                                                                 self.preconditioner            )
        # . Initialize target arrays. 
        self.qNR.CopyTo ( self.rhs )  # . RHS. 
        self.qNR.Set    ( 0.0      )  # . Initial guess at solution.
        # . Set up the solver.
        if self.solver is None:
            self.solver      = CGLinearEquationSolver ( )
            self.solverState = CGLinearEquationSolverState.FromTarget ( self, self.rhs, self.qNR, doPreconditioning = True )
        # . Solve. 
        self.report = self.solver.Solve ( self.solverState )
        if not self.report["Is Converged"]: self.warnings.append ( "CPHF calculation not converged." )
        # . Extract Z and convert to AO basis. 
        CICPHF_Transform ( self.numberNonRedundant ,
                           self.indicesNR          ,
                           self.qNR                ,
                           self.numberRedundant    ,
                           self.indicesR           ,
                           self.qR                 ,
                           self.orbitals           ,
                           False                   ,
                           self.work1              ,
                           self.zMatrix            )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
