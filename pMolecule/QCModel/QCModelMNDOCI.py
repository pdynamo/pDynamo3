"""A CI MNDO QC model."""

import math

from   enum                      import Enum
from   pCore                     import Align                 , \
                                        Clone                 , \
                                        logFile               , \
                                        LogFileActive
from   pScientific               import Units
from   pScientific.Arrays        import Array                 , \
                                        SparseSymmetricMatrix , \
                                        StorageType
from   pScientific.LinearAlgebra import EigenPairs
from  .CIConfigurationContainer  import CIConfigurationContainer
from  .CIFourIndexTransformation import CIFourIndexTransformation_Make
from  .CISparseSolver            import CISparseSolver
from  .CPHFSolver                import CPHFSolver
from  .ElectronicState           import OccupancyType         , \
                                        SpinMultiplicity
from  .FockConstruction          import FockConstruction_MakeFromTEIs
from  .QCModelMNDO               import QCModelMNDO           , \
                                        QCModelMNDOState
from  .QCModelError              import QCModelError
from ..EnergyModel               import EnergyClosurePriority

# . Note that the CI module is inappropriate if there are QC/QC image interactions.

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
class CIDiagonalization ( Enum ):
    """The type of CI diagonalization."""
    Dense  = 1
    Sparse = 2

class CIMethod ( Enum ):
    """The type of CI method."""
    Doubles        = 1
    Full           = 2
    Singles        = 3
    SinglesDoubles = 4
    User           = 5

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCModelMNDOCIState ( QCModelMNDOState ):
    """A QC MNDO model state."""

    _attributable = dict ( QCModelMNDOState._attributable )
    _attributable.update ( { "activeAlpha"     : 0    ,
                             "activeBeta"      : 0    ,
                             "configurations"  : None ,
                             "coreElectrons"   : 0    ,
                             "coreOrbitals"    : 0    } ) # . Dimensionless.

    def SummaryItems ( self ):
        """Summary items."""
        items = super ( QCModelMNDOCIState, self ).SummaryItems ( )
        items.extend ( [ ( "CI Matrix Sparsity" , "{:.1%}".format ( self.configurations.CIMatrixSparsity ( )[-1] / 100.0 ) ) ,
                         ( "Configurations"     , "{:d}".format ( len ( self.configurations ) ) ) ,
                         ( "Core Electrons"     , "{:d}".format ( self.coreElectrons          ) ) ,
                         ( "Core Orbitals"      , "{:d}".format ( self.coreOrbitals           ) ) ] )
        return items

    @property
    def spinDensity ( self ):
        """The spin density for property calculations."""
        return self.target.scratch.ci.Get ( "spinDensity", None )

    @property
    def totalDensity ( self ):
        """The total density for property calculations."""
        return self.target.scratch.ci.dTotal

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCModelMNDOCI ( QCModelMNDO ):
    """A CI MNDO QC model."""

    _attributable = dict ( QCModelMNDO._attributable )
    _classLabel   = "MNDO CI QC Model"
    _stateObject  = QCModelMNDOCIState
    _summarizable = dict ( QCModelMNDO._summarizable )
    _attributable.update ( { "activeElectrons"              : 0                       ,
                             "activeOrbitals"               : 0                       ,
                             "checkSparseDiagonalization"   : False                   ,
                             "ciDiagonalization"            : CIDiagonalization.Dense ,
                             "ciMethod"                     : CIMethod.Full           ,
                             "degeneracyTolerance"          : 1.0e-03                 , # . Hartrees.
                             "fractionalOccupancyTolerance" : 1.0e-06                 , # . Dimensionless.
                             "identifyRootSpin"             : True                    ,
                             "microStates"                  : None                    , # . User input.
                             "minimalMultiplicity"          : 1                       ,
                             "multiplicity"                 : 1                       ,
                             "numberOfStates"               : 20                      ,
                             "requiredRoot"                 : 0                       ,
                             "requiredSpin"                 : 0.0                     ,
                             "spinTolerance"                : 1.0e-02                 } ) # . Dimensionless.
    _summarizable.update ( { "activeElectrons"              :   "Active Electrons"                          ,
                             "activeOrbitals"               :   "Active Orbitals"                           ,
                             "checkSparseDiagonalization"   :   "Check Sparse Diagonalization"              ,
                             "ciDiagonalization"            :   "CI Diagonalization"                        ,
                             "ciMethod"                     :   "CI Method"                                 ,
                             "degeneracyTolerance"          : ( "Degeneracy Tolerance"           , "{:g}" ) ,
                             "fractionalOccupancyTolerance" : ( "Fractional Occupancy Tolerance" , "{:g}" ) ,
                             "identifyRootSpin"             :   "Identify Root Spin"                        ,
                             "minimalMultiplicity"          :   "Minimal Multiplicity"                      ,
                             "multiplicity"                 :   "Multiplicity"                              ,
                             "numberOfStates"               :   "Number Of States"                          ,
                             "requiredRoot"                 :   "Required Root"                             ,
                             "spinTolerance"                : ( "Spin Tolerance"                 , "{:g}" ) } )

    def BuildModel ( self, target, qcSelection = None ):
        """Build the model."""
        state = super ( QCModelMNDOCI, self ).BuildModel ( target, qcSelection = qcSelection )
        self.CheckInput         ( target )
        self.CheckVersusSystem  ( target )
        self.MakeConfigurations ( target )

    def CheckInput ( self, target ):
        """Check input options."""
        # . For generation of configurations use the minimal multiplicity as paired electrons can contribute to high-spin states.
        # . For root checking, use the full multiplicity.
        # . Microstates are checked only when generating configurations.
        state             = target.qcState
        multiplicity      = min ( self.minimalMultiplicity, self.multiplicity )
        state.activeAlpha = ( self.activeElectrons + multiplicity - 1 ) // 2
        state.activeBeta  = ( self.activeElectrons - multiplicity + 1 ) // 2
        self.requiredSpin = 0.25 * float ( ( self.multiplicity - 1 ) * ( self.multiplicity + 1 ) )
        isOK = ( self.activeElectrons     >= 0                                        ) and \
               ( self.activeElectrons     <= ( 2 * self.activeOrbitals              ) ) and \
               ( self.minimalMultiplicity >  0                                        ) and \
               ( self.multiplicity        >  0                                        ) and \
               ( self.requiredRoot        >= 0                                        ) and \
               ( self.activeElectrons     == ( state.activeAlpha + state.activeBeta ) ) and \
               ( state.activeAlpha        >= 0                                        ) and \
               ( state.activeBeta         >= 0                                        ) and \
               ( state.activeAlpha        <= self.activeOrbitals                      )
        if not isOK: raise QCModelError ( "Invalid CI model input options." )

    def CheckVersusSystem ( self, target ):
        """Check versus system parameters."""
        state               = target.qcState
        electrons           = int ( round ( state.alphaCharge + state.betaCharge ) )
        orbitals            = len ( state.orbitalBases )
        state.coreElectrons = electrons - self.activeElectrons
        state.coreOrbitals  = state.coreElectrons // 2
        activeElectrons     = electrons - ( 2 * state.coreOrbitals )
        virtualOrbitals     = orbitals - ( state.coreOrbitals + self.activeOrbitals )
        isOK = ( target.electronicState.isSpinRestricted                                  ) and \
               ( target.electronicState.occupancyType != OccupancyType.FractionalVariable ) and \
               ( state.coreElectrons   >=   0                                             ) and \
               ( state.coreElectrons   == ( 2 * state.coreOrbitals                      ) ) and \
               ( self.activeElectrons  == activeElectrons                                 ) and \
               ( virtualOrbitals       >=   0                                             )
        if not isOK: raise QCModelError ( "The CI model is incompatible with the system." )

    def CIEnergy ( self, target ):
        """The CI energy."""
        # . Set up the calculation.
        ( node, eHF, occupancies, orbitals, oneElectronMatrix, twoElectronIntegrals, energyTerms, warnings ) = self.CIInitialize ( target )
        state = target.qcState
        # . Determine the core quantities.
        # . Density.
        node.dCore.MakeFromEigenSystem ( state.coreOrbitals, occupancies, orbitals )
        # . Energies and Fock matrix.
        eOE = node.dCore.TraceOfProduct ( oneElectronMatrix )
        eTE = FockConstruction_MakeFromTEIs ( node.dCore, None, twoElectronIntegrals, 1.0, node.fCore, None )
        node.fCore.Add ( oneElectronMatrix )
        # . Baseline energy = eCore - eHF.
        eCore     = eOE + eTE + state.energyBaseLine
        eBaseLine = eCore * Units.Energy_Hartrees_To_Kilojoules_Per_Mole - eHF
        # . Fock matrix in the active MO basis.
        node.fCore.Transform ( node.activeMOs, node.fCoreMO )
        # . Calculate the CI matrix (column-major).
        moTEIs = CIFourIndexTransformation_Make ( node.activeMOs, twoElectronIntegrals, node )
        state.configurations.MakeCIMatrix ( node.fCoreMO, moTEIs, node.ciMatrixDense, node.ciMatrixSparse )
        # . Diagonalization.
        if self.ciDiagonalization == CIDiagonalization.Dense:
            EigenPairs ( node.ciMatrixDense, node.ciEnergies, node.ciVectors, columnMajor = True, preserveInput = False, upper = self.numberOfStates )
        else:
            node.ciMatrixSparse.MakeDiagonalPreconditioner ( node.preconditioner )
            report = node.sparseSolver.Solve ( )
            if self.checkSparseDiagonalization:
                EigenPairs ( node.ciMatrixDense, node.ciEnergiesReference, None, columnMajor = True, preserveInput = False, upper = self.numberOfStates )
                report.update ( node.sparseSolver.CheckSolution ( referenceEigenvalues = node.ciEnergiesReference ) )
            target.scratch.qcEnergyReport.update ( report )
        # . Convert energies to proper units and adjust with respect to HF energy.
        node.ciEnergies.Scale ( Units.Energy_Hartrees_To_Kilojoules_Per_Mole )
        node.ciEnergies.Increment ( eBaseLine )
        # . Find the state spins.
        state.configurations.StateSpins ( node.ciVectors, node.ciSpins ) # . ciSpins determines n.
        # . Find the required root.
        rootNumber = 0
        if self.identifyRootSpin:
            goodSpins = [ i for ( i, s ) in enumerate ( node.ciSpins ) if math.fabs ( s - self.requiredSpin ) < self.spinTolerance ]
            if self.requiredRoot < len ( goodSpins ): rootNumber = goodSpins[self.requiredRoot]
            else: node.warnings.append ( "Unable to identify root of required spin." )
        else: rootNumber = self.requiredRoot
        # . Save the final data.
        ciEnergy = node.ciEnergies[rootNumber  ]
        ciVector = node.ciVectors [rootNumber,:]
        energyTerms["QC CI"] = ciEnergy
        node.ciEnergy = ciEnergy
        node.ciRoot   = rootNumber
        node.ciVector = ciVector
        # . Make the CI densities.
        # . dTotal and dSpin are general.
        # . onePDMHF, onePDMMO(t), onePDM and dCore are only needed for the gradients.
        # . onePDMMOs is scratch.
        state.configurations.MakeDensities ( ciVector, node.onePDMMO, node.onePDMMOs, node.twoPDM )
        node.onePDMMO.Transform  ( node.activeMOs, node.onePDM, useTranspose = True )
        node.onePDMMOs.Transform ( node.activeMOs, node.dSpin , useTranspose = True )
        node.onePDM.CopyTo ( node.dTotal ) ; node.dTotal.Add ( node.dCore )

    def CIInitialize ( self, target ):
        """Set up a CI energy calculation."""
        # . Get already calculated quantities.
        scratch              = target.scratch
        energyTerms          = scratch.energyTerms
        eHF                  = energyTerms["QC Electronic"]
        source               = target.scratch.orbitalsP
        energies             = source.energies
        occupancies          = source.occupancies
        orbitals             = source.orbitals
        oneElectronMatrix    = scratch.oneElectronMatrix
        twoElectronIntegrals = scratch.twoElectronIntegrals
        # . Check the HF orbitals for degenerate energies and fractional occupancies at boundaries of active space.
        state    = target.qcState
        nActive  = self.activeOrbitals
        nCore    = state.coreOrbitals
        warnings = []
        if nCore > 0:
            if math.fabs ( energies[nCore] - energies[nCore-1] ) < self.degeneracyTolerance: 
                warnings.append ( "Degenerate orbitals at boundary of core and active spaces." )
            if math.fabs ( occupancies[nCore-1] - 2.0 ) > self.fractionalOccupancyTolerance:
                warnings.append ( "Fractionally occupied core orbitals." )
        n = nCore + self.activeOrbitals
        if n < len ( energies ):
            if math.fabs ( energies[n] - energies[n-1] ) < self.degeneracyTolerance: 
                warnings.append ( "Degenerate orbitals at boundary of active and virtual spaces." )
            if math.fabs ( occupancies[n] ) > self.fractionalOccupancyTolerance:
                warnings.append ( "Fractionally occupied virtual orbitals." )
        # . Set up the CI.
        doSetUp = ( not hasattr ( scratch, "ci" ) )
        node    = scratch.GetSetNode ( "ci" )
        if doSetUp:
            nBasis          = len ( state.orbitalBases )
            nConfigurations = len ( state.configurations )
            node.twoPDM     = Array.WithExtent  ( nActive, storageType = StorageType.DoubleSymmetric )
            node.dCore      = Array.WithExtent  ( nBasis , storageType = StorageType.Symmetric       )
            node.dSpin      = Array.WithExtent  ( nBasis , storageType = StorageType.Symmetric       )
            node.dTotal     = Array.WithExtent  ( nBasis , storageType = StorageType.Symmetric       )
            node.fCore      = Array.WithExtent  ( nBasis , storageType = StorageType.Symmetric       )
            node.fCoreMO    = Array.WithExtent  ( nActive, storageType = StorageType.Symmetric       )
            node.onePDM     = Array.WithExtent  ( nBasis , storageType = StorageType.Symmetric       )
            node.onePDMMO   = Array.WithExtent  ( nActive, storageType = StorageType.Symmetric       )
            node.onePDMMOs  = Array.WithExtent  ( nActive, storageType = StorageType.Symmetric       )
            node.ciEnergies = Array.WithExtent  ( self.numberOfStates                                )
            node.ciSpins    = Array.WithExtent  ( self.numberOfStates                                )
            node.ciVectors  = Array.WithExtents ( self.numberOfStates , nConfigurations              )
            # . Diagonalization.
            ciMatrixDense  = None
            ciMatrixSparse = None
            if self.ciDiagonalization == CIDiagonalization.Sparse:
                ( nonZero, sparsity ) = state.configurations.CIMatrixSparsity ( )
                ciMatrixSparse        = SparseSymmetricMatrix.WithExtentAndSize ( nConfigurations, nonZero )
                node.preconditioner   = Array.WithShape ( [ nConfigurations ] )
                sparseSolver          = CISparseSolver.FromArrays ( node.ciEnergies, node.ciVectors, ciMatrixSparse, node.preconditioner )
                node.sparseSolver     = sparseSolver
                if self.checkSparseDiagonalization:
                    node.ciEnergiesReference = Array.WithShape ( [ self.numberOfStates ] )
            if ( self.ciDiagonalization == CIDiagonalization.Dense ) or self.checkSparseDiagonalization:
                ciMatrixDense = Array.WithExtent ( nConfigurations, storageType = StorageType.Symmetric )
            node.ciMatrixDense  = ciMatrixDense
            node.ciMatrixSparse = ciMatrixSparse
        # . Save some other quantities.
        node.activeMOs = orbitals[:,nCore:nCore+nActive]     # . Alias only.
        node.onePDMHF = target.scratch.onePDMP.density # . Alias only.
        node.warnings = warnings
        # . Adjustments for QC/MM models.
        # . This is a fudge and needs cleaning up and generalizing.
        # . One-PDM is used temporarily as scratch space instead of destroying One-Electron Matrix.
        multipoleOrder = scratch.Get ( "qcmmMultipoleOrder", None )
        qcmmPotentials = scratch.Get ( "qcmmPotentials"    , None )
        if qcmmPotentials is not None:
            oneElectronMatrix.CopyTo ( node.onePDM )
            eHF += energyTerms["QC/MM Electrostatic"]
            # . Assume density-based QC/MM model.
            if multipoleOrder is None:
                node.onePDM.Add ( qcmmPotentials )
            # . Assume multipole-based QC/MM model.
            else:
                n    = len ( state.nuclearCharges )
                eHF -= ( qcmmPotentials[0:n].Dot ( state.nuclearCharges ) * Units.Energy_Hartrees_To_Kilojoules_Per_Mole ) # . To avoid double-counting of nuclei.
                self.multipoleEvaluator.FockMultipoleDerivatives ( target, multipoleOrder, qcmmPotentials, node.onePDM )
            oneElectronMatrix = node.onePDM
        return ( node, eHF, occupancies, orbitals, oneElectronMatrix, twoElectronIntegrals, energyTerms, warnings )

    def CIVectorsTable ( self, target, log = logFile, numberOfVectors = 20, vectorsPerRow = 8 ):
        """Output a table of CI vectors."""
        if LogFileActive ( log ):
            # . Gather data.
            node        = target.scratch.ci
            state       = target.qcState
            c           = state.configurations.numberOfConfigurations
            m           = state.configurations.numberOfActiveOrbitals
            n           = min ( self.numberOfStates, numberOfVectors )
            ciEnergies  = node.ciEnergies
            ciSpins     = node.ciSpins
            ciVectors   = node.ciVectors
            microStates = state.configurations.microStates
            # . Write the vectors.
            if n > 0:
                # . Find the number of tables to write.
                nTables = n // vectorsPerRow
                if ( n % vectorsPerRow > 0 ): nTables += 1
                # . Write the tables.
                for iTable in range ( nTables ):
                    # . Number of vectors to write.
                    nColumns = vectorsPerRow
                    if iTable == ( nTables - 1 ): nColumns = n - iTable * vectorsPerRow
                    # . Find width of first column.
                    nWidth = max ( 10, 2 * ( m + 1 ) )
                    # . Start tables.
                    columns = [ 6, nWidth ] + nColumns * [ 15 ]
                    table   = log.GetTable ( columns = columns )
                    table.Start ( )
                    table.Title ( "CI Vectors" )
                    table.Heading ( "" )
                    table.Heading ( "" )
                    for i in range ( nColumns ): table.Heading ( "{:d}".format ( i+1+iTable*vectorsPerRow ) )
                    table.Entry ( "Energy", align = Align.Left )
                    table.Entry ( "" )
                    for i in range ( nColumns ):
                        table.Entry ( "{:.6f}".format ( ciEnergies[i+iTable*vectorsPerRow] ) )
                    table.Entry ( "<S^2>", align = Align.Left )
                    table.Entry ( "" )
                    for i in range ( nColumns ):
                        table.Entry ( "{:.6f}".format ( ciSpins[i+iTable*vectorsPerRow] ) )
                    table.Entry ( "S", align = Align.Left )
                    table.Entry ( "" )
                    for i in range ( nColumns ):
                        table.Entry ( "{:.6f}".format ( 0.5 * ( -1.0 + math.sqrt ( 1.0 + 4.0 * ciSpins[i+iTable*vectorsPerRow] ) ) ) )
                    for i in range ( c ):
                        table.Entry ( "{:d}".format ( i+1 ) )
                        table.Entry ( microStates[i][0:m] + " " + microStates[i][m:] )
                        for j in range ( nColumns ):
                            table.Entry ( "{:.6f}".format ( ciVectors[j+iTable*vectorsPerRow,i] ) )
                    table.Stop ( )

    def CIWavefunctionSummary ( self, target, log = logFile, numberOfCoefficients = 100 ):
        """Write a summary of the CI wavefunction."""
        def _Sort ( x ): return math.fabs ( x[0] )
        if LogFileActive ( log ):
            # . Gather data.
            state       = target.qcState
            m           = state.configurations.numberOfActiveOrbitals
            n           = min ( state.configurations.numberOfConfigurations, numberOfCoefficients )
            microStates = state.configurations.microStates
            ciVector    = target.scratch.ci.ciVector
            items       = sorted ( [ ( c, i ) for ( i, c ) in enumerate ( ciVector ) ], key = _Sort, reverse = True )
            # . Write the table.
            doneFundamental = False
            columns = [ 7, 7, 14 ] + 2 * [ max ( len ( "Alphas" ), m ) + 2 ]
            table   = log.GetTable ( columns = columns )
            table.Start ( )
            table.Title ( "CI Wavefunction Coefficients" )
            table.Heading ( "Order"       )
            table.Heading ( "Conf."       )
            table.Heading ( "Coefficient" )
            table.Heading ( "Alphas"      )
            table.Heading ( "Betas"       )
            for ( order, ( c, index ) ) in enumerate ( items[0:n] ):
                table.Entry   ( "{:d}".format     ( order  ) )
                table.Entry   ( "{:d}".format     ( index  ) )
                table.Entry   ( "{:12.6f}".format ( c      ) )
                table.Entry   ( "{:s}".format     ( microStates[index][0:m] ) )
                table.Entry   ( "{:s}".format     ( microStates[index][m: ] ) )
                if index == 0: doneFundamental = True
            if not doneFundamental:
                table.Heading ( "HF Ground State", columnSpan = len ( columns ) )
                table.Entry   ( "{:d}".format     ( n           ) )
                table.Entry   ( "{:d}".format     ( 0           ) )
                table.Entry   ( "{:12.6f}".format ( ciVector[0] ) )
                table.Entry   ( "{:s}".format     ( microStates[0][0:m] ) )
                table.Entry   ( "{:s}".format     ( microStates[0][m: ] ) )
            table.Stop ( )

    def EnergyClosureGradients ( self, target ):
        """Gradient energy closure."""
        def a ( ): self.integralEvaluator.ResonanceGradients            ( target, doCI = True )
        def b ( ): self.integralEvaluator.ElectronNuclearTEIGradientsCI ( target )
        return [ ( EnergyClosurePriority.QCGradients, a, "QC Resonance Gradients"            ) ,
                 ( EnergyClosurePriority.QCGradients, b, "QC Electron-Nuclear/TEI Gradients" ) ] # . Note superclass closures not used.

    def EnergyClosures ( self, target ):
        """Return energy closures."""
        def a ( ): self.CIEnergy           ( target )
        def b ( ): self.SolveCPHFEquations ( target )
        closures = super ( QCModelMNDOCI, self ).EnergyClosures ( target )
        closures.extend ( [ ( EnergyClosurePriority.QCPostEnergy  , a, "QC CI Energy"         ) ,
                            ( EnergyClosurePriority.QCPreGradients, b, "QC CI CPHF Equations" ) ] )
        return closures

    def EnergyFinalize ( self, target ):
        """Energy finalization."""
        super ( QCModelMNDOCI, self ).EnergyFinalize ( target )
        self.EnergySummary ( target, log = target.scratch.log )

    def EnergySummary ( self, target, log = logFile, numberAbove = 10 ):
        """Summary of CI energy calculation."""
        if LogFileActive ( log ):
            node = target.scratch.ci
            # . Details of configurations.
            ciEnergies = node.ciEnergies
            ciSpins    = node.ciSpins
            table      = log.GetTable ( columns = [ 8 ] + 4 * [ 18 ] )
            table.Start ( )
            table.Title ( "CI Configurations" )
            table.Heading ( "Conf."  )
            table.Heading ( "Energy" )
            table.Heading ( "Type"   )
            table.Heading ( "<S^2>"  )
            table.Heading ( "S"      )
            for i in range ( min ( self.numberOfStates, node.ciRoot + numberAbove ) ):
                e     = ciEnergies[i]
                x     = ciSpins   [i]
                y     = 0.5 * ( -1.0 + math.sqrt ( 1.0 + 4.0 * x ) )
                j     = int ( round ( 2.0 * y + 1.0 ) )
                try   : label = SpinMultiplicity ( j ).name
                except: label = "{:d}".format    ( j )
                table.Entry ( "{:d}".format   ( i ) )
                table.Entry ( "{:.6f}".format ( e ) )
                table.Entry ( label )
                table.Entry ( "{:.5f}".format ( x ) )
                table.Entry ( "{:.5f}".format ( y ) )
            table.Stop ( )
            # . Warnings.
            warnings = node.warnings
            if len ( warnings ) > 0:
                n = max ( [ len ( s ) for s in warnings ] )
                table = log.GetTable ( columns = [ n + 2 ] )
                table.Start ( )
                table.Title ( "*** CI Warnings ***" )
                for s in sorted ( warnings ): table.Entry ( s, align = Align.Center )
                table.Stop ( )
            # . Final energy.
            log.Paragraph ( "Selected State CI Energy = {:.8g} kJ/mol.".format ( node.ciEnergy ) )

    def MakeConfigurations ( self, target ):
        """Make the configurations."""
        state = target.qcState
        if self.ciMethod == CIMethod.User:
            state.configurations = CIConfigurationContainer.FromMicroStates    ( self.microStates, self.activeElectrons, self.activeOrbitals )
        elif self.ciMethod == CIMethod.Doubles:
            state.configurations = CIConfigurationContainer.MakeDoubles        ( state.activeAlpha, state.activeBeta, self.activeOrbitals )
        elif self.ciMethod == CIMethod.Full:
            state.configurations = CIConfigurationContainer.MakeFull           ( state.activeAlpha, state.activeBeta, self.activeOrbitals )
        elif self.ciMethod == CIMethod.Singles:
            state.configurations = CIConfigurationContainer.MakeSingles        ( state.activeAlpha, state.activeBeta, self.activeOrbitals )
        elif self.ciMethod == CIMethod.SinglesDoubles:
            state.configurations = CIConfigurationContainer.MakeSinglesDoubles ( state.activeAlpha, state.activeBeta, self.activeOrbitals )
        if len ( state.configurations ) <= 0: raise QCModelError ( "There are no CI configurations." )
        if self.numberOfStates <= 0: self.numberOfStates = len ( state.configurations )
        else:                        self.numberOfStates = min ( len ( state.configurations ), self.numberOfStates )
        if self.numberOfStates <= self.requiredRoot: raise QCModelError ( "Root number incompatible with the number of configurations." )

    def ModifyElectronicState ( self, target ):
        """Modify the electronic state."""
        super ( QCModelMNDOCI, self ).ModifyElectronicState ( target )
        self.CheckVersusSystem ( target )

    def SolveCPHFEquations ( self, target ):
        """Solve the CPHF equations."""
        scratch = target.scratch
        if scratch.doGradients:
            node   = scratch.ci
            solver = node.Get ( "cphfSolver", None )
            if solver is None:
                solver = CPHFSolver.FromTarget ( target )
                node.cphfSolver = solver
            solver.SetUp ( )
            solver.Solve ( )
            # . Save the Z-matrix and a temporary density matrix needed by the gradients.
            dTotalZ = solver.work1
            zMatrix = solver.zMatrix
            node.dTotal.CopyTo ( dTotalZ )
            dTotalZ.Add ( zMatrix )
            node.dTotalZ = dTotalZ
            node.zMatrix = zMatrix
            node.warnings.extend ( solver.warnings )
            report = scratch.qcEnergyReport
            report["CPHF Degenerate Redundant"] = solver.numberDegenerateRedundant
            report["CPHF Non-Redundant"       ] = solver.numberNonRedundant
            report["CPHF Redundant"           ] = solver.numberRedundant
            report["CPHF Converged"           ] = solver.report["Is Converged"]
            report["CPHF Iterations"          ] = solver.report["Iterations"  ]
            # . For QC/MM multipole models adjust the multipoles for the QC/MM gradients.
            # . This is a fudge and the multipoles can no longer be relied upon.
            multipoleOrder = scratch.Get ( "qcmmMultipoleOrder", None )
            multipoles     = scratch.Get ( "qcmmMultipoles"    , None )
            if ( multipoleOrder is not None ) and ( multipoles is not None ):
                self.multipoleEvaluator.FockMultipoles ( target, multipoleOrder, multipoles, density = dTotalZ )

    def StateCharacters ( self, target, rotation, mapping, includeCoreOrbitals = False, selection = None ):
        """Determine the characters of CI states under a rotation."""
        state = target.qcState
        if selection is None: indices = range ( self.numberOfStates )
        else:                 indices = selection
        if includeCoreOrbitals: start = 0
        else:                   start = state.coreOrbitals
        stop                  = state.coreOrbitals + self.activeOrbitals
        orbitals              = target.scratch.orbitalsP.orbitals
        n                     = len ( state.configurations )
        o                     = stop - start
        characters            = Array.WithExtent ( len ( indices ) )
        outState              = Array.WithExtent ( n )
        orbitalTransformation = Array.WithShape ( [ o, o ] )
        rotatedOrbitals       = Array.WithShape ( [ orbitals.shape[0], o ] )
        stateTransformation   = Array.WithShape ( [ n, n ] )
        rotations             = state.orbitalBases.RotationMatrices ( rotation )
        for ( i, s ) in enumerate ( range ( start, stop ) ):
            self.RotateOrbital ( rotations, mapping, state.orbitalBases.centerFunctionPointers, orbitals[:,s], rotatedOrbitals[:,i] )
        orbitalTransformation.MatrixMultiply ( orbitals[:,start:stop], rotatedOrbitals, xTranspose = True )
        state.configurations.Characters ( orbitalTransformation, stateTransformation, coreOrbitals = state.coreOrbitals, includeCoreOrbitals = includeCoreOrbitals )
        for ( i, s ) in enumerate ( indices ):
            inState = target.scratch.ci.ciVectors[s,:]
            stateTransformation.VectorMultiply ( inState, outState )
            characters[i] = inState.Dot ( outState )
        return characters

    def TransitionDipoles ( self, target, root = None ):
        """Calculate the transition dipoles for the required root."""
        # . Units are atomic units.
        self.TransitionDipolesMake ( self, target )
        return self.TransitionDipolesExtract ( self, target, root )

    def TransitionDipolesExtract ( self, target, root ):
        """Extract the transition dipoles for root."""
        tds   = None
        node  = target.scratch.ci
        tdXYZ = node.Get ( "transitionDipoles", None )
        if tdXYZ is not None:
            s   = self.numberOfStates
            tds = Array.WithExtent ( s )
            tds.Set ( 0.0 )
            if ( root is None ) or ( root < 0 ) or ( root >= s ): root = node.Get ( "ciRoot", 0 )
            for m in tdXYZ:
                for i in s:
                    v       = m[i,root]
                    tds[i] += ( v * v )
        return tds

    def TransitionDipolesMake ( self, target ):
        """Make the transition dipoles."""
        # . Dense only for the moment.
        # . In this case the transition dipole matrix between basis functions is diagonal.
        state     = target.qcState
        a         = state.configurations.numberOfActiveOrbitals
        nAtoms    = len ( state.atomicNumbers  )
        n         = len ( state.configurations )
        s         = self.numberOfStates
        indices   = state.orbitalBases.centerFunctionPointers
        node      = target.scratch.ci
        ciVectors = node.ciVectors
        xyz       = target.scratch.qcCoordinates3AU
        tdXYZ     = node.Get ( "transitionDipoles", None )
        if tdXYZ is None:
            tdXYZ = [ Array.WithExtent ( s, storageType = StorageType.Symmetric ) for c in range ( 3 ) ]
            node.transitionDipoles = tdXYZ
        work     = Array.WithExtent ( indices[-1], storageType = StorageType.Symmetric )
        tdMatrix = Array.WithExtent ( n          , storageType = StorageType.Symmetric )
        tdMOs    = Array.WithExtent ( a          , storageType = StorageType.Symmetric )
        for ( c, m ) in enumerate ( tdXYZ ):
            work.Set ( 0.0 )
            for i in range ( nAtoms ):
                r = xyz[i,c]
                for b in range ( indices[i], indices[i+1] ): work[b,b] = r
            work.Transform ( activeMOs, tdMOs ) # . Active MOs stored as columns.
            state.configurations.TransitionDipoles ( tdMOs, tdMatrix )
            tdMatrix.Transform ( ciVectors, m, useTranpose = True ) # . CI vectors stored as rows.

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
