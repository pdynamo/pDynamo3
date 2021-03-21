"""Classes for pDynamo's in-built QC models."""

from   pCore                         import Clone                 , \
                                            logFile               , \
                                            LogFileActive
from   pScientific                   import Units
from   pScientific.Arrays            import Array                 , \
                                            StorageType
from   pScientific.Geometry3         import Vector3
from  .DIISSCFConverger              import DIISSCFConverger
from  .ElectronicState               import SpinType
from  .FockConstruction              import FockConstruction_MakeFromTEIs
from  .OrthogonalizingTransformation import OrthogonalizingTransformation_Make
from  .QCDefinitions                 import BasisRepresentation   , \
                                            ChargeModel           , \
                                            FockClosurePriority   , \
                                            OrthogonalizationType
from  .QCModel                       import QCModel               , \
                                            QCModelState
from  .QCModelError                  import QCModelError
from  .QCOnePDM                      import QCOnePDM
from  .QCOrbitals                    import QCOrbitals
from ..EnergyModel                   import EnergyClosurePriority

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCModelBaseState ( QCModelState ):
    """A QC base model state."""

    _attributable = dict ( QCModelState._attributable )
    _unpicklable  = set  ( QCModelState._unpicklable  )
    _attributable.update ( { "alphaCharge"    : 0.0  ,
                             "betaCharge"     : 0.0  ,
                             "energyBaseLine" : 0.0  ,
                             "fockClosures"   : None ,
                             "fockModels"     : dict ,
                             "nuclearCharges" : None ,
                             "orbitalBases"   : None } )
    _unpicklable.add ( "fockClosures" )

    def _UpdateFockClosures ( self ):
        """Update the list of fock closures."""
        closures = []
        for value in self.fockModels.values ( ):
            for ( priority, closure ) in value.FockClosures ( self.target ):
                closures.append ( ( priority, closure ) )
        self.fockClosures = [ closure for ( priority, closure ) in sorted ( closures, key = lambda x: x[0] ) ]

    def AddFockModel ( self, key, value ):
        """Add a Fock model."""
        if value is None: self.fockModels.pop ( key, None )
        else:             self.fockModels[key] = value
        self._UpdateFockClosures ( )

    def SummaryItems ( self ):
        """Summary items."""
        items = super ( QCModelBaseState, self ).SummaryItems ( )
        items.extend ( [ ( "Alpha Charge"            , "{:.3f}".format ( self.alphaCharge          ) ) ,
                         ( "Beta Charge"             , "{:.3f}".format ( self.betaCharge           ) ) ,
                         ( "Energy Base Line"        , "{:.3f}".format ( self.energyBaseLine       ) ) ,
                         ( "Orbital Basis Functions" , "{:d}".format   ( len ( self.orbitalBases ) ) ) ] )
        return items

    @property
    def spinDensity ( self ):
        """The spin density for property calculations."""
        onePDM = self.target.scratch.Get ( "onePDMQ", None )
        if onePDM is None: return None
        else:              return onePDM.density

    @property
    def totalDensity ( self ):
        """The total density for property calculations."""
        onePDM = self.target.scratch.Get ( "onePDMP", None )
        if onePDM is None: return None
        else:              return onePDM.density

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCModelBase ( QCModel ):
    """Base class for pDynamo's in-built QC models."""

    _attributable = dict ( QCModel._attributable )
    _chargeModels = {}
    _classLabel   = "Base QC Model"
    _stateObject  = QCModelBaseState
    _summarizable = dict ( QCModel._summarizable )
    _attributable.update ( { "exchangeScaling"       : 1.0                                    ,
                             "orthogonalizationType" : OrthogonalizationType.Symmetric        ,
                             "converger"             : None                                   ,
                             "integralEvaluator"     : None                                   ,
                             "multipoleEvaluator"    : None                                   } )
    _summarizable.update ( { "converger"             :   "SCF Converger"                      ,
                             "exchangeScaling"       : ( "Exchange Scaling"      , "{:.3f}" ) ,
                             "orthogonalizationType" :   "Orthogonalization Type"             } )

    def AtomicCharges ( self, target, chargeModel = None ):
        """Atomic charges."""
        evaluator = self.GetChargeModelEvaluator ( chargeModel )
        state     = getattr ( target, self.__class__._stateName )
        charges   = Array.WithExtent ( len ( state.atomicNumbers ) )
        evaluator.FockMultipoles ( target, 0, charges, density = state.totalDensity )
        return charges

    def AtomicSpins ( self, target, chargeModel = None ):
        """Atomic spins."""
        evaluator = self.GetChargeModelEvaluator ( chargeModel )
        state     = getattr ( target, self.__class__._stateName )
        spins     = Array.WithExtent ( len ( state.atomicNumbers ) )
        if state.spinDensity is None:
            spins.Set ( 0.0 )
        else:
            evaluator.FockMultipoles ( target, 0, spins, density = state.spinDensity, withNuclei = False )
            spins.Scale ( -1.0 )
        return spins

    def BondOrders ( self, target, chargeModel = None ):
        """Bond orders."""
        evaluator    = self.GetChargeModelEvaluator ( chargeModel )
        state        = getattr ( target, self.__class__._stateName )
        n            = len ( state.atomicNumbers )
        bOs          = Array.WithExtent ( n, storageType = StorageType.Symmetric ) ; bOs.Set ( 0.0 )
        freeValence  = Array.WithExtent ( n ) ; freeValence.Set  ( 0.0 )
        totalValence = Array.WithExtent ( n ) ; totalValence.Set ( 0.0 )
        for ( p, ( density, valence ) ) in enumerate ( [ ( state.totalDensity, totalValence ) ,
                                                         ( state.spinDensity , freeValence  ) ] ):
            if density is not None:
                evaluator.BondOrders ( target, density, bOs )
                for a in range ( n ):
                    sum = 0.0
                    for b in range ( n ):
                        if b != a: sum += bOs[a,b]
                    valence[a] = sum
                if p == 1:
                    freeValence.Scale ( -1.0 )
                    freeValence.Add   ( totalValence )
        return ( bOs, freeValence, totalValence )

    def BuildModel ( self, target, qcSelection = None ):
        """Build the model."""
        state = super ( QCModelBase, self ).BuildModel ( target, qcSelection = qcSelection )
        self.GetParameters ( target )
        ( state.alphaCharge, state.betaCharge ) = target.electronicState.Verify ( sum ( state.nuclearCharges ) )
        # . Ensure the basis and evaluator representations are consistent.
        state.orbitalBases.basisRepresentation = self.integralEvaluator.basisRepresentation
        state.AddFockModel ( self.__class__._classLabel, self                   )
        state.AddFockModel ( "Electronic State"        , target.electronicState )
        if self.converger is None: self.converger = DIISSCFConverger ( )
        return state

    def DipoleMoment ( self, target, center = None, dipole = None ):
        """The dipole moment in Debyes."""
        cX = cY = cZ = 0.0
        dX = dY = dZ = 0.0
        if center is not None:
            cX = center[0] ; cY = center[1] ; cZ = center[2]
        state        = getattr ( target, self.__class__._stateName )
        scratch      = target.scratch
        coordinates3 = scratch.qcCoordinates3AU
        dTotal       = state.totalDensity
        for ( i, q ) in enumerate ( state.nuclearCharges ):
            dX += q * ( coordinates3[i,0] - cX )
            dY += q * ( coordinates3[i,1] - cY )
            dZ += q * ( coordinates3[i,2] - cZ )
        ( iX, iY, iZ ) = self.integralEvaluator.DipoleIntegrals ( target, center = center )
        dX -= dTotal.TraceOfProduct ( iX )
        dY -= dTotal.TraceOfProduct ( iY )
        dZ -= dTotal.TraceOfProduct ( iZ )
        if dipole is None: dipole = Vector3.Null ( )
        dipole[0] += ( dX * Units.Dipole_Atomic_Units_To_Debyes )
        dipole[1] += ( dY * Units.Dipole_Atomic_Units_To_Debyes )
        dipole[2] += ( dZ * Units.Dipole_Atomic_Units_To_Debyes )
        return dipole

    def Energy ( self, target ):
        """Calculate the QC energy."""
        # . It is assumed that the last converger call to the Fock closures correspond to the best energy.
        # . More sophisticated handling would be preferable here.
        report   = self.converger.Iterate ( target, log = target.scratch.log )
        error    = report.get ( "Error", None )
        qcReport = target.scratch.qcEnergyReport
        qcReport["SCF Converged" ] = report["Converged" ]
        qcReport["SCF Iterations"] = report["Iterations"]
        if error is not None: raise QCModelError ( "Converger error: {:s}".format ( error ) )

    def EnergyFinalize ( self, target ):
        """Energy finalization."""
        super ( QCModelBase, self ).EnergyFinalize ( target )
        log = target.scratch.log
        if log is not None:
            report = target.scratch.qcEnergyReport
            if hasattr ( target.scratch, "onePDMQ" ): # . Maybe should always do these even if not printing?
                ( S2, Sz ) = target.scratch.onePDMP.SpinExpectationValues ( target.scratch.onePDMQ, overlap = target.scratch.Get ( "overlapMatrix", None ) )
                report["<Sz>"] = ( Sz, "{:.3f}" )
                report["<S2>"] = ( S2, "{:.3f}" )
            items = []
            for ( key, value ) in report.items ( ):
                if value is not None:
                    if   isinstance ( value, tuple ): valueString = value[1].format ( value[0] ) # . 2-tuple of ( value, format ).
                    elif isinstance ( value, float ): valueString = "{:.6f}".format ( value    ) # . Suitable for energy in atomic units.
                    elif isinstance ( value, str   ): valueString = value
                    else:                             valueString = str ( value )
                    items.append ( ( key, valueString ) )
            log.SummaryOfItems ( items, title = "QC Energy Report" )

    def EnergyInitialize ( self, target ):
        """Energy initialization."""
        super ( QCModelBase, self ).EnergyInitialize ( target )
        state   = getattr ( target, self.__class__._stateName )
        if not hasattr ( state, "fockClosures" ) or ( state.fockClosures is None ): state._UpdateFockClosures ( )
        self.SetUpDensities ( target )
        n       = len ( state.orbitalBases )
        scratch = target.scratch
        oem     = scratch.Get ( "oneElectronMatrix", None )
        if ( oem is None ) or ( oem.rows != n ):
            oem = Array.WithExtent ( n, storageType = StorageType.Symmetric )
            scratch.oneElectronMatrix = oem
        oem.Set ( 0.0 )
        scratch.qcEnergyReport = { }

    def FockOne ( self, target ):
        """The one-electron contribution to the Fock matrices."""
        scratch = target.scratch
        dTotal  = scratch.onePDMP.density
        fTotal  = scratch.onePDMP.fock
        oem     = scratch.oneElectronMatrix
        eOE     = dTotal.TraceOfProduct ( oem )
        fTotal.Add ( oem )
        scratch.qcEnergyReport["One-Electron Energy"] = eOE
        eOE    += target.qcState.energyBaseLine                                   # . Include baseline.
        eQCE    = scratch.qcEnergyReport.pop ( "QC Electronic Accumulator", 0.0 ) # . Remove accumulator.
        scratch.energyTerms["QC Electronic"] = ( ( eOE + eQCE ) * Units.Energy_Hartrees_To_Kilojoules_Per_Mole )
        return eOE

    def FockClosures ( self, target ):
        """Fock closures."""
        def a ( ):
            return self.FockOne ( target )
        def b ( ):
            return self.FockTwo ( target )
        return [ ( FockClosurePriority.VeryLow , a ) ,
                 ( FockClosurePriority.VeryHigh, b ) ]

    # . Accumulator and Fock matrices initialized here.
    def FockTwo ( self, target ):
        """The two-electron contribution to the Fock matrices."""
        scratch = target.scratch
        doSpin  = hasattr ( scratch, "onePDMQ" )
        dTotal  = scratch.onePDMP.density
        fTotal  = scratch.onePDMP.fock
        if doSpin:
            dSpin = scratch.onePDMQ.density
            fSpin = scratch.onePDMQ.fock
        else:
            dSpin = None
            fSpin = None
        eTE = FockConstruction_MakeFromTEIs ( dTotal                       ,
                                              dSpin                        ,
                                              scratch.twoElectronIntegrals ,
                                              self.exchangeScaling         ,
                                              fTotal                       ,
                                              fSpin                        )
        scratch.qcEnergyReport["QC Electronic Accumulator"] = eTE
        scratch.qcEnergyReport["Two-Electron Energy"      ] = eTE
        return eTE

    def GetChargeModelEvaluator ( self, chargeModel ):
        """Get the charge model evaluator."""
        if chargeModel is None:
            return self.multipoleEvaluator
        elif chargeModel in self.__class__._chargeModels:
            evaluator = self.__class__._chargeModels[chargeModel]
            if isinstance ( self.multipoleEvaluator, evaluator ):
                return self.multipoleEvaluator
            else:
                return evaluator ( )
        else: raise QCModelError ( "Invalid charge model for this QC model." )

    # . X - orthogonalizer, Y = S * X - inverse orthogonalizer such that Y^T * X = I.
    # . For orthogonal density and orbitals from non-orthogonal quantities need Po = Y^T * Pn * Y and A = Y^T * C.
    # . The inverses are Pn = X * Po * X^T and C = X * A.
    def GetNonOrthogonalDensity ( self, target, spinType = None ):
        """Get the density in the non-orthogonal basis."""
        state = getattr ( target, self.__class__._stateName )
        if   spinType in ( None, SpinType.Total ): density = state.totalDensity
        elif spinType ==         SpinType.Spin   : density = state.spinDensity
        else: raise QCModelError ( "Unknown density spin type: {:s}.".format ( repr ( spinType ) ) )
        return density

    def GetNonOrthogonalOrbitals ( self, target, indices = None, spinType = None ):
        """Get the orbitals in the non-orthogonal basis."""
        if   spinType in ( None, SpinType.Alpha, SpinType.Total ): orbitals = target.scratch.orbitalsP.orbitals
        elif spinType ==                         SpinType.Beta   : orbitals = target.scratch.orbitalsQ.orbitals
        else: raise QCModelError ( "Unknown orbital spin type: {:s}.".format ( repr ( spinType ) ) )
        if ( indices is not None ) and ( len ( indices ) > 0 ):
            indices = sorted ( set ( indices ) )
            if ( indices[0] >= 0 ) and ( indices[-1] < orbitals.shape[1] ):
                if len ( indices ) < orbitals.shape[1]:
                    old      = orbitals
                    orbitals = Array.WithShape ( [ old.shape[0], len ( indices ) ] )
                    for ( i, s ) in enumerate ( indices ):
                        old[:,s].CopyTo ( orbitals[:,i] )
            else: raise QCModelError ( "Orbital indices out of range." )
        return orbitals

    # . This is not general as it assumes the overlap matrix is in the W-basis.
    def GetOrthogonalizer ( self, target, orthogonalizationType = None, doInverse = True, doLoewdin = True ):
        """Get the orthogonalizing transformation and, optionally, its inverse."""
        state = getattr ( target, self.__class__._stateName )
        if orthogonalizationType is None: orthogonalizationType = self.orthogonalizationType
        # . Transform the overlap to the proper basis: Sa = U^T * Sw * U.
        w2a = state.orbitalBases.w2a
        Sw  = target.scratch.overlapMatrix
        Sa  = Array.WithExtent ( w2a.shape[1], storageType = StorageType.Symmetric )
        Sw.Transform ( w2a, Sa )
        # . Make the transformation: Xa^T * Sa * Xa = I.
        ( X, eigenValues, eigenVectors ) = OrthogonalizingTransformation_Make ( Sa ,
                                                                                doCanonical   = ( orthogonalizationType is OrthogonalizationType.Canonical ) ,
                                                                                preserveInput = False )
        # . Save the eigenvalues and vectors of Sa.
        target.scratch.overlapEigenValues  = eigenValues
        target.scratch.overlapEigenVectors = eigenVectors
        if ( X.columns != X.rows ) and hasattr ( target.scratch, "qcEnergyReport" ):
            target.scratch.qcEnergyReport["Basis Linear Dependence"] = abs ( X.columns - X.rows )
        # . Transform back to the working basis and save: Xw = U * Xa => Xw^T * Sw * Xw = I.
        orthogonalizer = Array.WithShape ( [ w2a.shape[0], X.shape[1] ] )
        orthogonalizer.MatrixMultiply ( w2a, X )
        target.scratch.orthogonalizer = orthogonalizer
        # . Make the inverse directly from the orthogonalizer: Yw = Sw * Xw => Yw^T * Xw = I.
        if doInverse:
            inverseOrthogonalizer = Array.WithShape ( orthogonalizer.shape )
            Sw.PostMatrixMultiply ( orthogonalizer, inverseOrthogonalizer )
            target.scratch.inverseOrthogonalizer = inverseOrthogonalizer
        # . Make the Loewdin transformation.
        if doLoewdin:
            if ( orthogonalizationType is OrthogonalizationType.Symmetric ) and ( target.scratch.Get ( "inverseOrthogonalizer", None ) is not None ):
                loewdinT = inverseOrthogonalizer
            else:
                n           = eigenValues.extent
                sHalf       = Array.WithExtent ( n, storageType = StorageType.Symmetric )
                squareRoots = Clone ( eigenValues )
                squareRoots.Power ( 0.5 )
                sHalf.MakeFromEigenSystem ( n, squareRoots, eigenVectors )
                #ArrayPrint2D ( sHalf, itemFormat = "{:.6f}", title = "Sa^(1/2)" )
                loewdinT    = Array.WithShape ( orthogonalizer.shape )
                sHalf.PreMatrixMultiply ( state.orbitalBases.a2w, loewdinT )
            target.scratch.loewdinT = loewdinT

    def GetParameters ( self, target ):
        """Get the parameters for the model."""
        pass

    # . Input grid point coordinates for the following methods should be in atomic units.
    def GridPointDensities ( self, target, gridPoints, spinType = None ):
        """Densities at grid points."""
        density = self.GetNonOrthogonalDensity ( target, spinType = spinType )
        bValues = Array.WithShape ( [ density.shape[0], gridPoints.shape[0] ] )
        values  = Array.WithExtent ( gridPoints.shape[0] )
        self.gridPointEvaluator.GridValues ( target, gridPoints, bValues )
        density.DiagonalOfTransform ( bValues, values )
        return values

    def GridPointOrbitals ( self, target, gridPoints, orbitalIndices, spinType = None ):
        """Orbital values at grid points."""
        values = None
        if ( orbitalIndices is not None ) and ( len ( orbitalIndices ) > 0 ):
            orbitals = self.GetNonOrthogonalOrbitals ( target, indices = orbitalIndices, spinType = spinType )
            bValues  = Array.WithShape ( [ orbitals.shape[0], gridPoints.shape[0] ] )
            values   = Array.WithShape ( [ gridPoints.shape[0], len ( orbitalIndices ) ] )
            self.gridPointEvaluator.GridValues ( target, gridPoints, bValues )
            values.MatrixMultiply ( bValues, orbitals, xTranspose = True )
        else: raise QCModelError ( "Missing orbital index specification." )
        return values

    def GridPointPotentials ( self, target, gridPoints, spinType = None ):
        """Electrostatic potentials at grid points."""
        density = self.GetNonOrthogonalDensity ( target, spinType = spinType )
        values  = Array.WithExtent ( gridPoints.shape[0] )
        values.Set ( 0.0 )
        self.gridPointEvaluator.ElectronPotentials ( target, density, gridPoints, values )
        if spinType in ( None, SpinType.Total ):
            self.gridPointEvaluator.NuclearPotentials ( target, gridPoints, values )
        return values

    def ModifyElectronicState ( self, target ):
        """Modify the electronic state."""
        super ( QCModelBase, self ).ModifyElectronicState ( target )
        state = getattr ( target, self.__class__._stateName )
        ( state.alphaCharge, state.betaCharge ) = target.electronicState.Verify ( sum ( state.nuclearCharges ) )
        state.AddFockModel ( "Electronic State", target.electronicState )

    def OrbitalCharacters ( self, target, rotation, mapping, selection = None, useOrbitalsP = True ):
        """Determine the characters of a selected set of orbitals under a rotation."""
        state          = getattr ( target, self.__class__._stateName )
        orbitalIndices = state.orbitalBases.centerFunctionPointers
        if useOrbitalsP: orbitals = target.scratch.orbitalsP.orbitals
        else:            orbitals = target.scratch.orbitalsQ.orbitals
        if selection is None: indices = range ( orbitals.shape[-1] )
        else:                 indices = selection
        characters = Array.WithExtent ( len ( indices ) )
        orbitalOut = Array.WithExtent ( len ( state.orbitalBases ) )
        overlap    = target.scratch.Get ( "overlapMatrix", None )
        if overlap is None:
            sOrbitals = orbitals
        else:
            # . Clearly this could be more efficient by storing it for use with other rotations ...
            sOrbitals = Array.WithShape ( orbitals.shape )
            overlap.PostMatrixMultiply ( orbitals, sOrbitals ) 
        rotations = state.orbitalBases.RotationMatrices ( rotation )
        for ( i, s ) in enumerate ( indices ):
            self.RotateOrbital ( rotations, mapping, orbitalIndices, orbitals[:,s], orbitalOut )
            characters[i] = orbitalOut.Dot ( sOrbitals[:,s] )
        return characters

    @staticmethod
    def RotateOrbital ( rotations, mapping, orbitalIndices, orbitalIn, orbitalOut ):
        """Rotate an orbital."""
        for ( a, R ) in enumerate ( rotations ):
            b      = mapping[a]
            aStart = orbitalIndices[a  ]
            aStop  = orbitalIndices[a+1]
            bStart = orbitalIndices[b  ]
            bStop  = orbitalIndices[b+1]
            R.VectorMultiply ( orbitalIn[aStart:aStop], orbitalOut[bStart:bStop] )

    def SetUpDensities ( self, target ):
        """Set up the one-particle densities."""
        state   = getattr ( target, self.__class__._stateName )
        extent  = len ( state.orbitalBases )
        scratch = target.scratch
        items   = [ ( "onePDMP", SpinType.Total, state.alphaCharge + state.betaCharge ) ]
        if not target.electronicState.isSpinRestricted:
            items.append ( ( "onePDMQ", SpinType.Spin, state.alphaCharge - state.betaCharge ) )
        for ( name, spinType, charge ) in items:
            item = scratch.Get ( name, None )
            if ( item                is None     ) or \
               ( item.numberOrbitals != extent   ) or \
               ( item.spinType       != spinType ) or \
               ( item.totalCharge    != charge   ):
                scratch.Set ( name, QCOnePDM.FromDiagonalGuess ( extent, spinType, charge ) )

    def SetUpOrbitals ( self, target ):
        """Set up the orbital sets."""
        state                = getattr ( target, self.__class__._stateName )
        bases                = state.orbitalBases
        numberBasisFunctions = len ( bases ) # . A-basis or W-basis depending on integral evaluator.
        numberOrbitals       = bases._numberOfFunctions[BasisRepresentation.Actual] # . Always A-basis.
        handlers             = target.electronicState.OccupancyHandlers ( state.alphaCharge, state.betaCharge )
        scratch              = target.scratch
        for ( name, handler ) in zip ( ( "orbitalsP", "orbitalsQ" ), handlers ):
            item = scratch.Get ( name, None )
            if ( item                is None           ) or \
               ( item.numberOrbitals != numberOrbitals ) or \
               ( not isinstance ( item.occupancyHandler, handler.__class__ ) ): # . More reliable checks here ...
                scratch.Set ( name, QCOrbitals.WithExtents ( numberBasisFunctions, numberOrbitals, handler ) )

    @property
    def gridPointEvaluator ( self ):
        """The grid point evaluator for property calculations."""
        return self.integralEvaluator

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
