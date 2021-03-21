"""MNDO integral evaluator."""

from  pCore         import logFile              , \
                           LogFileActive        , \
                           RawObjectConstructor
from  pScientific   import Units
from .QCDefinitions import BasisRepresentation
from .QCModelError  import QCModelError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class MNDOIntegralEvaluator:

    def CoreCoreEnergy ( self, target ):
        """The core-core energy."""
        cdef Coordinates3            coordinates3
        cdef Coordinates3            gradients3
        cdef MNDOParametersContainer parameters
        cdef CReal                   energy      = 0.0
        cdef CRealArray2D           *cGradients3 = NULL
        scratch      = target.scratch
        coordinates3 = scratch.qcCoordinates3AU
        gradients3   = scratch.Get ( "qcGradients3AU", None )
        parameters   = target.qcState.mndoParameters
        if gradients3 is not None: cGradients3 = gradients3.cObject
        energy = MNDO_CoreCoreEnergy ( parameters.cObject   ,
                                       coordinates3.cObject ,
                                       cGradients3          )
        scratch.energyTerms["QC Core-Core"] = ( energy * Units.Energy_Hartrees_To_Kilojoules_Per_Mole )

    def DipoleIntegrals ( self, target, Vector3 center = None ):
        """The dipole integrals."""
        cdef IntegerArray1D          basisIndices
        cdef Coordinates3            coordinates3
        cdef MNDOParametersContainer parameters
        cdef SymmetricMatrix         dX
        cdef SymmetricMatrix         dY
        cdef SymmetricMatrix         dZ
        cdef CInteger                n
        cdef CRealArray1D           *cCenter = NULL
        cdef CStatus                 cStatus = CStatus_OK
        cdef CSymmetricMatrix       *cDX     = NULL
        cdef CSymmetricMatrix       *cDY     = NULL
        cdef CSymmetricMatrix       *cDZ     = NULL
        scratch      = target.scratch
        coordinates3 = scratch.qcCoordinates3AU
        basisIndices = target.qcState.orbitalBases._centerFunctionPointers[self.basisRepresentation]
        parameters   = target.qcState.mndoParameters
        if center is not None: cCenter = center.cObject
        n  = len ( target.qcState.orbitalBases )
        dX = SymmetricMatrix.WithExtent ( n )
        dY = SymmetricMatrix.WithExtent ( n )
        dZ = SymmetricMatrix.WithExtent ( n )
        MNDO_DipoleIntegrals ( parameters.cObject   ,
                               basisIndices.cObject ,
                               coordinates3.cObject ,
                               cCenter              ,
                               dX.cObject           ,
                               dY.cObject           ,
                               dZ.cObject           ,
                               &cStatus             )
        if cStatus != CStatus_OK: raise QCModelError ( "Error calculating dipole integrals." )
        return ( dX, dY, dZ )

    def ElectronNuclearTEIGradients ( self, target ):
        """The electron-nuclear and TEI gradients."""
        cdef IntegerArray1D          basisIndices
        cdef Coordinates3            coordinates3
        cdef Coordinates3            gradients3
        cdef MNDOParametersContainer parameters
        cdef SymmetricMatrix         dSpin
        cdef SymmetricMatrix         dTotal
        cdef CSymmetricMatrix       *cDSpin = NULL
        scratch = target.scratch
        if scratch.doGradients:
            basisIndices = target.qcState.orbitalBases._centerFunctionPointers[self.basisRepresentation]
            coordinates3 = scratch.qcCoordinates3AU
            dTotal       = scratch.onePDMP.density
            gradients3   = scratch.qcGradients3AU
            parameters   = target.qcState.mndoParameters
            if hasattr ( scratch, "onePDMQ" ):
                dSpin  = scratch.onePDMQ.density
                cDSpin = dSpin.cObject
            MNDO_ElectronNuclearTEIGradients ( parameters.cObject   ,
                                               basisIndices.cObject ,
                                               coordinates3.cObject ,
                                               dTotal.cObject       ,
                                               cDSpin               ,
                                               gradients3.cObject   )

    def ElectronNuclearTEIGradientsCI ( self, target ):
        """The CI electron-nuclear and TEI gradients."""
        cdef IntegerArray1D          basisIndices
        cdef Coordinates3            coordinates3
        cdef Coordinates3            gradients3
        cdef DoubleSymmetricMatrix   twoPDM
        cdef MNDOParametersContainer parameters
        cdef RealArray2D             orbitals
        cdef SymmetricMatrix         dCore
        cdef SymmetricMatrix         dHF
        cdef SymmetricMatrix         dTotalZ
        cdef SymmetricMatrix         onePDM
        cdef SymmetricMatrix         zMatrix
        scratch = target.scratch
        if scratch.doGradients:
            basisIndices = target.qcState.orbitalBases._centerFunctionPointers[self.basisRepresentation]
            coordinates3 = scratch.qcCoordinates3AU
            dHF          = scratch.onePDMP.density
            orbitals     = scratch.orbitalsP.orbitals
            gradients3   = scratch.qcGradients3AU
            parameters   = target.qcState.mndoParameters
            node         = scratch.ci
            dCore        = node.dCore
            dTotalZ      = node.dTotalZ
            onePDM       = node.onePDM
            twoPDM       = node.twoPDM
            zMatrix      = node.zMatrix
            MNDO_ElectronNuclearTEIGradientsCI ( target.qcModel.activeOrbitals       ,
                                                 target.qcState.coreOrbitals         ,
                                                 len ( target.qcState.orbitalBases ) ,
                                                 parameters.cObject                  , 
                                                 basisIndices.cObject                , 
                                                 coordinates3.cObject                , 
                                                 twoPDM.cObject                      , 
                                                 orbitals.cObject                    , 
                                                 dCore.cObject                       , 
                                                 dHF.cObject                         , 
                                                 dTotalZ.cObject                     , 
                                                 onePDM.cObject                      , 
                                                 zMatrix.cObject                     , 
                                                 gradients3.cObject                  ) 

    def ElectronNuclearTEIIntegrals ( self, target ):
        """The electron-nuclear and TEI integrals."""
        cdef BlockStorage            twoElectronIntegrals
        cdef IntegerArray1D          basisIndices
        cdef Coordinates3            coordinates3
        cdef MNDOParametersContainer parameters
        cdef SymmetricMatrix         oneElectronMatrix
        cdef CBlockStorage          *cTwoElectronIntegrals = NULL
        scratch           = target.scratch
        basisIndices      = target.qcState.orbitalBases._centerFunctionPointers[self.basisRepresentation]
        coordinates3      = scratch.qcCoordinates3AU
        oneElectronMatrix = scratch.oneElectronMatrix
        parameters        = target.qcState.mndoParameters
        MNDO_ElectronNuclearTEIIntegrals ( parameters.cObject        ,
                                           basisIndices.cObject      ,
                                           coordinates3.cObject      ,
                                           oneElectronMatrix.cObject ,
                                           &cTwoElectronIntegrals    )
        if cTwoElectronIntegrals == NULL: raise QCModelError ( "Error calculating MNDO two-electron integrals." )
        twoElectronIntegrals         = BlockStorage.Raw ( )
        twoElectronIntegrals.cObject = cTwoElectronIntegrals
        twoElectronIntegrals.isOwner = True
        scratch.twoElectronIntegrals = twoElectronIntegrals

    def ResonanceGradients ( self, target, doCI = False ):
        """The resonance gradients."""
        cdef IntegerArray1D          basisIndices
        cdef Coordinates3            coordinates3
        cdef Coordinates3            gradients3
        cdef GaussianBasisContainer  bases
        cdef MNDOParametersContainer parameters
        cdef SymmetricMatrix         dTotal
        scratch = target.scratch
        if scratch.doGradients:
            bases        = target.qcState.orbitalBases
            basisIndices = bases._centerFunctionPointers[self.basisRepresentation]
            coordinates3 = scratch.qcCoordinates3AU
            gradients3   = scratch.qcGradients3AU
            parameters   = target.qcState.mndoParameters
            if doCI: dTotal = scratch.ci.dTotalZ
            else:    dTotal = scratch.onePDMP.density
            MNDO_ResonanceGradients ( parameters.cObject   ,
                                      bases.cObject        ,
                                      basisIndices.cObject ,
                                      coordinates3.cObject ,
                                      dTotal.cObject       ,
                                      gradients3.cObject   )

    def ResonanceIntegrals ( self, target ):
        """The resonance integrals."""
        cdef IntegerArray1D          basisIndices
        cdef Coordinates3            coordinates3
        cdef GaussianBasisContainer  bases
        cdef MNDOParametersContainer parameters
        cdef SymmetricMatrix         oneElectronMatrix
        bases             = target.qcState.orbitalBases
        basisIndices      = bases._centerFunctionPointers[self.basisRepresentation]
        coordinates3      = target.scratch.qcCoordinates3AU
        parameters        = target.qcState.mndoParameters
        oneElectronMatrix = target.scratch.oneElectronMatrix
        MNDO_ResonanceIntegrals ( parameters.cObject        ,
                                  bases.cObject             ,
                                  basisIndices.cObject      ,
                                  coordinates3.cObject      ,
                                  oneElectronMatrix.cObject )

    @property
    def basisRepresentation ( self ): return BasisRepresentation.Actual

