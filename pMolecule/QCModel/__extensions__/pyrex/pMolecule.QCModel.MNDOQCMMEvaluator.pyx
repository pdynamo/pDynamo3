"""MNDO QC/MM density evaluator."""

from .QCDefinitions import BasisRepresentation
from .QCModelError  import QCModelError

# . Everything in standard units except one-electron matrix.

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . QC/MM term labels and the number of terms as a function of the number of orbitals.
_MNDOQCMMTermLabels = ( "Core Signed"   ,
                        "Core Unsigned" ,
                        "S/S"           ,
                        "Pz/S"          ,
                        "Pz/Pz"         ,
                        "Px/Px"         ,
                        "Dz2/S"         ,
                        "Dz2/Pz"        ,
                        "Dz2/Dz2"       ,
                        "Dxz/Px"        ,
                        "Dxz/Dxz"       ,
                        "Dx2y2/Dx2y2"   )
                        
_MNDOQCMMTerms      = { 0 : 2 , 1 : 3 , 4 : 6 , 9 : 12 }

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class MNDOQCMMEvaluator:

    def Gradients ( self                                  ,
                    IntegerArray1D  atomIndices  not None ,
                    SymmetricMatrix dTotal       not None ,
                    BlockStorage    integrals    not None ,
                    RealArray2D     qcGradients3 not None ,
                    RealArray2D     mmGradients3 not None ):
        """Gradients."""
        cdef CStatus cStatus = CStatus_OK
        MNDO_QCMMGradients ( atomIndices.cObject  ,
                             dTotal.cObject       ,
                             integrals.cObject    ,
                             qcGradients3.cObject ,
                             mmGradients3.cObject ,
                             &cStatus             )
        if cStatus != CStatus_OK: raise QCModelError ( "Error evaluating QC/MM gradients." )

    def Integrals ( self                                                ,
                    MNDOParametersContainer parameters         not None ,
                    IntegerArray1D          basisIndices       not None ,
                    CubicSplineContainer    splines                     ,
                    CReal                   electrostaticScale          ,
                    RealArray2D             qcCoordinates3     not None ,
                    RealArray2D             mmCoordinates3     not None ,
                    RealArray1D             mmCharges          not None ,
                    PairList                pairList           not None ,
                    SymmetricMatrix         oneElectronMatrix  not None ,
                    Coordinates3            qcGradients3                ,
                    Coordinates3            mmGradients3                ,
                    CReal                   cutOff = -1.0               ): # . -1.0 is infinite.
        """Integrals."""
        cdef BlockStorage           integrals
        cdef CBlockStorage         *cIntegrals    = NULL
        cdef CBlockStorage        **pIntegrals    = NULL
        cdef CCubicSplineContainer *cSplines      = NULL
        cdef CReal                  energy
        cdef CRealArray2D          *cMMGradients3 = NULL
        cdef CRealArray2D          *cQCGradients3 = NULL
        cdef CStatus                cStatus       = CStatus_OK
        if ( mmGradients3 is not None ) and ( qcGradients3 is not None ):
            cMMGradients3 = mmGradients3.cObject
            cQCGradients3 = qcGradients3.cObject
            pIntegrals    = &cIntegrals
        if splines is not None:
            cSplines = splines.cObject
        energy = MNDO_QCMMIntegrals ( parameters.cObject        ,
                                      basisIndices.cObject      ,
                                      cSplines                  ,
                                      cutOff                    ,
                                      electrostaticScale        ,
                                      qcCoordinates3.cObject    ,
                                      mmCoordinates3.cObject    ,
                                      mmCharges.cObject         ,
                                      pairList.cObject          ,
                                      oneElectronMatrix.cObject ,
                                      cQCGradients3             ,
                                      cMMGradients3             ,
                                      pIntegrals                ,
                                      &cStatus                  )
        if cStatus != CStatus_OK: raise QCModelError ( "Error evaluating QC/MM integrals." )
        if cIntegrals != NULL:
            integrals         = BlockStorage.Raw ( )
            integrals.cObject = cIntegrals
            integrals.isOwner = True
        else:
            integrals = None
        return ( energy, integrals )

    def IntegralValue ( self                               ,
                        MNDOParameters parameters not None ,
                        CReal          R                   ,
                        CInteger       index               ,
                                       derivative = False  ):
        """Return values of the core-charge and local frame integrals or derivatives (in atomic units)."""
        cdef CReal *pF = NULL
        cdef CReal *pG = NULL
        cdef CReal  value
        n = _MNDOQCMMTerms[parameters.numberOfOrbitals]
        if ( index < 0 ) or ( index >= n ): raise QCModelError ( "Term index out-of-range: {:d}/{:d}.".format ( index, n ) )
        if derivative: pG = &value
        else:          pF = &value
        MNDOIntegralsMM_Values ( parameters.cObject, R, index, pF, pG )
        return value

    @property
    def basisRepresentation ( self ): return BasisRepresentation.Actual
