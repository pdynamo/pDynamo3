from pCore.CPrimitiveTypes                        cimport CBoolean                 , \
                                                          CFalse                   , \
                                                          CInteger                 , \
                                                          CReal                    , \
                                                          CTrue
from pCore.PairList                               cimport CPairList                , \
                                                          PairList
from pCore.Status                                 cimport CStatus                  , \
                                                          CStatus_OK
from pMolecule.QCModel.GaussianBases.BlockStorage cimport BlockStorage             , \
                                                          CBlockStorage
from pMolecule.QCModel.MNDOParameters             cimport CMNDOParameters          , \
                                                          MNDOParameters
from pMolecule.QCModel.MNDOParametersContainer    cimport CMNDOParametersContainer , \
                                                          MNDOParametersContainer
from pScientific.Arrays.IntegerArray1D            cimport IntegerArray1D           , \
                                                          CIntegerArray1D      
from pScientific.Arrays.RealArray1D               cimport CRealArray1D             , \
                                                          RealArray1D
from pScientific.Arrays.RealArray2D               cimport CRealArray2D             , \
                                                          RealArray2D
from pScientific.Arrays.SymmetricMatrix           cimport CSymmetricMatrix         , \
                                                          SymmetricMatrix
from pScientific.Geometry3.Coordinates3           cimport Coordinates3
from pScientific.Splines.CubicSplineContainer     cimport CCubicSplineContainer    , \
                                                          CubicSplineContainer

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "MNDOIntegralsMM.h":


    cdef void MNDOIntegralsMM_Values ( CMNDOParameters *qData ,
                                       CReal            R     ,
                                       CInteger         index ,
                                       CReal           *f     ,
                                       CReal           *g     )

cdef extern from "MNDOQCMM.h":

    cdef void  MNDO_QCMMGradients ( CIntegerArray1D          *atomIndices         ,
                                    CSymmetricMatrix         *dTotal              ,
                                    CBlockStorage            *integrals           ,
                                    CRealArray2D             *qcGradients3        ,
                                    CRealArray2D             *mmGradients3        ,
                                    CStatus                  *status              )
    cdef CReal MNDO_QCMMIntegrals ( CMNDOParametersContainer *parameters          ,
                                    CIntegerArray1D          *basisIndices        ,
                                    CCubicSplineContainer    *splines             ,
                                    CReal                     cutOff              ,
                                    CReal                     eScale              ,
                                    CRealArray2D             *qcCoordinates3      ,
                                    CRealArray2D             *mmCoordinates3      ,
                                    CRealArray1D             *mmCharges           ,
                                    CPairList                *pairList            ,
                                    CSymmetricMatrix         *oneElectronMatrix   ,
                                    CRealArray2D             *qcGradients3        ,
                                    CRealArray2D             *mmGradients3        ,
                                    CBlockStorage           **derivativeIntegrals ,
                                    CStatus                  *status              )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class MNDOQCMMEvaluator:
    pass
