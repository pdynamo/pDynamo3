from pCore.CPrimitiveTypes                           cimport CBoolean                    , \
                                                             CFalse                      , \
                                                             CInteger                    , \
                                                             CReal                       , \
                                                             CTrue
from pCore.PairList                                  cimport CPairList                   , \
                                                             PairList
from pCore.Status                                    cimport CStatus                     , \
                                                             CStatus_OK
from pMolecule.NBModel.ImagePairListContainer        cimport CImagePairListContainer     , \
                                                             ImagePairListContainer
from pMolecule.QCModel.BlockStorage                  cimport BlockStorage                , \
                                                             CBlockStorage
from pMolecule.QCModel.BlockStorageContainer         cimport BlockStorageContainer       , \
                                                             CBlockStorageContainer
from pMolecule.QCModel.MNDOParameters                cimport CMNDOParameters             , \
                                                             MNDOParameters
from pMolecule.QCModel.MNDOParametersContainer       cimport CMNDOParametersContainer    , \
                                                             MNDOParametersContainer
from pMolecule.QCModel.MNDOQCMMEvaluator             cimport MNDOQCMMEvaluator
from pScientific.Arrays.IntegerArray1D               cimport IntegerArray1D              , \
                                                             CIntegerArray1D      
from pScientific.Arrays.RealArray1D                  cimport CRealArray1D                , \
                                                             RealArray1D
from pScientific.Arrays.RealArray2D                  cimport CRealArray2D                , \
                                                             RealArray2D
from pScientific.Arrays.SymmetricMatrix              cimport CSymmetricMatrix            , \
                                                             SymmetricMatrix
from pScientific.Geometry3.Coordinates3              cimport Coordinates3
from pScientific.Splines.CubicSplineContainer        cimport CCubicSplineContainer       , \
                                                             CubicSplineContainer
from pScientific.Symmetry.SymmetryParameters         cimport CSymmetryParameters         , \
                                                             SymmetryParameters
from pScientific.Symmetry.SymmetryParameterGradients cimport CSymmetryParameterGradients , \
                                                             SymmetryParameterGradients

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "MNDOQCMMImage.h":

    void  MNDOQCMMImage_QCMMGradientsImage  ( CIntegerArray1D             *atomIndices                ,
                                              CSymmetricMatrix            *dTotal                     ,
                                              CRealArray2D                *coordinates3A              ,
                                              CRealArray2D                *coordinates3B              ,
                                              CSymmetryParameters         *symmetryParameters         ,
                                              CImagePairListContainer     *imagePairLists             ,
                                              CBlockStorageContainer      *integralContainer          ,
                                              CRealArray2D                *gradients3A                ,
                                              CRealArray2D                *gradients3B                ,
                                              CSymmetryParameterGradients *symmetryParameterGradients ,
                                              CStatus                     *status                     )
    CReal MNDOQCMMImage_QCMMPotentialsImage ( CMNDOParametersContainer    *parameters                 ,
                                              CIntegerArray1D             *basisIndices               ,
                                              CCubicSplineContainer       *splines                    ,
                                              CReal                        cutOff                     ,
                                              CReal                        eScale                     ,
                                              CRealArray2D                *qcCoordinates3             ,
                                              CRealArray2D                *mmCoordinates3             ,
                                              CSymmetryParameters         *symmetryParameters         ,
                                              CRealArray1D                *mmCharges                  ,
                                              CImagePairListContainer     *imagePairLists             ,
                                              CSymmetricMatrix            *oneElectronMatrix          ,
                                              CRealArray2D                *qcGradients3               ,
                                              CRealArray2D                *mmGradients3               ,
                                              CSymmetryParameterGradients *symmetryParameterGradients ,
                                              CBlockStorageContainer      *integralContainer          ,
                                              CStatus                     *status                     )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class MNDOQCMMImageEvaluator ( MNDOQCMMEvaluator ):
    pass
