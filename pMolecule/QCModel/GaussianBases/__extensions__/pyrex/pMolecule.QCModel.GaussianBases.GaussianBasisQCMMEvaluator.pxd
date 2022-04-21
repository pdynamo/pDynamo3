from pCore.CPrimitiveTypes                                  cimport CBoolean                , \
                                                                    CCardinal               , \
                                                                    CFalse                  , \
                                                                    CInteger                , \
                                                                    CReal                   , \
                                                                    CTrue
from pCore.Selection                                        cimport CSelection              , \
                                                                    Selection
from pCore.Status                                           cimport CStatus                 , \
                                                                    CStatus_OK
from pMolecule.QCModel.GaussianBases.GaussianBasisContainer cimport CGaussianBasisContainer , \
                                                                    GaussianBasisContainer
from pScientific.Arrays.IntegerArray1D                      cimport CIntegerArray1D         , \
                                                                    IntegerArray1D
from pScientific.Arrays.RealArray1D                         cimport CRealArray1D            , \
                                                                    RealArray1D
from pScientific.Arrays.RealArray2D                         cimport CRealArray2D            , \
                                                                    RealArray2D
from pScientific.Arrays.SymmetricMatrix                     cimport CSymmetricMatrix        , \
                                                                    SymmetricMatrix
from pScientific.Geometry3.Coordinates3                     cimport Coordinates3

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "GaussianBasisContainerIntegrals_f1Cp1.h":

    cdef void GaussianBasisContainerIntegrals_f1Cm1R1 ( CGaussianBasisContainer *self              ,
                                                        CRealArray1D            *charges           ,
                                                        CRealArray1D            *widthsE           ,
                                                        CRealArray1D            *widthsN           ,
                                                        CRealArray2D            *coordinates3      ,
                                                        CRealArray2D            *coordinates3G     ,
                                                        CSelection              *selectionG        ,
                                                        CRealArray1D            *fitCoefficients   ,
                                                        CRealArray2D            *gradients3        ,
                                                        CRealArray2D            *gradients3G       ,
                                                        CStatus                 *status            )
    cdef void GaussianBasisContainerIntegrals_f1Cm1V  ( CGaussianBasisContainer *self              ,
                                                        CRealArray1D            *charges           ,
                                                        CRealArray1D            *widthsE           ,
                                                        CRealArray1D            *widthsN           ,
                                                        CRealArray2D            *coordinates3      ,
                                                        CRealArray2D            *coordinates3G     ,
                                                        CSelection              *selectionG        ,
                                                        CRealArray1D            *integrals         ,
                                                        CStatus                 *status            )
    cdef void GaussianBasisContainerIntegrals_f1Cp1V  ( CGaussianBasisContainer *self              ,
                                                        CRealArray1D            *widthsE           ,
                                                        CRealArray1D            *widthsN           ,
                                                        CRealArray2D            *coordinates3      ,
                                                        CRealArray2D            *coordinates3G     ,
                                                        CSelection              *selectionG        ,
                                                        CRealArray1D            *fitCoefficients   ,
                                                        CRealArray1D            *potentials        ,
                                                        CStatus                 *status            )

cdef extern from "GaussianBasisContainerIntegrals_f2Cp1.h":

    cdef void GaussianBasisContainerIntegrals_f2Cm1R1 ( CGaussianBasisContainer *self              ,
                                                        CRealArray1D            *charges           ,
                                                        CRealArray1D            *widthsE           ,
                                                        CRealArray1D            *widthsN           ,
                                                        CRealArray2D            *coordinates3      ,
                                                        CRealArray2D            *coordinates3G     ,
                                                        CSelection              *selectionG        ,
                                                        CSymmetricMatrix        *density           ,
                                                        CRealArray2D            *gradients3        ,
                                                        CRealArray2D            *gradients3G       ,
                                                        CStatus                 *status            )
    cdef void GaussianBasisContainerIntegrals_f2Cm1V  ( CGaussianBasisContainer *self              ,
                                                        CRealArray1D            *charges           ,
                                                        CRealArray1D            *widthsE           ,
                                                        CRealArray1D            *widthsN           ,
                                                        CRealArray2D            *coordinates3      ,
                                                        CRealArray2D            *coordinates3G     ,
                                                        CSelection              *selectionG        ,
                                                        CSymmetricMatrix        *oneElectronMatrix ,
                                                        CStatus                 *status            )
    cdef void GaussianBasisContainerIntegrals_f2Cp1V  ( CGaussianBasisContainer *self              ,
                                                        CRealArray1D            *widthsE           ,
                                                        CRealArray1D            *widthsN           ,
                                                        CRealArray2D            *coordinates3      ,
                                                        CRealArray2D            *coordinates3G     ,
                                                        CSelection              *selectionG        ,
                                                        CSymmetricMatrix        *density           ,
                                                        CRealArray1D            *potentials        ,
                                                        CStatus                 *status            )

cdef extern from "GaussianBasisContainerIntegrals_p1Cq1.h":

    cdef CReal GaussianBasisContainer_m1Cn1ER1        ( CRealArray1D            *chargesI          ,
                                                        CRealArray1D            *chargesJ          ,
                                                        CRealArray2D            *coordinates3I     ,
                                                        CRealArray2D            *coordinates3J     ,
                                                        CSelection              *selectionI        ,
                                                        CSelection              *selectionJ        ,
                                                        CRealArray1D            *widthsEI          ,
                                                        CRealArray1D            *widthsEJ          ,
                                                        CRealArray1D            *widthsNI          ,
                                                        CRealArray1D            *widthsNJ          ,
                                                        CRealArray2D            *gradients3I       ,
                                                        CRealArray2D            *gradients3J       )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class GaussianBasisQCMMEvaluator:
    pass
