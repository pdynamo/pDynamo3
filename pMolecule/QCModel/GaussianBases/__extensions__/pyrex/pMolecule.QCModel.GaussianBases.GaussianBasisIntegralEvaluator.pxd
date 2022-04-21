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
from pMolecule.QCModel.GaussianBases.BlockStorage           cimport BlockStorage            , \
                                                                    CBlockStorage
from pMolecule.QCModel.GaussianBases.GaussianBasis          cimport CGaussianBasisOperator
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
from pScientific.Geometry3.Vector3                          cimport Vector3

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

cdef extern from "GaussianBasisContainerIntegrals_f1Xg2.h":

    cdef void GaussianBasisContainerIntegrals_f1Xg2i    ( CGaussianBasisContainer *self              ,
                                                          CGaussianBasisContainer *other             ,
                                                          CRealArray2D            *coordinates3      ,
                                                          CGaussianBasisOperator   operator          ,
                                                          CBlockStorage           *fitIntegrals      ,
                                                          CStatus                 *status            )
    cdef void GaussianBasisContainerIntegrals_f1Xg2R1   ( CGaussianBasisContainer *self              ,
                                                          CGaussianBasisContainer *other             ,
                                                          CRealArray2D            *coordinates3      ,
                                                          CSymmetricMatrix        *density           ,
                                                          CRealArray1D            *xVector           ,
                                                          CGaussianBasisOperator   operator          ,
                                                          CRealArray2D            *gradients3        ,
                                                          CStatus                 *status            )

cdef extern from "GaussianBasisContainerIntegrals_f1Op1.h":

    cdef void GaussianBasisContainerIntegrals_f1Op1i    ( CGaussianBasisContainer *self              ,
                                                          CRealArray2D            *coordinates3      ,
                                                          CRealArray2D            *rG                ,
                                                          CRealArray2D            *values            ,
                                                          CStatus                 *status            )

cdef extern from "GaussianBasisContainerIntegrals_f1X.h":

    cdef void GaussianBasisContainerIntegrals_f1Di      ( CGaussianBasisContainer *self              ,
                                                          CRealArray2D            *coordinates3      ,
                                                          CRealArray1D            *center            ,
                                                          CRealArray1D            *dipoleX           ,
                                                          CRealArray1D            *dipoleY           ,
                                                          CRealArray1D            *dipoleZ           ,
                                                          CStatus                 *status            )
    cdef void GaussianBasisContainerIntegrals_f1Oi      ( CGaussianBasisContainer *self              ,
                                                          CRealArray1D            *overlap           ,
                                                          CStatus                 *status            )
    cdef void GaussianBasisContainerIntegrals_f1Qi      ( CGaussianBasisContainer *self              ,
                                                          CRealArray2D            *coordinates3      ,
                                                          CRealArray1D            *center            ,
                                                          CRealArray1D            *qXX               ,
                                                          CRealArray1D            *qYY               ,
                                                          CRealArray1D            *qZZ               ,
                                                          CRealArray1D            *qXY               ,
                                                          CRealArray1D            *qXZ               ,
                                                          CRealArray1D            *qYZ               ,
                                                          CStatus                 *status            )

cdef extern from "GaussianBasisContainerIntegrals_f1Xg1.h":

    cdef void GaussianBasisContainerIntegrals_f1Af1i    ( CGaussianBasisContainer *self              ,
                                                          CRealArray2D            *coordinates3      ,
                                                          CSymmetricMatrix        *integrals         ,
                                                          CStatus                 *status            )
    cdef void GaussianBasisContainerIntegrals_f1Cf1i    ( CGaussianBasisContainer *self              ,
                                                          CRealArray2D            *coordinates3      ,
                                                          CSymmetricMatrix        *integrals         ,
                                                          CStatus                 *status            )
    cdef void GaussianBasisContainerIntegrals_f1Df1i    ( CGaussianBasisContainer *self              ,
                                                          CRealArray2D            *coordinates3      ,
                                                          CRealArray1D            *center            ,
                                                          CSymmetricMatrix        *dipoleX           ,
                                                          CSymmetricMatrix        *dipoleY           ,
                                                          CSymmetricMatrix        *dipoleZ           ,
                                                          CStatus                 *status            )
    cdef void GaussianBasisContainerIntegrals_f1KOf1i   ( CGaussianBasisContainer *self              ,
                                                          CRealArray2D            *coordinates3      ,
                                                          CSymmetricMatrix        *kinetic           ,
                                                          CSymmetricMatrix        *overlap           ,
                                                          CStatus                 *status            )
    cdef void GaussianBasisContainerIntegrals_f1KOf1R1  ( CGaussianBasisContainer *self              ,
                                                          CRealArray2D            *coordinates3      ,
                                                          CSymmetricMatrix        *kDensity          ,
                                                          CSymmetricMatrix        *oDensity          ,
                                                          CRealArray2D            *gradients3        ,
                                                          CStatus                 *status            )
    cdef void GaussianBasisContainerIntegrals_f1Of1i    ( CGaussianBasisContainer *self              ,
                                                          CRealArray2D            *coordinates3      ,
                                                          CSymmetricMatrix        *integrals         ,
                                                          CStatus                 *status            )
    cdef void GaussianBasisContainerIntegrals_f1Qf1i    ( CGaussianBasisContainer *self              ,
                                                          CRealArray2D            *coordinates3      ,
                                                          CRealArray1D            *center            ,
                                                          CSymmetricMatrix        *qXX               ,
                                                          CSymmetricMatrix        *qYY               ,
                                                          CSymmetricMatrix        *qZZ               ,
                                                          CSymmetricMatrix        *qXY               ,
                                                          CSymmetricMatrix        *qXZ               ,
                                                          CSymmetricMatrix        *qYZ               ,
                                                          CStatus                 *status            )
    cdef void GaussianBasisContainerIntegrals_f1Xf1R1   ( CGaussianBasisContainer *self              ,
                                                          CRealArray2D            *coordinates3      ,
                                                          CRealArray1D            *aVector           ,
                                                          CRealArray1D            *xVector           ,
                                                          CGaussianBasisOperator   operator          ,
                                                          CRealArray2D            *gradients3        ,
                                                          CStatus                 *status            )

cdef extern from "GaussianBasisContainerIntegrals_f2Cf2.h":

    cdef void GaussianBasisContainerIntegrals_f2Cf2i    ( CGaussianBasisContainer *self              ,
                                                          CRealArray2D            *coordinates3      ,
                                                          CBlockStorage           *teis              ,
                                                          CStatus                 *status            )
    cdef void GaussianBasisContainerIntegrals_f2Cf2R1   ( CGaussianBasisContainer *self              ,
                                                          CRealArray2D            *coordinates3      ,
                                                          CSymmetricMatrix        *dTotal            ,
                                                          CSymmetricMatrix        *dSpin             ,
                                                          CBoolean                 doCoulomb         ,
                                                          CReal                    exchangeScaling   ,
                                                          CRealArray2D            *gradients3        ,
                                                          CStatus                 *status            )

cdef extern from "GaussianBasisContainerIntegrals_f2Cp1.h":

    cdef void GaussianBasisContainerIntegrals_f2Cm1R1   ( CGaussianBasisContainer *self              ,
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
    cdef void GaussianBasisContainerIntegrals_f2Cm1V    ( CGaussianBasisContainer *self              ,
                                                          CRealArray1D            *charges           ,
                                                          CRealArray1D            *widthsE           ,
                                                          CRealArray1D            *widthsN           ,
                                                          CRealArray2D            *coordinates3      ,
                                                          CRealArray2D            *coordinates3G     ,
                                                          CSelection              *selectionG        ,
                                                          CSymmetricMatrix        *oneElectronMatrix ,
                                                          CStatus                 *status            )
    cdef void GaussianBasisContainerIntegrals_f2Cp1V    ( CGaussianBasisContainer *self              ,
                                                          CRealArray1D            *widthsE           ,
                                                          CRealArray1D            *widthsN           ,
                                                          CRealArray2D            *coordinates3      ,
                                                          CRealArray2D            *coordinates3G     ,
                                                          CSelection              *selectionG        ,
                                                          CSymmetricMatrix        *density           ,
                                                          CRealArray1D            *potentials        ,
                                                          CStatus                 *status            )

cdef extern from "GaussianBasisContainerIntegrals_p1Cq1.h":

    cdef CReal GaussianBasisContainer_m1Cn1ER1          ( CRealArray1D            *chargesI          ,
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
    cdef void GaussianBasisContainer_m1Cp1V             ( CRealArray1D            *chargesJ          ,
                                                          CRealArray2D            *coordinates3I     ,
                                                          CRealArray2D            *coordinates3J     ,
                                                          CSelection              *selectionI        ,
                                                          CSelection              *selectionJ        ,
                                                          CRealArray1D            *widthsEI          ,
                                                          CRealArray1D            *widthsEJ          ,
                                                          CRealArray1D            *widthsNI          ,
                                                          CRealArray1D            *widthsNJ          ,
                                                          CRealArray1D            *potentialsI       )          

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class GaussianBasisIntegralEvaluator:
    pass
