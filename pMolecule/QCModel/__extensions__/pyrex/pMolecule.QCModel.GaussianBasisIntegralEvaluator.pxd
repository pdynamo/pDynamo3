from pCore.CPrimitiveTypes                     cimport CBoolean                , \
                                                       CCardinal               , \
                                                       CFalse                  , \
                                                       CInteger                , \
                                                       CReal                   , \
                                                       CTrue
from pCore.Selection                           cimport CSelection              , \
                                                       Selection
from pCore.Status                              cimport CStatus                 , \
                                                       CStatus_OK
from pMolecule.QCModel.BlockStorage            cimport BlockStorage            , \
                                                       CBlockStorage
from pMolecule.QCModel.GaussianBasisContainer  cimport CGaussianBasisContainer , \
                                                       GaussianBasisContainer
from pScientific.Arrays.IntegerArray1D         cimport CIntegerArray1D         , \
                                                       IntegerArray1D
from pScientific.Arrays.RealArray1D            cimport CRealArray1D            , \
                                                       RealArray1D
from pScientific.Arrays.RealArray2D            cimport CRealArray2D            , \
                                                       RealArray2D
from pScientific.Arrays.SymmetricMatrix        cimport CSymmetricMatrix        , \
                                                       SymmetricMatrix
from pScientific.Geometry3.Coordinates3        cimport Coordinates3
from pScientific.Geometry3.Vector3             cimport Vector3

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "GaussianBasisContainerIntegrals_b0e0n2.h":

    cdef CReal GaussianBasisContainer_NuclearNuclearEnergy              ( CRealArray1D            *chargesI          ,
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
    cdef void GaussianBasisContainer_NuclearNuclearPotentials           ( CRealArray1D            *chargesJ          ,
                                                                          CRealArray2D            *coordinates3I     ,
                                                                          CRealArray2D            *coordinates3J     ,
                                                                          CSelection              *selectionI        ,
                                                                          CSelection              *selectionJ        ,
                                                                          CRealArray1D            *widthsEI          ,
                                                                          CRealArray1D            *widthsEJ          ,
                                                                          CRealArray1D            *widthsNI          ,
                                                                          CRealArray1D            *widthsNJ          ,
                                                                          CRealArray1D            *potentialsI       )          

cdef extern from "GaussianBasisContainerIntegrals_b1e0n1.h":

    cdef void GaussianBasisContainerIntegrals_Grid                      ( CGaussianBasisContainer *self              ,
                                                                          CIntegerArray1D         *basisIndices      ,
                                                                          CRealArray2D            *coordinates3      ,
                                                                          CRealArray2D            *rG                ,
                                                                          CRealArray2D            *values            ,
                                                                          CStatus                 *status            )

cdef extern from "GaussianBasisContainerIntegrals_b1e1n0.h":

    cdef void GaussianBasisContainerIntegrals_SelfOverlap               ( CGaussianBasisContainer *self              ,
                                                                          CIntegerArray1D         *basisIndices      ,
                                                                          CRealArray1D            *selfOverlap       ,
                                                                          CStatus                 *status            )

cdef extern from "GaussianBasisContainerIntegrals_b2e1n0.h":

    cdef void GaussianBasisContainerIntegrals_2Coulomb                  ( CGaussianBasisContainer *self              ,
                                                                          CIntegerArray1D         *basisIndices      ,
                                                                          CRealArray2D            *coordinates3      ,
                                                                          CSymmetricMatrix        *integrals         ,
                                                                          CStatus                 *status            )
    cdef void GaussianBasisContainerIntegrals_2CoulombD                 ( CGaussianBasisContainer *self              ,
                                                                          CIntegerArray1D         *basisIndices      ,
                                                                          CRealArray2D            *coordinates3      ,
                                                                          CRealArray1D            *fPotential        ,
                                                                          CRealArray1D            *wVector           ,
                                                                          CRealArray2D            *gradients3        ,
                                                                          CStatus                 *status            )
    cdef void GaussianBasisContainerIntegrals_2Overlap                  ( CGaussianBasisContainer *self              ,
                                                                          CIntegerArray1D         *basisIndices      ,
                                                                          CRealArray2D            *coordinates3      ,
                                                                          CSymmetricMatrix        *integrals         ,
                                                                          CStatus                 *status            )
    cdef void GaussianBasisContainerIntegrals_Dipole                    ( CGaussianBasisContainer *self              ,
                                                                          CIntegerArray1D         *basisIndices      ,
                                                                          CRealArray2D            *coordinates3      ,
                                                                          CRealArray1D            *center            ,
                                                                          CSymmetricMatrix        *dipoleX           ,
                                                                          CSymmetricMatrix        *dipoleY           ,
                                                                          CSymmetricMatrix        *dipoleZ           ,
                                                                          CStatus                 *status            )
    cdef void GaussianBasisContainerIntegrals_Kinetic2Overlap           ( CGaussianBasisContainer *self              ,
                                                                          CIntegerArray1D         *basisIndices      ,
                                                                          CRealArray2D            *coordinates3      ,
                                                                          CSymmetricMatrix        *kinetic           ,
                                                                          CSymmetricMatrix        *overlap           ,
                                                                          CStatus                 *status            )
    cdef void GaussianBasisContainerIntegrals_Kinetic2OverlapD          ( CGaussianBasisContainer *self              ,
                                                                          CIntegerArray1D         *basisIndices      ,
                                                                          CRealArray2D            *coordinates3      ,
                                                                          CSymmetricMatrix        *kDensity          ,
                                                                          CSymmetricMatrix        *oDensity          ,
                                                                          CRealArray2D            *gradients3        ,
                                                                          CStatus                 *status            )

cdef extern from "GaussianBasisContainerIntegrals_b2e1n1.h":

    cdef void GaussianBasisContainerIntegrals_ElectronNuclear           ( CGaussianBasisContainer *self              ,
                                                                          CIntegerArray1D         *basisIndices      ,
                                                                          CRealArray1D            *charges           ,
                                                                          CRealArray1D            *widthsE           ,
                                                                          CRealArray1D            *widthsN           ,
                                                                          CRealArray2D            *coordinates3      ,
                                                                          CRealArray2D            *coordinates3G     ,
                                                                          CSelection              *selectionG        ,
                                                                          CSymmetricMatrix        *oneElectronMatrix ,
                                                                          CStatus                 *status            )
    cdef void GaussianBasisContainerIntegrals_ElectronNuclearD          ( CGaussianBasisContainer *self              ,
                                                                          CIntegerArray1D         *basisIndices      ,
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
    cdef void GaussianBasisContainerIntegrals_ElectronNuclearPotentials ( CGaussianBasisContainer *self              ,
                                                                          CIntegerArray1D         *basisIndices      ,
                                                                          CRealArray1D            *widthsE           ,
                                                                          CRealArray1D            *widthsN           ,
                                                                          CRealArray2D            *coordinates3      ,
                                                                          CRealArray2D            *coordinates3G     ,
                                                                          CSelection              *selectionG        ,
                                                                          CSymmetricMatrix        *density           ,
                                                                          CRealArray1D            *potentials        ,
                                                                          CStatus                 *status            )

cdef extern from "GaussianBasisContainerIntegrals_b3e2n0.h":

    cdef void GaussianBasisContainerIntegrals_ElectronFit               ( CGaussianBasisContainer *self              ,
                                                                          CIntegerArray1D         *selfIndices       ,
                                                                          CGaussianBasisContainer *other             ,
                                                                          CIntegerArray1D         *otherIndices      ,
                                                                          CRealArray2D            *coordinates3      ,
                                                                          CBlockStorage           *fitIntegrals      ,
                                                                          CStatus                 *status            )
    cdef void GaussianBasisContainerIntegrals_ElectronFitD              ( CGaussianBasisContainer *self              ,
                                                                          CIntegerArray1D         *selfIndices       ,
                                                                          CGaussianBasisContainer *other             ,
                                                                          CIntegerArray1D         *otherIndices      ,
                                                                          CRealArray2D            *coordinates3      ,
                                                                          CSymmetricMatrix        *sDensity          ,
                                                                          CRealArray1D            *oDensity          ,
                                                                          CRealArray2D            *gradients3        ,
                                                                          CStatus                 *status            )

cdef extern from "GaussianBasisContainerIntegrals_b4e2n0.h":

    cdef void GaussianBasisContainerIntegrals_TEIs                      ( CGaussianBasisContainer *self              ,
                                                                          CIntegerArray1D         *basisIndices      ,
                                                                          CRealArray2D            *coordinates3      ,
                                                                          CBlockStorage           *teis              ,
                                                                          CStatus                 *status            )
    cdef void GaussianBasisContainerIntegrals_TEIsD                     ( CGaussianBasisContainer *self              ,
                                                                          CIntegerArray1D         *basisIndices      ,
                                                                          CRealArray2D            *coordinates3      ,
                                                                          CSymmetricMatrix        *dTotal            ,
                                                                          CSymmetricMatrix        *dSpin             ,
                                                                          CBoolean                 doCoulomb         ,
                                                                          CReal                    exchangeScaling   ,
                                                                          CRealArray2D            *gradients3        ,
                                                                          CStatus                 *status            )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class GaussianBasisIntegralEvaluator:
    pass
