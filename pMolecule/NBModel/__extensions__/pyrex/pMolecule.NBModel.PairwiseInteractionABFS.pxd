from pCore.CPrimitiveTypes                           cimport CBoolean                    , \
                                                             CFalse                      , \
                                                             CInteger                    , \
                                                             CReal                       , \
                                                             CTrue
from pCore.PairList                                  cimport CPairList                   , \
                                                             PairList
from pCore.Status                                    cimport CStatus                     , \
                                                             CStatus_OK
from pMolecule.MMModel.LJParameterContainer          cimport CLJParameterContainer       , \
                                                             LJParameterContainer
from pMolecule.NBModel.ImagePairListContainer        cimport CImagePairListContainer     , \
                                                             ImagePairListContainer
from pMolecule.NBModel.PairwiseInteraction           cimport PairwiseInteraction
from pScientific.Arrays.IntegerArray1D               cimport CIntegerArray1D              , \
                                                             IntegerArray1D
from pScientific.Arrays.RealArray1D                  cimport CRealArray1D                 , \
                                                             RealArray1D
from pScientific.Arrays.RealArray2D                  cimport CRealArray2D
from pScientific.Arrays.SymmetricMatrix              cimport CSymmetricMatrix             , \
                                                             SymmetricMatrix
from pScientific.Geometry3.Coordinates3              cimport Coordinates3
from pScientific.Symmetry.SymmetryParameters         cimport CSymmetryParameters         , \
                                                             SymmetryParameters
from pScientific.Symmetry.SymmetryParameterGradients cimport CSymmetryParameterGradients , \
                                                             SymmetryParameterGradients

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "PairwiseInteractionABFS.h":

    ctypedef struct CPairwiseInteractionABFS "PairwiseInteractionABFS":
        CReal dampingCutOff
        CReal innerCutOff
        CReal outerCutOff

    cdef  CPairwiseInteractionABFS *PairwiseInteractionABFS_Allocate            ( CStatus                     *status                     )
    cdef  CPairwiseInteractionABFS *PairwiseInteractionABFS_Clone               ( CPairwiseInteractionABFS    *self                       ,
                                                                                  CStatus                     *status                     )
    cdef  void                      PairwiseInteractionABFS_Deallocate          ( CPairwiseInteractionABFS   **self                       )
    cdef  void                      PairwiseInteractionABFS_Interactions        ( CPairwiseInteractionABFS    *self                       ,
                                                                                  CRealArray1D                *r                          ,
                                                                                  CRealArray1D                *electrostatic              ,
                                                                                  CRealArray1D                *lennardJonesA              ,
                                                                                  CRealArray1D                *lennardJonesB              )
    cdef  void                      PairwiseInteractionABFS_MMMMEnergy          ( CPairwiseInteractionABFS    *self                       ,
                                                                                  CRealArray1D                *chargesI                   ,
                                                                                  CRealArray1D                *chargesJ                   ,
                                                                                  CIntegerArray1D             *ljTypesI                   ,
                                                                                  CIntegerArray1D             *ljTypesJ                   ,
                                                                                  CLJParameterContainer       *ljParameters               ,
                                                                                  CReal                        electrostaticScale         ,
                                                                                  CReal                        lennardJonesScale          ,
                                                                                  CRealArray2D                *coordinates3I              ,
                                                                                  CRealArray2D                *coordinates3J              ,
                                                                                  CPairList                   *pairList                   ,
                                                                                  CReal                       *eElectrostatic             ,
                                                                                  CReal                       *eLennardJones              ,
                                                                                  CRealArray2D                *gradients3I                ,
                                                                                  CRealArray2D                *gradients3J                ,
                                                                                  CStatus                     *status                     )
    cdef  void                      PairwiseInteractionABFS_MMMMEnergyImage     ( CPairwiseInteractionABFS    *self                       ,
                                                                                  CRealArray1D                *charges                    ,
                                                                                  CIntegerArray1D             *ljTypes                    ,
                                                                                  CLJParameterContainer       *ljParameters               ,
                                                                                  CReal                        electrostaticScale         ,
                                                                                  CRealArray2D                *coordinates3               ,
                                                                                  CSymmetryParameters         *symmetryParameters         ,
                                                                                  CImagePairListContainer     *imagePairLists             ,
                                                                                  CReal                       *eElectrostatic             ,
                                                                                  CReal                       *eLennardJones              ,
                                                                                  CRealArray2D                *gradients3                 ,
                                                                                  CSymmetryParameterGradients *symmetryParameterGradients ,
                                                                                  CStatus                     *status                     )
    cdef  void                      PairwiseInteractionABFS_MMMMEnergyMI        ( CPairwiseInteractionABFS    *self                       ,
                                                                                  CRealArray1D                *chargesI                   ,
                                                                                  CRealArray1D                *chargesJ                   ,
                                                                                  CIntegerArray1D             *ljTypesI                   ,
                                                                                  CIntegerArray1D             *ljTypesJ                   ,
                                                                                  CLJParameterContainer       *ljParameters               ,
                                                                                  CReal                        electrostaticScale         ,
                                                                                  CReal                        lennardJonesScale          ,
                                                                                  CRealArray2D                *coordinates3I              ,
                                                                                  CRealArray2D                *coordinates3J              ,
                                                                                  CSymmetryParameters         *symmetryParameters         ,
                                                                                  CPairList                   *pairList                   ,
                                                                                  CReal                       *eElectrostatic             ,
                                                                                  CReal                       *eLennardJones              ,
                                                                                  CRealArray2D                *gradients3I                ,
                                                                                  CRealArray2D                *gradients3J                ,
                                                                                  CSymmetryParameterGradients *symmetryParameterGradients ,
                                                                                  CStatus                     *status                     )
    cdef  void                      PairwiseInteractionABFS_QCMMGradients       ( CPairwiseInteractionABFS    *self                       ,
                                                                                  CRealArray1D                *chargesQ                   ,
                                                                                  CRealArray1D                *chargesM                   ,
                                                                                  CReal                        electrostaticScale         ,
                                                                                  CRealArray2D                *coordinates3Q              ,
                                                                                  CRealArray2D                *coordinates3M              ,
                                                                                  CPairList                   *pairList                   ,
                                                                                  CRealArray2D                *gradients3Q                ,
                                                                                  CRealArray2D                *gradients3M                ,
                                                                                  CStatus                     *status                     )
    cdef  void                      PairwiseInteractionABFS_QCMMGradientsImage  ( CPairwiseInteractionABFS    *self                       ,
                                                                                  CRealArray1D                *chargesA                   ,
                                                                                  CRealArray1D                *chargesB                   ,
                                                                                  CReal                        electrostaticScale         ,
                                                                                  CRealArray2D                *coordinates3A              ,
                                                                                  CRealArray2D                *coordinates3B              ,
                                                                                  CSymmetryParameters         *symmetryParameters         ,
                                                                                  CImagePairListContainer     *imagePairLists             ,
                                                                                  CRealArray2D                *gradients3A                ,
                                                                                  CRealArray2D                *gradients3B                ,
                                                                                  CSymmetryParameterGradients *symmetryParameterGradients ,
                                                                                  CStatus                     *status                     )
    cdef  void                      PairwiseInteractionABFS_QCMMGradientsMI     ( CPairwiseInteractionABFS    *self                       ,
                                                                                  CRealArray1D                *chargesQ                   ,
                                                                                  CRealArray1D                *chargesM                   ,
                                                                                  CReal                        electrostaticScale         ,
                                                                                  CRealArray2D                *coordinates3Q              ,
                                                                                  CRealArray2D                *coordinates3M              ,
                                                                                  CSymmetryParameters         *symmetryParameters         ,
                                                                                  CPairList                   *pairList                   ,
                                                                                  CRealArray2D                *gradients3Q                ,
                                                                                  CRealArray2D                *gradients3M                ,
                                                                                  CSymmetryParameterGradients *symmetryParameterGradients ,
                                                                                  CStatus                     *status                     )
    cdef  void                      PairwiseInteractionABFS_QCMMPotentials      ( CPairwiseInteractionABFS    *self                       ,
                                                                                  CRealArray1D                *chargesM                   ,
                                                                                  CReal                        electrostaticScale         ,
                                                                                  CRealArray2D                *coordinates3Q              ,
                                                                                  CRealArray2D                *coordinates3M              ,
                                                                                  CPairList                   *pairList                   ,
                                                                                  CRealArray1D                *potentials                 ,
                                                                                  CStatus                     *status                     )
    cdef  void                      PairwiseInteractionABFS_QCMMPotentialsImage ( CPairwiseInteractionABFS    *self                       ,
                                                                                  CRealArray1D                *charges                    ,
                                                                                  CReal                        electrostaticScale         ,
                                                                                  CRealArray2D                *coordinates3A              ,
                                                                                  CRealArray2D                *coordinates3B              ,
                                                                                  CSymmetryParameters         *symmetryParameters         ,
                                                                                  CImagePairListContainer     *imagePairLists             ,
                                                                                  CRealArray1D                *potentials                 ,
                                                                                  CStatus                     *status                     )
    cdef  void                      PairwiseInteractionABFS_QCMMPotentialsMI    ( CPairwiseInteractionABFS    *self                       ,
                                                                                  CRealArray1D                *chargesM                   ,
                                                                                  CReal                        electrostaticScale         ,
                                                                                  CRealArray2D                *coordinates3Q              ,
                                                                                  CRealArray2D                *coordinates3M              ,
                                                                                  CSymmetryParameters         *symmetryParameters         ,
                                                                                  CPairList                   *pairList                   ,
                                                                                  CRealArray1D                *potentials                 ,
                                                                                  CStatus                     *status                     )
#   cdef  void                      PairwiseInteractionABFS_QCQCGradients       ( CPairwiseInteractionABFS    *self                       ,
#                                                                                 CRealArray1D                *charges                    ,
#                                                                                 CReal                        electrostaticScale         ,
#                                                                                 CRealArray2D                *coordinates3I              ,
#                                                                                 CRealArray2D                *coordinates3J              ,
#                                                                                 CPairList                   *pairList                   ,
#                                                                                 CRealArray2D                *gradients3I                ,
#                                                                                 CRealArray2D                *gradients3J                ,
#                                                                                 CStatus                     *status                     )
    cdef  void                      PairwiseInteractionABFS_QCQCGradientsImage  ( CPairwiseInteractionABFS    *self                       ,
                                                                                  CRealArray1D                *charges                    ,
                                                                                  CReal                        electrostaticScale         ,
                                                                                  CRealArray2D                *coordinates3               ,
                                                                                  CSymmetryParameters         *symmetryParameters         ,
                                                                                  CImagePairListContainer     *imagePairLists             ,
                                                                                  CRealArray2D                *gradients3                 ,
                                                                                  CSymmetryParameterGradients *symmetryParameterGradients ,
                                                                                  CStatus                     *status                     )
#   cdef  void                      PairwiseInteractionABFS_QCQCPotentials      ( CPairwiseInteractionABFS    *self                       ,
#                                                                                 CReal                        electrostaticScale         ,
#                                                                                 CRealArray2D                *coordinates3I              ,
#                                                                                 CRealArray2D                *coordinates3J              ,
#                                                                                 CPairList                   *pairList                   ,
#                                                                                 CSymmetricMatrix            *potentials                 ,
#                                                                                 CStatus                     *status                     )
    cdef  void                      PairwiseInteractionABFS_QCQCPotentialsImage ( CPairwiseInteractionABFS    *self                       ,
                                                                                  CReal                        electrostaticScale         ,
                                                                                  CRealArray2D                *coordinates3               ,
                                                                                  CSymmetryParameters         *symmetryParameters         ,
                                                                                  CImagePairListContainer     *imagePairLists             ,
                                                                                  CSymmetricMatrix            *potentials                 ,
                                                                                  CStatus                     *status                     )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class PairwiseInteractionABFS ( PairwiseInteraction ):

    cdef CPairwiseInteractionABFS   *cObject
    cdef public object               isOwner
