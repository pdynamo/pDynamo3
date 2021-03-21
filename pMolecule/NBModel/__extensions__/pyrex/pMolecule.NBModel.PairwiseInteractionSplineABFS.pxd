from pCore.CPrimitiveTypes                           cimport CBoolean                     , \
                                                             CFalse                       , \
                                                             CInteger                     , \
                                                             CReal                        , \
                                                             CTrue
from pCore.PairList                                  cimport CPairList                    , \
                                                             PairList
from pCore.Status                                    cimport CStatus                      , \
                                                             CStatus_OK
from pMolecule.MMModel.LJParameterContainer          cimport CLJParameterContainer        , \
                                                             LJParameterContainer
from pMolecule.NBModel.ImagePairListContainer        cimport CImagePairListContainer      , \
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
from pScientific.Splines.CubicSpline                 cimport CCubicSpline                 , \
                                                             CubicSpline
from pScientific.Symmetry.SymmetryParameters         cimport CSymmetryParameters          , \
                                                             SymmetryParameters
from pScientific.Symmetry.SymmetryParameterGradients cimport CSymmetryParameterGradients  , \
                                                             SymmetryParameterGradients

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "PairwiseInteractionSpline.h":

    ctypedef struct CPairwiseInteractionSpline "PairwiseInteractionSpline":
        CReal         cutOff
        CReal         cutOff2
        CCubicSpline *electrostaticSpline
        CCubicSpline *lennardJonesASpline
        CCubicSpline *lennardJonesBSpline

    cdef  CPairwiseInteractionSpline *PairwiseInteractionSpline_Allocate                  ( CStatus                     *status                     )
    cdef  void                        PairwiseInteractionSpline_AssignElectrostaticSpline ( CPairwiseInteractionSpline  *self                       ,
                                                                                            CCubicSpline                *spline                     ,
                                                                                            CReal                        cutOff                     )
    cdef  void                        PairwiseInteractionSpline_AssignLennardJonesASpline ( CPairwiseInteractionSpline  *self                       ,
                                                                                            CCubicSpline                *spline                     ,
                                                                                            CReal                        cutOff                     )
    cdef  void                        PairwiseInteractionSpline_AssignLennardJonesBSpline ( CPairwiseInteractionSpline  *self                       ,
                                                                                            CCubicSpline                *spline                     ,
                                                                                            CReal                        cutOff                     )
    cdef  CPairwiseInteractionSpline *PairwiseInteractionSpline_Clone                     ( CPairwiseInteractionSpline  *self                       ,
                                                                                            CStatus                     *status                     )
    cdef  void                        PairwiseInteractionSpline_Deallocate                ( CPairwiseInteractionSpline **self                       )
    cdef  void                        PairwiseInteractionSpline_DeassignSplines           ( CPairwiseInteractionSpline  *self                       )
    cdef  void                        PairwiseInteractionSpline_Interactions              ( CPairwiseInteractionSpline  *self                       ,
                                                                                            CRealArray1D                *r                          ,
                                                                                            CRealArray1D                *electrostatic              ,
                                                                                            CRealArray1D                *lennardJonesA              ,
                                                                                            CRealArray1D                *lennardJonesB              )
    cdef  void                        PairwiseInteractionSpline_MMMMEnergy                ( CPairwiseInteractionSpline  *self                       ,  
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
    cdef  void                        PairwiseInteractionSpline_MMMMEnergyImage           ( CPairwiseInteractionSpline  *self                       ,  
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
    cdef  void                      PairwiseInteractionSpline_MMMMEnergyMI                ( CPairwiseInteractionSpline  *self                       ,
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
    cdef  void                        PairwiseInteractionSpline_QCMMGradients             ( CPairwiseInteractionSpline  *self                       ,  
                                                                                            CRealArray1D                *chargesQ                   ,  
                                                                                            CRealArray1D                *chargesM                   ,  
                                                                                            CReal                        electrostaticScale         ,  
                                                                                            CRealArray2D                *coordinates3Q              ,  
                                                                                            CRealArray2D                *coordinates3M              ,  
                                                                                            CPairList                   *pairList                   ,  
                                                                                            CRealArray2D                *gradients3Q                ,  
                                                                                            CRealArray2D                *gradients3M                ,  
                                                                                            CStatus                     *status                     ) 
    cdef  void                        PairwiseInteractionSpline_QCMMGradientsImage        ( CPairwiseInteractionSpline  *self                       ,  
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
    cdef  void                      PairwiseInteractionSpline_QCMMGradientsMI             ( CPairwiseInteractionSpline  *self                       ,
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
    cdef  void                        PairwiseInteractionSpline_QCMMPotentials            ( CPairwiseInteractionSpline  *self                       ,  
                                                                                            CRealArray1D                *chargesM                   ,  
                                                                                            CReal                        electrostaticScale         ,  
                                                                                            CRealArray2D                *coordinates3Q              ,  
                                                                                            CRealArray2D                *coordinates3M              ,  
                                                                                            CPairList                   *pairList                   ,  
                                                                                            CRealArray1D                *potentials                 ,  
                                                                                            CStatus                     *status                     ) 
    cdef  void                        PairwiseInteractionSpline_QCMMPotentialsImage       ( CPairwiseInteractionSpline  *self                       ,  
                                                                                            CRealArray1D                *charges                    ,  
                                                                                            CReal                        electrostaticScale         ,  
                                                                                            CRealArray2D                *coordinates3A              ,  
                                                                                            CRealArray2D                *coordinates3B              ,  
                                                                                            CSymmetryParameters         *symmetryParameters         ,  
                                                                                            CImagePairListContainer     *imagePairLists             ,  
                                                                                            CRealArray1D                *potentials                 ,  
                                                                                            CStatus                     *status                     ) 
    cdef  void                        PairwiseInteractionSpline_QCMMPotentialsMI          ( CPairwiseInteractionSpline  *self                       ,
                                                                                            CRealArray1D                *chargesM                   ,
                                                                                            CReal                        electrostaticScale         ,
                                                                                            CRealArray2D                *coordinates3Q              ,
                                                                                            CRealArray2D                *coordinates3M              ,
                                                                                            CSymmetryParameters         *symmetryParameters         ,
                                                                                            CPairList                   *pairList                   ,
                                                                                            CRealArray1D                *potentials                 ,
                                                                                            CStatus                     *status                     )
#   cdef  void                        PairwiseInteractionSpline_QCQCGradients             ( CPairwiseInteractionSpline  *self                       ,  
#                                                                                           CRealArray1D                *charges                    ,  
#                                                                                           CReal                        electrostaticScale         ,  
#                                                                                           CRealArray2D                *coordinates3I              ,  
#                                                                                           CRealArray2D                *coordinates3J              ,  
#                                                                                           CPairList                   *pairList                   ,  
#                                                                                           CRealArray2D                *gradients3I                ,  
#                                                                                           CRealArray2D                *gradients3J                ,  
#                                                                                           CStatus                     *status                     ) 
    cdef  void                        PairwiseInteractionSpline_QCQCGradientsImage        ( CPairwiseInteractionSpline  *self                       ,  
                                                                                            CRealArray1D                *charges                    ,  
                                                                                            CReal                        electrostaticScale         ,  
                                                                                            CRealArray2D                *coordinates3               ,  
                                                                                            CSymmetryParameters         *symmetryParameters         ,  
                                                                                            CImagePairListContainer     *imagePairLists             ,  
                                                                                            CRealArray2D                *gradients3                 ,  
                                                                                            CSymmetryParameterGradients *symmetryParameterGradients ,  
                                                                                            CStatus                     *status                     ) 
#   cdef  void                        PairwiseInteractionSpline_QCQCPotentials            ( CPairwiseInteractionSpline  *self                       ,  
#                                                                                           CReal                        electrostaticScale         ,  
#                                                                                           CRealArray2D                *coordinates3I              ,  
#                                                                                           CRealArray2D                *coordinates3J              ,  
#                                                                                           CPairList                   *pairList                   ,  
#                                                                                           CSymmetricMatrix            *potentials                 ,  
#                                                                                           CStatus                     *status                     ) 
    cdef  void                        PairwiseInteractionSpline_QCQCPotentialsImage       ( CPairwiseInteractionSpline  *self                       ,  
                                                                                            CReal                        electrostaticScale         ,  
                                                                                            CRealArray2D                *coordinates3               ,  
                                                                                            CSymmetryParameters         *symmetryParameters         ,  
                                                                                            CImagePairListContainer     *imagePairLists             ,  
                                                                                            CSymmetricMatrix            *potentials                 ,  
                                                                                            CStatus                     *status                     ) 

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class PairwiseInteractionSplineABFS ( PairwiseInteraction ):

    cdef CPairwiseInteractionSpline *cObject
    cdef public object               isOwner

    cdef public CubicSpline          electrostaticSpline
    cdef public CubicSpline          lennardJonesASpline
    cdef public CubicSpline          lennardJonesBSpline
    cdef public object               dampingCutOff
    cdef public object               electrostaticModel
    cdef public object               innerCutOff
    cdef public object               integrator
    cdef public object               outerCutOff
    cdef public object               pointDensity
    cdef public object               width1
    cdef public object               width2
