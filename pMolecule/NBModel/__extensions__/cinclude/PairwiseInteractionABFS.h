# ifndef _PAIRWISEINTERACTIONABFS
# define _PAIRWISEINTERACTIONABFS

# include "Coordinates3.h"
# include "ImagePairListContainer.h"
# include "IntegerArray1D.h"
# include "LJParameterContainer.h"
# include "PairList.h"
# include "Real.h"
# include "RealArray1D.h"
# include "Status.h"
# include "SymmetryParameterGradients.h"
# include "SymmetryParameters.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The ABFS pairwise interaction type. */
typedef struct {
    Real dampingCutOff ;
    Real innerCutOff   ;
    Real outerCutOff   ;
} PairwiseInteractionABFS ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern PairwiseInteractionABFS *PairwiseInteractionABFS_Allocate            (       Status                     *status                     ) ;
extern PairwiseInteractionABFS *PairwiseInteractionABFS_Clone               (       PairwiseInteractionABFS    *self                       ,
                                                                                    Status                     *status                     ) ;
extern void                     PairwiseInteractionABFS_Deallocate          (       PairwiseInteractionABFS   **self                       ) ;
extern void                     PairwiseInteractionABFS_Interactions        ( const PairwiseInteractionABFS    *self                       ,
                                                                              const RealArray1D                *r                          ,
                                                                                    RealArray1D                *electrostatic              ,
                                                                                    RealArray1D                *lennardJonesA              ,
                                                                                    RealArray1D                *lennardJonesB              ) ;
extern void                     PairwiseInteractionABFS_MMMMEnergy          ( const PairwiseInteractionABFS    *self                       ,
                                                                              const RealArray1D                *chargesI                   ,
                                                                              const RealArray1D                *chargesJ                   ,
                                                                                    IntegerArray1D             *ljTypesI                   ,
                                                                                    IntegerArray1D             *ljTypesJ                   ,
                                                                              const LJParameterContainer       *ljParameters               ,
                                                                              const Real                        electrostaticScale         ,
                                                                              const Real                        lennardJonesScale          ,
                                                                              const Coordinates3               *coordinates3I              ,
                                                                              const Coordinates3               *coordinates3J              ,
                                                                                    PairList                   *pairList                   ,
                                                                                    Real                       *eElectrostatic             ,
                                                                                    Real                       *eLennardJones              ,
                                                                                    Coordinates3               *gradients3I                ,
                                                                                    Coordinates3               *gradients3J                ,
                                                                                    Status                     *status                     ) ;
extern void                     PairwiseInteractionABFS_MMMMEnergyImage     ( const PairwiseInteractionABFS    *self                       ,
                                                                              const RealArray1D                *charges                    ,
                                                                                    IntegerArray1D             *ljTypes                    ,
                                                                              const LJParameterContainer       *ljParameters               ,
                                                                              const Real                        electrostaticScale         ,
                                                                                    Coordinates3               *coordinates3               ,
                                                                                    SymmetryParameters         *symmetryParameters         ,
                                                                                    ImagePairListContainer     *imagePairLists             ,
                                                                                    Real                       *eElectrostatic             ,
                                                                                    Real                       *eLennardJones              ,
                                                                                    Coordinates3               *gradients3                 ,
                                                                                    SymmetryParameterGradients *symmetryParameterGradients ,
                                                                                    Status                     *status                     ) ;
extern void                     PairwiseInteractionABFS_MMMMEnergyMI        ( const PairwiseInteractionABFS    *self                       ,
                                                                              const RealArray1D                *chargesI                   ,
                                                                              const RealArray1D                *chargesJ                   ,
                                                                                    IntegerArray1D             *ljTypesI                   ,
                                                                                    IntegerArray1D             *ljTypesJ                   ,
                                                                              const LJParameterContainer       *ljParameters               ,
                                                                              const Real                        electrostaticScale         ,
                                                                              const Real                        lennardJonesScale          ,
                                                                              const Coordinates3               *coordinates3I              ,
                                                                              const Coordinates3               *coordinates3J              ,
                                                                                    SymmetryParameters         *symmetryParameters         ,
                                                                                    PairList                   *pairList                   ,
                                                                                    Real                       *eElectrostatic             ,
                                                                                    Real                       *eLennardJones              ,
                                                                                    Coordinates3               *gradients3I                ,
                                                                                    Coordinates3               *gradients3J                ,
                                                                                    SymmetryParameterGradients *symmetryParameterGradients ,
                                                                                    Status                     *status                     ) ;
extern void                     PairwiseInteractionABFS_QCMMGradients       ( const PairwiseInteractionABFS    *self                       ,
                                                                              const RealArray1D                *chargesQ                   ,
                                                                              const RealArray1D                *chargesM                   ,
                                                                              const Real                        electrostaticScale         ,
                                                                              const Coordinates3               *coordinates3Q              ,
                                                                              const Coordinates3               *coordinates3M              ,
                                                                                    PairList                   *pairList                   ,
                                                                              const Coordinates3               *gradients3Q                ,
                                                                              const Coordinates3               *gradients3M                ,
                                                                                    Status                     *status                     ) ;
extern void                     PairwiseInteractionABFS_QCMMGradientsImage  ( const PairwiseInteractionABFS    *self                       ,
                                                                              const RealArray1D                *chargesA                   ,
                                                                              const RealArray1D                *chargesB                   ,
                                                                              const Real                        electrostaticScale         ,
                                                                                    Coordinates3               *coordinates3A              ,
                                                                                    Coordinates3               *coordinates3B              ,
                                                                                    SymmetryParameters         *symmetryParameters         ,
                                                                                    ImagePairListContainer     *imagePairLists             ,
                                                                                    Coordinates3               *gradients3A                ,
                                                                                    Coordinates3               *gradients3B                ,
                                                                                    SymmetryParameterGradients *symmetryParameterGradients ,
                                                                                    Status                     *status                     ) ;
extern void                     PairwiseInteractionABFS_QCMMGradientsMI     ( const PairwiseInteractionABFS    *self                       ,
                                                                              const RealArray1D                *chargesQ                   ,
                                                                              const RealArray1D                *chargesM                   ,
                                                                              const Real                        electrostaticScale         ,
                                                                              const Coordinates3               *coordinates3Q              ,
                                                                              const Coordinates3               *coordinates3M              ,
                                                                                    SymmetryParameters         *symmetryParameters         ,
                                                                                    PairList                   *pairList                   ,
                                                                              const Coordinates3               *gradients3Q                ,
                                                                              const Coordinates3               *gradients3M                ,
                                                                                    SymmetryParameterGradients *symmetryParameterGradients ,
                                                                                    Status                     *status                     ) ;
extern void                     PairwiseInteractionABFS_QCMMPotentials      ( const PairwiseInteractionABFS    *self                       ,
                                                                              const RealArray1D                *chargesM                   ,
                                                                              const Real                        electrostaticScale         ,
                                                                              const Coordinates3               *coordinates3Q              ,
                                                                              const Coordinates3               *coordinates3M              ,
                                                                                    PairList                   *pairList                   ,
                                                                                    RealArray1D                *potentials                 ,
                                                                                    Status                     *status                     ) ;
extern void                     PairwiseInteractionABFS_QCMMPotentialsImage ( const PairwiseInteractionABFS    *self                       ,
                                                                              const RealArray1D                *charges                    ,
                                                                              const Real                        electrostaticScale         ,
                                                                                    Coordinates3               *coordinates3A              ,
                                                                                    Coordinates3               *coordinates3B              ,
                                                                                    SymmetryParameters         *symmetryParameters         ,
                                                                                    ImagePairListContainer     *imagePairLists             ,
                                                                                    RealArray1D                *potentials                 ,
                                                                                    Status                     *status                     ) ;
extern void                     PairwiseInteractionABFS_QCMMPotentialsMI    ( const PairwiseInteractionABFS    *self                       ,
                                                                              const RealArray1D                *chargesM                   ,
                                                                              const Real                        electrostaticScale         ,
                                                                              const Coordinates3               *coordinates3Q              ,
                                                                              const Coordinates3               *coordinates3M              ,
                                                                                    SymmetryParameters         *symmetryParameters         ,
                                                                                    PairList                   *pairList                   ,
                                                                                    RealArray1D                *potentials                 ,
                                                                                    Status                     *status                     ) ;
extern void                     PairwiseInteractionABFS_QCQCGradients       ( const PairwiseInteractionABFS    *self                       ,
                                                                              const RealArray1D                *charges                    ,
                                                                              const Real                        electrostaticScale         ,
                                                                              const Coordinates3               *coordinates3I              ,
                                                                              const Coordinates3               *coordinates3J              ,
                                                                                    PairList                   *pairList                   ,
                                                                              const Coordinates3               *gradients3I                ,
                                                                              const Coordinates3               *gradients3J                ,
                                                                                    Status                     *status                     ) ;
extern void                     PairwiseInteractionABFS_QCQCGradientsImage  ( const PairwiseInteractionABFS    *self                       ,
                                                                              const RealArray1D                *charges                    ,
                                                                              const Real                        electrostaticScale         ,
                                                                                    Coordinates3               *coordinates3               ,
                                                                                    SymmetryParameters         *symmetryParameters         ,
                                                                                    ImagePairListContainer     *imagePairLists             ,
                                                                                    Coordinates3               *gradients3                 ,
                                                                                    SymmetryParameterGradients *symmetryParameterGradients ,
                                                                                    Status                     *status                     ) ;
extern void                     PairwiseInteractionABFS_QCQCPotentials      ( const PairwiseInteractionABFS    *self                       ,
                                                                              const Real                        electrostaticScale         ,
                                                                              const Coordinates3               *coordinates3I              ,
                                                                              const Coordinates3               *coordinates3J              ,
                                                                                    PairList                   *pairList                   ,
                                                                                    SymmetricMatrix            *potentials                 ,
                                                                                    Status                     *status                     ) ;
extern void                     PairwiseInteractionABFS_QCQCPotentialsImage ( const PairwiseInteractionABFS    *self                       ,
                                                                              const Real                        electrostaticScale         ,
                                                                                    Coordinates3               *coordinates3               ,
                                                                                    SymmetryParameters         *symmetryParameters         ,
                                                                                    ImagePairListContainer     *imagePairLists             ,
                                                                                    SymmetricMatrix            *potentials                 ,
                                                                                    Status                     *status                     ) ;

# endif
