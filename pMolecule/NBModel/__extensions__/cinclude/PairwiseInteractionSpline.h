# ifndef _PAIRWISEINTERACTIONSPLINE
# define _PAIRWISEINTERACTIONSPLINE

# include "Coordinates3.h"
# include "CubicSpline.h"
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
/* . The spline pairwise interaction type. */
typedef struct {
    Real         cutOff              ;
    Real         cutOff2             ;
    CubicSpline *electrostaticSpline ;
    CubicSpline *lennardJonesASpline ;
    CubicSpline *lennardJonesBSpline ;
} PairwiseInteractionSpline ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern PairwiseInteractionSpline *PairwiseInteractionSpline_Allocate                  (       Status                     *status                     ) ;
extern void                       PairwiseInteractionSpline_AssignElectrostaticSpline (       PairwiseInteractionSpline  *self                       ,
                                                                                              CubicSpline                *spline                     ,
                                                                                              Real                        cutOff                     ) ;
extern void                       PairwiseInteractionSpline_AssignLennardJonesASpline (       PairwiseInteractionSpline  *self                       ,
                                                                                              CubicSpline                *spline                     ,
                                                                                              Real                        cutOff                     ) ;
extern void                       PairwiseInteractionSpline_AssignLennardJonesBSpline (       PairwiseInteractionSpline  *self                       ,
                                                                                              CubicSpline                *spline                     ,
                                                                                              Real                        cutOff                     ) ;
extern PairwiseInteractionSpline *PairwiseInteractionSpline_Clone                     (       PairwiseInteractionSpline  *self                       ,
                                                                                              Status                     *status                     ) ;
extern void                       PairwiseInteractionSpline_Deallocate                (       PairwiseInteractionSpline **self                       ) ;
extern void                       PairwiseInteractionSpline_DeassignSplines           (       PairwiseInteractionSpline  *self                       ) ;
extern void                       PairwiseInteractionSpline_Interactions              ( const PairwiseInteractionSpline  *self                       ,
                                                                                        const RealArray1D                *r                          ,
                                                                                              RealArray1D                *electrostatic              ,
                                                                                              RealArray1D                *lennardJonesA              ,
                                                                                              RealArray1D                *lennardJonesB              ) ;
extern void                       PairwiseInteractionSpline_MMMMEnergy                ( const PairwiseInteractionSpline  *self                       ,   
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
extern void                       PairwiseInteractionSpline_MMMMEnergyImage           ( const PairwiseInteractionSpline  *self                       ,   
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
extern void                       PairwiseInteractionSpline_MMMMEnergyMI              ( const PairwiseInteractionSpline  *self                       ,
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
extern void                       PairwiseInteractionSpline_QCMMGradients             ( const PairwiseInteractionSpline  *self                       ,   
                                                                                        const RealArray1D                *chargesQ                   ,   
                                                                                        const RealArray1D                *chargesM                   ,   
                                                                                        const Real                        electrostaticScale         ,   
                                                                                        const Coordinates3               *coordinates3Q              ,   
                                                                                        const Coordinates3               *coordinates3M              ,   
                                                                                              PairList                   *pairList                   ,   
                                                                                        const Coordinates3               *gradients3Q                ,   
                                                                                        const Coordinates3               *gradients3M                ,   
                                                                                              Status                     *status                     ) ; 
extern void                       PairwiseInteractionSpline_QCMMGradientsImage        ( const PairwiseInteractionSpline  *self                       ,   
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
extern void                       PairwiseInteractionSpline_QCMMGradientsMI           ( const PairwiseInteractionSpline  *self                       ,
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
extern void                       PairwiseInteractionSpline_QCMMPotentials            ( const PairwiseInteractionSpline  *self                       ,   
                                                                                        const RealArray1D                *chargesM                   ,   
                                                                                        const Real                        electrostaticScale         ,   
                                                                                        const Coordinates3               *coordinates3Q              ,   
                                                                                        const Coordinates3               *coordinates3M              ,   
                                                                                              PairList                   *pairList                   ,   
                                                                                              RealArray1D                *potentials                 ,   
                                                                                              Status                     *status                     ) ; 
extern void                       PairwiseInteractionSpline_QCMMPotentialsImage       ( const PairwiseInteractionSpline  *self                       ,   
                                                                                        const RealArray1D                *charges                    ,   
                                                                                        const Real                        electrostaticScale         ,   
                                                                                              Coordinates3               *coordinates3A              ,   
                                                                                              Coordinates3               *coordinates3B              ,   
                                                                                              SymmetryParameters         *symmetryParameters         ,   
                                                                                              ImagePairListContainer     *imagePairLists             ,   
                                                                                              RealArray1D                *potentials                 ,   
                                                                                              Status                     *status                     ) ; 
extern void                       PairwiseInteractionSpline_QCMMPotentialsMI          ( const PairwiseInteractionSpline  *self                       ,
                                                                                        const RealArray1D                *chargesM                   ,
                                                                                        const Real                        electrostaticScale         ,
                                                                                        const Coordinates3               *coordinates3Q              ,
                                                                                        const Coordinates3               *coordinates3M              ,
                                                                                              SymmetryParameters         *symmetryParameters         ,
                                                                                              PairList                   *pairList                   ,
                                                                                              RealArray1D                *potentials                 ,
                                                                                              Status                     *status                     ) ;
extern void                       PairwiseInteractionSpline_QCQCGradients             ( const PairwiseInteractionSpline  *self                       ,   
                                                                                        const RealArray1D                *charges                    ,   
                                                                                        const Real                        electrostaticScale         ,   
                                                                                        const Coordinates3               *coordinates3I              ,   
                                                                                        const Coordinates3               *coordinates3J              ,   
                                                                                              PairList                   *pairList                   ,   
                                                                                        const Coordinates3               *gradients3I                ,   
                                                                                        const Coordinates3               *gradients3J                ,   
                                                                                              Status                     *status                     ) ; 
extern void                       PairwiseInteractionSpline_QCQCGradientsImage        ( const PairwiseInteractionSpline  *self                       ,   
                                                                                        const RealArray1D                *charges                    ,   
                                                                                        const Real                        electrostaticScale         ,   
                                                                                              Coordinates3               *coordinates3               ,   
                                                                                              SymmetryParameters         *symmetryParameters         ,   
                                                                                              ImagePairListContainer     *imagePairLists             ,   
                                                                                              Coordinates3               *gradients3                 ,   
                                                                                              SymmetryParameterGradients *symmetryParameterGradients ,   
                                                                                              Status                     *status                     ) ; 
extern void                       PairwiseInteractionSpline_QCQCPotentials            ( const PairwiseInteractionSpline  *self                       ,   
                                                                                        const Real                        electrostaticScale         ,   
                                                                                        const Coordinates3               *coordinates3I              ,   
                                                                                        const Coordinates3               *coordinates3J              ,   
                                                                                              PairList                   *pairList                   ,   
                                                                                              SymmetricMatrix            *potentials                 ,   
                                                                                              Status                     *status                     ) ; 
extern void                       PairwiseInteractionSpline_QCQCPotentialsImage       ( const PairwiseInteractionSpline  *self                       ,   
                                                                                        const Real                        electrostaticScale         ,   
                                                                                              Coordinates3               *coordinates3               ,   
                                                                                              SymmetryParameters         *symmetryParameters         ,   
                                                                                              ImagePairListContainer     *imagePairLists             ,   
                                                                                              SymmetricMatrix            *potentials                 ,   
                                                                                              Status                     *status                     ) ; 

# endif
