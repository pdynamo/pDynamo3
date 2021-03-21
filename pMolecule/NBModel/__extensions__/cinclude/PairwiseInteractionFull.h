# ifndef _PAIRWISEINTERACTIONFULL
# define _PAIRWISEINTERACTIONFULL

# include "Coordinates3.h"
# include "IntegerArray1D.h"
# include "LJParameterContainer.h"
# include "PairList.h"
# include "Real.h"
# include "RealArray1D.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The full pairwise interaction type. */
typedef struct {
    Real dampingCutOff ;
    Real alpha1        ;
    Real alpha2        ;
    Real alpha3        ;
    Real alpha6        ;
    Real alpha12       ;
    Real beta1         ;
    Real beta2         ;
    Real beta3         ;
    Real beta6         ;
    Real beta12        ;
} PairwiseInteractionFull  ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern PairwiseInteractionFull *PairwiseInteractionFull_Allocate            (       Status                   *status             ) ;
extern PairwiseInteractionFull *PairwiseInteractionFull_Clone               (       PairwiseInteractionFull  *self               ,
                                                                                    Status                   *status             ) ;
extern void                     PairwiseInteractionFull_Deallocate          (       PairwiseInteractionFull **self               ) ;
extern void                     PairwiseInteractionFull_InitializeDependent (       PairwiseInteractionFull  *self               ) ;
extern void                     PairwiseInteractionFull_Interactions        ( const PairwiseInteractionFull  *self               ,
                                                                              const RealArray1D              *r                  ,
                                                                                    RealArray1D              *lennardJonesA      ,
                                                                                    RealArray1D              *lennardJonesB      ,
                                                                                    RealArray1D              *multipole0         ,
                                                                                    RealArray1D              *multipole1         ,
                                                                                    RealArray1D              *multipole2         ) ;
extern void                     PairwiseInteractionFull_MMMMEnergy          ( const PairwiseInteractionFull  *self               ,
                                                                              const RealArray1D              *chargesI           ,
                                                                              const RealArray1D              *chargesJ           ,
                                                                                    IntegerArray1D           *ljTypesI           ,
                                                                                    IntegerArray1D           *ljTypesJ           ,
                                                                              const LJParameterContainer     *ljParameters       ,
                                                                              const Real                      electrostaticScale ,
                                                                              const Real                      lennardJonesScale  ,
                                                                              const Coordinates3             *coordinates3I      ,
                                                                              const Coordinates3             *coordinates3J      ,
                                                                                    PairList                 *pairList           ,
                                                                                    Real                     *eElectrostatic     ,
                                                                                    Real                     *eLennardJones      ,
                                                                                    Coordinates3             *gradients3I        ,
                                                                                    Coordinates3             *gradients3J        ,
                                                                                    Status                   *status             ) ;
extern void                     PairwiseInteractionFull_QCMMGradients       ( const PairwiseInteractionFull  *self               ,
                                                                              const Integer                   multipoleOrder     ,
                                                                              const RealArray1D              *multipolesQ        ,
                                                                              const RealArray1D              *chargesM           ,
                                                                              const Real                      electrostaticScale ,
                                                                              const Coordinates3             *coordinates3Q      ,
                                                                              const Coordinates3             *coordinates3M      ,
                                                                                    PairList                 *pairList           ,
                                                                              const Coordinates3             *gradients3Q        ,
                                                                              const Coordinates3             *gradients3M        ,
                                                                                    Status                   *status             ) ;
extern void                     PairwiseInteractionFull_QCMMPotentials      ( const PairwiseInteractionFull  *self               ,
                                                                              const Integer                   multipoleOrder     ,
                                                                              const RealArray1D              *chargesM           ,
                                                                              const Real                      electrostaticScale ,
                                                                              const Coordinates3             *coordinates3Q      ,
                                                                              const Coordinates3             *coordinates3M      ,
                                                                                    PairList                 *pairList           ,
                                                                                    RealArray1D              *potentials         ,
                                                                                    Status                   *status             ) ;
# endif
