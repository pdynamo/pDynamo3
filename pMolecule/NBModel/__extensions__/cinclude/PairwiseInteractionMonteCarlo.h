# ifndef _PAIRWISEINTERACTIONMONTECARLO
# define _PAIRWISEINTERACTIONMONTECARLO

# include "Boolean.h"
# include "Coordinates3.h"
# include "Integer.h"
# include "IntegerArray1D.h"
# include "LJParameterContainer.h"
# include "PairwiseInteractionMonteCarlo.h"
# include "Real.h"
# include "RealArray1D.h"
# include "Selection.h"
# include "SelectionContainer.h"
# include "SymmetryParameters.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The Monte Carlo pairwise interaction type. */
typedef struct {
    Integer  isolateScale ;
    Real     buffer       ;
    Real     chargeScale  ;
    Real     cutOff       ;
    Real     epsilonScale ;
    Real     sigmaScale   ;
    Real     underFlowL   ;
    Real     underFlowQ   ;
} PairwiseInteractionMonteCarlo ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern PairwiseInteractionMonteCarlo *PairwiseInteractionMonteCarlo_Allocate                          (       Status                         *status             ) ;
extern PairwiseInteractionMonteCarlo *PairwiseInteractionMonteCarlo_Clone                             ( const PairwiseInteractionMonteCarlo  *self               ,
                                                                                                              Status                         *status             ) ;
extern void                           PairwiseInteractionMonteCarlo_Deallocate                        (       PairwiseInteractionMonteCarlo **self               ) ;
extern void                           PairwiseInteractionMonteCarlo_Interactions                      ( const PairwiseInteractionMonteCarlo  *self               ,
                                                                                                        const RealArray1D                    *r                  ,
                                                                                                              RealArray1D                    *electrostatic      ,
                                                                                                              RealArray1D                    *lennardJonesA      ,
                                                                                                              RealArray1D                    *lennardJonesB      ) ;
extern Real                           PairwiseInteractionMonteCarlo_MMMMEnergy                        ( const PairwiseInteractionMonteCarlo  *self               ,
                                                                                                              RealArray1D                    *charges            ,
                                                                                                              IntegerArray1D                 *ljTypes            ,
                                                                                                              LJParameterContainer           *ljParameters       ,
                                                                                                        const Real                            electrostaticScale ,
                                                                                                        const Real                            lennardJonesScale  ,
                                                                                                              SelectionContainer             *isolates           ,
                                                                                                              BooleanBlock                   *isFree             ,
                                                                                                              Coordinates3                   *coordinates3       ,
                                                                                                              SymmetryParameters             *symmetryParameters ,
                                                                                                              Real                           *eElectrostatic     ,
                                                                                                              Real                           *eLennardJones      ,
                                                                                                              Status                         *status             ) ;
extern Real                           PairwiseInteractionMonteCarlo_MMMMIsolateEnergy                 ( const PairwiseInteractionMonteCarlo  *self               ,
                                                                                                        const Integer                         isolate            ,
                                                                                                              RealArray1D                    *charges            ,
                                                                                                              IntegerArray1D                 *ljTypes            ,
                                                                                                              LJParameterContainer           *ljParameters       ,
                                                                                                        const Real                            electrostaticScale ,
                                                                                                        const Real                            lennardJonesScale  ,
                                                                                                              SelectionContainer             *isolates           ,
                                                                                                              BooleanBlock                   *isFree             ,
                                                                                                              Coordinates3                   *coordinates3       ,
                                                                                                              SymmetryParameters             *symmetryParameters ,
                                                                                                              Real                           *eElectrostatic     ,
                                                                                                              Real                           *eLennardJones      ,
                                                                                                              Status                         *status             ) ;
extern void                           PairwiseInteractionMonteCarlo_ScaleIsolateInteractionParameters (       PairwiseInteractionMonteCarlo  *self               ,
                                                                                                        const Integer                         isolate            ,
                                                                                                        const Real                            chargeScale        ,
                                                                                                        const Real                            epsilonScale       ,
                                                                                                        const Real                            sigmaScale         ) ;

# endif
