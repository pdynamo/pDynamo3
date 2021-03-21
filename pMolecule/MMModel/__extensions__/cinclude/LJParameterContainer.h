# ifndef _LJPARAMETERCONTAINER
# define _LJPARAMETERCONTAINER

# include "Integer.h"
# include "Real.h"

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
/* . The Lennard-Jones parameter container type. */
typedef struct {
    Integer  ntypes     ;
    Integer *tableindex ;
    Real    *epsilon    ;
    Real    *sigma      ;
    Real    *tableA     ;
    Real    *tableB     ;
} LJParameterContainer ;


/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
extern LJParameterContainer *LJParameterContainer_Allocate              ( const Integer                ntypes ) ;
extern LJParameterContainer *LJParameterContainer_Clone                 ( const LJParameterContainer  *self   ) ;
extern void                  LJParameterContainer_Deallocate            (       LJParameterContainer **self   ) ;
extern void                  LJParameterContainer_MakeEpsilonSigmaAMBER (       LJParameterContainer  *self   ) ;
extern void                  LJParameterContainer_MakeEpsilonSigmaOPLS  (       LJParameterContainer  *self   ) ;
extern void                  LJParameterContainer_MakeTableAMBER        (       LJParameterContainer  *self   ) ;
extern void                  LJParameterContainer_MakeTableOPLS         (       LJParameterContainer  *self   ) ;
extern LJParameterContainer *LJParameterContainer_MergeEpsilonSigma     ( const LJParameterContainer  *self   ,
                                                                          const LJParameterContainer  *other  ) ;
extern void                  LJParameterContainer_Scale                 (       LJParameterContainer  *self   ,
                                                                          const Real                   scale  ) ;

# endif
