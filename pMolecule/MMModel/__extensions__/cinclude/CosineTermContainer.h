# ifndef _COSINETERMCONTAINER
# define _COSINETERMCONTAINER

# include "CosineParameter.h"
# include "CosineTerm.h"
# include "Integer.h"
# include "Selection.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
typedef struct {
    Integer          nIndices    ;
    Integer          nParameters ;
    Integer          nTerms      ;
    CosineParameter *parameters  ;
    CosineTerm      *terms       ;
} CosineTermContainer ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void                 CosineTermContainer_ActivateTerms         (       CosineTermContainer  *self          ) ;
extern CosineTermContainer *CosineTermContainer_Allocate              ( const Integer               nIndices      ,
                                                                        const Integer               nTerms        ,
                                                                        const Integer               nParameters   ) ;
extern CosineTermContainer *CosineTermContainer_Clone                 ( const CosineTermContainer  *self          ) ;
extern void                 CosineTermContainer_DeactivateTerms       (       CosineTermContainer  *self          ,
                                                                              Selection            *selection     ) ;
extern void                 CosineTermContainer_Deallocate            (       CosineTermContainer **self          ) ;
extern Integer              CosineTermContainer_FindMaximumPeriod     (       CosineTermContainer  *self          ) ;
extern void                 CosineTermContainer_MakePowers            (       CosineTermContainer  *self          ) ;
extern Integer              CosineTermContainer_NumberOfInactiveTerms ( const CosineTermContainer  *self          ) ;
extern CosineTermContainer *CosineTermContainer_Prune                 (       CosineTermContainer  *self          ,
                                                                              Selection            *selection     ) ;
extern Integer              CosineTermContainer_UpperBound            (       CosineTermContainer  *self          ) ;

# endif
