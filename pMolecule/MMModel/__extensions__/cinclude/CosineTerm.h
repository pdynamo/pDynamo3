# ifndef _COSINETERM
# define _COSINETERM

# include "Boolean.h"
# include "Integer.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
typedef struct {
    Boolean  isActive ;
    Integer  type     ;
    Integer *indices  ;
} CosineTerm ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void CosineTerm_Allocate   (       CosineTerm *self     ,
                                    const Integer     nIndices ) ;
extern void CosineTerm_Clone      (       CosineTerm *self     ,
                                    const CosineTerm *other    ,
                                    const Integer     nIndices ) ;
extern void CosineTerm_Deallocate (       CosineTerm *self     ) ;
extern void CosineTerm_Initialize (       CosineTerm *self     ) ;

# endif
