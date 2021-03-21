# ifndef _CUBICSPLINECONTAINER
# define _CUBICSPLINECONTAINER

# include "Boolean.h"
# include "CubicSpline.h"
# include "Integer.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The container type. */
typedef struct {
    Boolean       isOwner  ;
    Integer       capacity ;
    CubicSpline **entries  ;
} CubicSplineContainer ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern CubicSplineContainer *CubicSplineContainer_Allocate   ( const Integer                capacity ,
                                                                     Status                *status   ) ;
extern CubicSplineContainer *CubicSplineContainer_Clone      ( const CubicSplineContainer  *self     ,
                                                                     Status                *status   ) ;
extern void                  CubicSplineContainer_Deallocate (       CubicSplineContainer **self     ) ;

# endif
