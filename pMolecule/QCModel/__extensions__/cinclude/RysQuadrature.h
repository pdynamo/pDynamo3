# ifndef _RYSQUADRATURE
# define _RYSQUADRATURE

# include "Integer.h"
# include "Real.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The maximum number of rys roots. */
# define _MAXRYS 9

/* . The Rys quadrature type. */
typedef struct
{
   Real roots  [_MAXRYS] ;
   Real weights[_MAXRYS] ;
} RysQuadrature ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void RysQuadrature_Roots ( RysQuadrature *roots, Integer nRoots, Real x ) ;

# endif
