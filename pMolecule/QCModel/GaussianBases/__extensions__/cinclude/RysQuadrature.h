# ifndef _RYSQUADRATURE
# define _RYSQUADRATURE

# include "Integer.h"
# include "GaussianBasis.h" /* . Only needed for MAXIMUM_ANGULAR_MOMENTUM. */
# include "Real.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The maximum number of rys roots. */
/* . Sufficient for the first derivatives of TEIs with an anti-Coulomb operator. */
/* . Current method good up to _MAXRYS = 12. After this some failures, especially in range X=10-100. */
# define _MAXRYS ( ( 4 * MAXIMUM_ANGULAR_MOMENTUM + 7 ) / 2 + 1 )

/* . The Rys quadrature type. */
typedef struct
{
   Real roots  [_MAXRYS] ;
   Real weights[_MAXRYS] ;
} RysQuadrature ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Integer RysQuadrature_MaximumRoots ( void ) ;
extern void    RysQuadrature_Roots        ( RysQuadrature *roots, Integer nRoots, Real x ) ;

# endif
