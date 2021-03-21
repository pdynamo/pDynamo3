/*==================================================================================================================================
! . The default real type (at least 64 bit floating point).
!=================================================================================================================================*/
# ifndef _REAL
# define _REAL

# include <float.h>
# include "Integer.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Type. */
typedef double Real ;

/* . Largest and smallest values. */
# define Real_Largest   DBL_MAX
# define Real_Smallest -DBL_MAX

/* . Rounding function to return the nearest integer from a real. */
# define Real_RoundToInteger( a ) ( ( ( a ) >= 0 ) ? ( Integer ) ( ( a ) + 0.5 ) : ( Integer ) ( ( a ) - 0.5 ) )

/* . The smallest positive value. */
# define Real_SmallestPositive DBL_MIN

/* . Specific values. */
# define Real_Initializer 0.0e+00
# define Real_Unknown     Real_Largest

/* . Other constants - these all need to be checked. */
# define Real_SafeMinimum DBL_MIN /* . Originally { return dlamch_ ( "S" ) ; }. */

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Real Real_UnitRoundOff ( void ) ;

# endif
