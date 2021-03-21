/*==================================================================================================================================
! . The default cardinal type (at least 32 bit).
!=================================================================================================================================*/
# ifndef _CARDINAL
# define _CARDINAL

# include <limits.h>

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Type. */
typedef unsigned int Cardinal ;

/* . Largest and smallest values. */
# define Cardinal_Largest  UINT_MAX
# define Cardinal_Smallest 0

/* . Specific values. */
# define Cardinal_Initializer Cardinal_Smallest
# define Cardinal_Unknown     Cardinal_Largest

# endif
