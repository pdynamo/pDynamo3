# ifndef _F2C
# define _F2C

/* . MJF CHANGE: long int -> int. */
/* . This is done so as to resolve the incompatibility between the cblas interface that
!    uses int and the f2c/clapack/blas interface that uses integer. On 32-bit machines
!    these are the same but on 64-bit machines different. This is the easiest fix but
!    may cause problems for interfacing with Fortran. An alternative would be to change
!    the cblas interface to integer.
*/

typedef          int  integer ;
typedef unsigned int uinteger ; /* . Only in lbitshft. */

/* . ORIGINAL. */
/*
typedef long int integer ;
typedef unsigned long uinteger ;
*/

typedef short int shortint ;
typedef float           real ;
typedef double    doublereal ;
typedef struct {       real r, i ; } complex ;
typedef struct { doublereal r, i ; } doublecomplex ;
typedef long int  logical ;
typedef short int shortlogical ;
typedef char      logical1 ;
typedef char      integer1 ;

#define TRUE_ (1)
#define FALSE_ (0)

typedef char    *address ;
typedef long int ftnint  ;
typedef long int ftnlen  ;

#define VOID void

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (doublereal)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (doublereal)min(a,b)
#define dmax(a,b) (doublereal)max(a,b)

/* . MJF - The following have been changed to remove compiler warnings.
! .  Only L_fp is used (in dgees, dgeesx, dgges, dggesx) and has been
! .  explicitly defined for these cases (the first two have 2 doublereal
! .  arguments and the last two, 3).
*/
typedef logical (*L_fp)(doublereal *a,...) ;

# endif
