/*==================================================================================================================================
! . Numerical macros.
!=================================================================================================================================*/
# ifndef _NUMERICALMACROS
# define _NUMERICALMACROS

/* . Even and odd. */
# define IsEven( n ) ( ( n % 2 ) == 0 )
# define IsOdd( n )  ( ( n % 2 ) != 0 )

/* . Maximum and minimum functions. */
# define Maximum( a, b ) ( ( a ) > ( b ) ? ( a ) : ( b ) )
# define Minimum( a, b ) ( ( a ) < ( b ) ? ( a ) : ( b ) )

/* . Modulo function returning positive number always. */
# define Modulo( a, b ) ( ( ( a < 0 ) ? ( ( a % b ) + b ) : a ) % b )

/* . Return the sign of a number. */
# define Sign( a ) ( ( a ) >= 0.0 ? 1 : -1 )

# endif
