# ifndef _ARRAY2D_N3MACROS
# define _ARRAY2D_N3MACROS

# include "Array2D_Macros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define Array2D_DecrementRowN3( self, i, xij, yij, zij ) \
        { \
	    Array2D_Item ( self, i, 0 ) -= xij ; \
	    Array2D_Item ( self, i, 1 ) -= yij ; \
	    Array2D_Item ( self, i, 2 ) -= zij ; \
        }

# define Array2D_DifferenceRowN3( self, i, j, xij, yij, zij ) \
        { \
	    xij = Array2D_Item ( self, i, 0 ) - Array2D_Item ( self, j, 0 ) ; \
	    yij = Array2D_Item ( self, i, 1 ) - Array2D_Item ( self, j, 1 ) ; \
	    zij = Array2D_Item ( self, i, 2 ) - Array2D_Item ( self, j, 2 ) ; \
        }

# define Array2D_GetRowN3( self, i, x, y, z ) \
        { \
	   x = Array2D_Item ( self, i, 0 ) ; \
	   y = Array2D_Item ( self, i, 1 ) ; \
	   z = Array2D_Item ( self, i, 2 ) ; \
        }

# define Array2D_IncrementRowN3( self, i, xij, yij, zij ) \
        { \
	    Array2D_Item ( self, i, 0 ) += xij ; \
	    Array2D_Item ( self, i, 1 ) += yij ; \
	    Array2D_Item ( self, i, 2 ) += zij ; \
        }

# define Array2D_ScaleRowN3( self, i, value ) \
        { \
	    Array2D_Item ( self, i, 0 ) *= value ; \
	    Array2D_Item ( self, i, 1 ) *= value ; \
	    Array2D_Item ( self, i, 2 ) *= value ; \
        }

# define Array2D_SetRowN3( self, i, x, y, z ) \
        { \
	    Array2D_Item ( self, i, 0 ) = x ; \
	    Array2D_Item ( self, i, 1 ) = y ; \
	    Array2D_Item ( self, i, 2 ) = z ; \
        }

# endif
