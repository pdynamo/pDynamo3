# ifndef _ARRAY2D_MACROS
# define _ARRAY2D_MACROS

# include "View2D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . A pointer to a column. */
# define Array2D_ColumnPointer( self, i ) Array2D_ItemPointer( self, 0, i )

/* . A pointer to the start of the data. */
# define Array2D_Data( self ) ( &((self)->data[0]) )

/* . An item. */
# define Array2D_Item( self, i, j ) ( (self)->data[View2D_ItemIndex( self, i, j )] )

/* . An item by index - needed by Transpose only. */
# define Array2D_ItemByIndex( self, ij ) ( (self)->data[ij] )

/* . A pointer to an item. */
# define Array2D_ItemPointer( self, i, j ) ( &((self)->data[View2D_ItemIndex( self, i, j )]) )

/* . A pointer to a row. */
# define Array2D_RowPointer( self, i ) Array2D_ItemPointer( self, i, 0 )

# endif
