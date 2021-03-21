# ifndef _ARRAY1D_MACROS
# define _ARRAY1D_MACROS

# include "View1D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . A pointer to the start of the data. */
# define Array1D_Data( self ) ( &((self)->data[0]) )

/* . An item. */
# define Array1D_Item( self, i ) ( (self)->data[View1D_ItemIndex( self, i )] )

/* . A pointer to an item. */
# define Array1D_ItemPointer( self, i ) ( &((self)->data[View1D_ItemIndex( self, i )]) )

# endif
