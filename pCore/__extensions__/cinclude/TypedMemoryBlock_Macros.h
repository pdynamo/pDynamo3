# ifndef _TYPEDMEMORYBLOCKMACROS
# define _TYPEDMEMORYBLOCKMACROS

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The block capacity. */
# define Block_Capacity( self ) ( ( (self) == NULL ) ? 0 : self->capacity )

/* . An item. */
# define Block_Item( self, i ) ( (self)->items[i] )

/* . All items. */
# define Block_Items( self ) ( (self)->items )

# endif
