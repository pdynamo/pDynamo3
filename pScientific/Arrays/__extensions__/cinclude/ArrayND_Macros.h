# ifndef _ARRAYND_MACROS
# define _ARRAYND_MACROS

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Extents. */
# define ArrayND_Extent0( self ) ( ( (self) == NULL ) || ( (self)->view == NULL ) ? 0 : (self)->view->extents[0] )
# define ArrayND_Extent1( self ) ( ( (self) == NULL ) || ( (self)->view == NULL ) ? 0 : (self)->view->extents[1] )
# define ArrayND_Extent2( self ) ( ( (self) == NULL ) || ( (self)->view == NULL ) ? 0 : (self)->view->extents[2] )
# define ArrayND_Extent3( self ) ( ( (self) == NULL ) || ( (self)->view == NULL ) ? 0 : (self)->view->extents[3] )

/* . The array block. */
# define ArrayND_Block( self ) ( ( (self) == NULL ) ? NULL : (self)->block )

/* . A pointer to the start of the data. */
# define ArrayND_Data( self ) ( &((self)->data[0]) )

/* . Items. */
/* . These can be used on arrays of higher dimension as all higher arrays indices will have the value 0. */
# define ArrayND_Item1D( self, i          ) ( (self)->data[ ViewND_Index1D ( (self)->view, i          ) ] )
# define ArrayND_Item2D( self, i, j       ) ( (self)->data[ ViewND_Index2D ( (self)->view, i, j       ) ] )
# define ArrayND_Item3D( self, i, j, k    ) ( (self)->data[ ViewND_Index3D ( (self)->view, i, j, k    ) ] )
# define ArrayND_Item4D( self, i, j, k, l ) ( (self)->data[ ViewND_Index4D ( (self)->view, i, j, k, l ) ] )

/* . Pointers to items. */
# define ArrayND_ItemPointer1D( self, i          ) ( &( ArrayND_Item1D( self, i          ) ) )
# define ArrayND_ItemPointer2D( self, i, j       ) ( &( ArrayND_Item2D( self, i, j       ) ) )
# define ArrayND_ItemPointer3D( self, i, j, k    ) ( &( ArrayND_Item3D( self, i, j, k    ) ) )
# define ArrayND_ItemPointer4D( self, i, j, k, l ) ( &( ArrayND_Item4D( self, i, j, k, l ) ) )

/* . The rank of the array. */
# define ArrayND_Rank( self ) ( ( (self) == NULL ) || ( (self)->view == NULL ) ? 0 : (self)->view->rank )

/* . The size of the array. */
# define ArrayND_Size( self ) ( ( (self) == NULL ) || ( (self)->view == NULL ) ? 0 : (self)->view->size )

/* . The array view. */
# define ArrayND_View( self ) ( ( (self) == NULL ) ? NULL : (self)->view )

# endif
