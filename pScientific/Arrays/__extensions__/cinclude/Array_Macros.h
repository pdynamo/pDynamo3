# ifndef _ARRAYMACROS
# define _ARRAYMACROS

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Offset is normally either 0, self->offset, or self->view->offset. */

/* . Assign a target block to an array without a reference. */
# define Array_AssignBlockWithoutReference( self, tBlock, offset ) \
        { \
            (self)->block = NULL    ; \
            (self)->data  = &((tBlock)->items[offset]) ; \
        }

/* . Assign a target block to an array with a reference. */
# define Array_AssignBlockWithReference( self, tBlock, offset ) \
        { \
            (self)->block         = (tBlock) ; \
            (tBlock)->references += 1        ; \
            (self)->data          = &((tBlock)->items[offset]) ; \
        }

/* . A pointer to the start of the data in block. */
# define Array_BlockDataPointer( self, offset ) ( &(self->block->items[offset]) )

/* . A pointer to the start of the data. */
# define Array_DataPointer( self ) ( &((self)->data[0]) )

# endif
