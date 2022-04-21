/*==================================================================================================================================
! . A container for block storages.
!=================================================================================================================================*/

# include "BlockStorageContainer.h"
# include "Memory.h"
# include "NumericalMacros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
BlockStorageContainer *BlockStorageContainer_Allocate ( const Integer capacity, Status *status )
{
    BlockStorageContainer *self = Memory_AllocateType ( BlockStorageContainer ) ;
    if ( self != NULL )
    {
        self->capacity = capacity ;
        self->entries  = NULL     ;
        self->isOwner  = True     ; /* . By default True (c.f. other containers). */
        if ( capacity > 0 )
        {
            self->entries = Memory_AllocateArrayOfReferences ( capacity, BlockStorage ) ;
            if ( self->entries == NULL ) BlockStorageContainer_Deallocate ( &self ) ;
            else
            {
                auto Integer  i ;
                for ( i = 0 ; i < capacity ; i++ ) self->entries[i] = NULL ;
            }
        }
    }
    if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void BlockStorageContainer_Deallocate ( BlockStorageContainer **self )
{
    if ( (*self) != NULL )
    {
        if ( (*self)->isOwner )
        {
            auto Integer  i ;
            for ( i = 0 ; i < (*self)->capacity ; i++ ) BlockStorage_Deallocate ( &((*self)->entries[i]) ) ;
        }
        Memory_Deallocate ( (*self)->entries ) ;
        Memory_Deallocate ( (*self)          ) ;
    }
}
