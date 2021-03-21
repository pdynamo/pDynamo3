/*==================================================================================================================================
! . A container for MNDO parameter sets.
!=================================================================================================================================*/

# include "MNDOParametersContainer.h"
# include "Memory.h"
# include "NumericalMacros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
MNDOParametersContainer *MNDOParametersContainer_Allocate ( const Integer  capacity, Status *status )
{
    MNDOParametersContainer *self = Memory_AllocateType ( MNDOParametersContainer ) ;
    if ( self != NULL )
    {
        self->capacity = capacity ;
        self->entries  = NULL     ;
        self->isOwner  = False    ;
        if ( capacity > 0 )
        {
            self->entries = Memory_AllocateArrayOfReferences ( capacity, MNDOParameters ) ;
            if ( self->entries == NULL ) MNDOParametersContainer_Deallocate ( &self ) ;
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
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
MNDOParametersContainer *MNDOParametersContainer_Clone ( const MNDOParametersContainer *self, Status *status )
{
    MNDOParametersContainer *clone = NULL ;
    if ( self != NULL )
    {
        clone = MNDOParametersContainer_Allocate ( self->capacity, status ) ;
        if ( clone != NULL )
        {
            auto Integer  i ;
            clone->isOwner = self->isOwner ;
            if ( self->isOwner )
            {
                for ( i = 0 ; i < self->capacity ; i++ )
                {
                    clone->entries[i] = MNDOParameters_Clone ( self->entries[i] ) ;
                    if ( clone->entries[i] == NULL ) { MNDOParametersContainer_Deallocate ( &clone ) ; break ; }
                }
            }
            else
            {
                for ( i = 0 ; i < self->capacity ; i++ ) { clone->entries[i] = self->entries[i] ; }
            }
        }
        if ( clone == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDOParametersContainer_Deallocate ( MNDOParametersContainer **self )
{
    if ( (*self) != NULL )
    {
        if ( (*self)->isOwner )
        {
            auto Integer  i ;
            for ( i = 0 ; i < (*self)->capacity ; i++ ) MNDOParameters_Deallocate ( &((*self)->entries[i]) ) ;
        }
        Memory_Deallocate ( (*self)->entries ) ;
        Memory_Deallocate ( (*self)          ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Largest basis.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer MNDOParametersContainer_LargestBasis ( const MNDOParametersContainer *self )
{
    Integer n = 0 ;
    if ( self != NULL )
    {
        auto Integer i ;
        for ( i = 0 ; i < self->capacity ; i++ ) { n = Maximum ( n, self->entries[i]->norbitals ) ; }
    }
    return n ;
}
