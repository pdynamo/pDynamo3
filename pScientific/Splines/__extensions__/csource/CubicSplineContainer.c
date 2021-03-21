/*==================================================================================================================================
! . A container for cubic splines.
!=================================================================================================================================*/

# include "CubicSplineContainer.h"
# include "Memory.h"
# include "NumericalMacros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
CubicSplineContainer *CubicSplineContainer_Allocate ( const Integer  capacity, Status *status )
{
    CubicSplineContainer *self = Memory_AllocateType ( CubicSplineContainer ) ;
    if ( self != NULL )
    {
        self->capacity = capacity ;
        self->entries  = NULL     ;
        self->isOwner  = False    ;
        if ( capacity > 0 )
        {
            self->entries = Memory_AllocateArrayOfReferences ( capacity, CubicSpline ) ;
            if ( self->entries == NULL ) CubicSplineContainer_Deallocate ( &self ) ;
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
CubicSplineContainer *CubicSplineContainer_Clone ( const CubicSplineContainer *self, Status *status )
{
    CubicSplineContainer *clone = NULL ;
    if ( self != NULL )
    {
        clone = CubicSplineContainer_Allocate ( self->capacity, status ) ;
        if ( clone != NULL )
        {
            auto Integer  i ;
            clone->isOwner = self->isOwner ;
            if ( self->isOwner )
            {
                for ( i = 0 ; i < self->capacity ; i++ )
                {
                    clone->entries[i] = CubicSpline_Clone ( self->entries[i], status ) ;
                    if ( clone->entries[i] == NULL ) { CubicSplineContainer_Deallocate ( &clone ) ; break ; }
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
void CubicSplineContainer_Deallocate ( CubicSplineContainer **self )
{
    if ( (*self) != NULL )
    {
        if ( (*self)->isOwner )
        {
            auto Integer  i ;
            for ( i = 0 ; i < (*self)->capacity ; i++ ) CubicSpline_Deallocate ( &((*self)->entries[i]) ) ;
        }
        Memory_Deallocate ( (*self)->entries ) ;
        Memory_Deallocate ( (*self)          ) ;
    }
}
