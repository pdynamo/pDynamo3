/*==================================================================================================================================
! . 1-D array view.
!=================================================================================================================================*/

# include "IntegerUtilities.h"
# include "Iterator1D.h"
# include "Memory.h"
# include "NumericalMacros.h"
# include "View1D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation of type.
!---------------------------------------------------------------------------------------------------------------------------------*/
View1D *View1D_Allocate ( Status *status )
{
    View1D *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( View1D ) ;
        if ( self != NULL )
        {
            View1D_InitializeFields ( self ) ;
        }
        else Status_Set ( status, Status_OutOfMemory ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check the view versus a capacity.
! . Both maximum and minimum indices are calculated and checked.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean View1D_CheckCapacity ( const View1D *self, const Integer capacity )
{
    Boolean isOK = False ;
    if ( self != NULL )
    {
        auto Integer m = 0, n = 0, p ;
        if ( self->size > 0 )
        {
            m = self->offset ;
            n = self->offset ;
            p = ( self->extent - 1 ) * self->stride ;
            if ( p < 0 ) m += p ;
            else         n += p ;
        }
        isOK = ( ( m >= 0 ) && ( n < capacity ) ) ;
    }
    return isOK ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check two views for conformability.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean View1D_CheckConformability ( const View1D *self, const View1D *other, Status *status )
{
    Boolean areConformable = False ;
    if ( ( self != NULL ) && ( other != NULL ) && Status_IsOK ( status ) )
    {
        areConformable = ( self->extent == other->extent ) ;
        if ( ! areConformable ) Status_Set ( status, Status_NonConformableArrays ) ;
    }
    return areConformable ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copying.
!---------------------------------------------------------------------------------------------------------------------------------*/
void View1D_CopyTo ( const View1D *self, View1D *other )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        other->extent = self->extent ;
        other->offset = self->offset ;
        other->size   = self->size   ;
        other->stride = self->stride ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void View1D_Deallocate ( View1D **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) ) Memory_Deallocate ( (*self) ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Extent.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer View1D_GetExtent ( const View1D *self ) { return ( ( (self) == NULL ) ? 0 : (self)->extent ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Offset.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer View1D_GetOffset ( const View1D *self ) { return ( ( (self) == NULL ) ? 0 : (self)->offset ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Stride.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer View1D_GetStride ( const View1D *self ) { return ( ( (self) == NULL ) ? 0 : (self)->stride ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization of a compact view given a shape.
!---------------------------------------------------------------------------------------------------------------------------------*/
void View1D_Initialize ( View1D *self, const Integer extent, Status *status )
{
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        if ( extent < 0 ) Status_Set ( status, Status_InvalidArgument ) ;
        self->extent = Maximum ( extent, 0 ) ;
        self->offset = 0 ;
        self->size   = self->extent ;
        self->stride = 1 ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the default iterator (flattened as much as possible).
!---------------------------------------------------------------------------------------------------------------------------------*/
Iterator *View1D_MakeIterator ( const View1D *self, Status *status )
{
    Iterator *iterator = NULL ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto Iterator1D *iterator1D ;
        iterator   = Iterator_Allocate   ( status ) ;
        iterator1D = Iterator1D_Allocate ( status ) ;
        Iterator1D_Initialize   ( iterator1D, self->offset, self->extent, self->stride, status ) ;
        Iterator1D_MakeIterator ( iterator1D, iterator ) ;
    }
    return iterator ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set the state of the view - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
void View1D_SetState (       View1D  *self   ,
                       const Integer  extent ,
                       const Integer  offset ,
                       const Integer  size   ,
                       const Integer  stride )
{
    if ( self != NULL )
    {
        self->extent = extent ;
        self->offset = offset ;
        self->size   = size   ;
        self->stride = stride ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . View.
!---------------------------------------------------------------------------------------------------------------------------------*/
void View1D_View ( const View1D  *self   ,
                   const Integer  start  ,
                   const Integer  extent ,
                   const Integer  stride ,
                         View1D  *view   ,
                         Status  *status )
{
    if ( ( self != NULL ) && ( view != NULL ) && Status_IsOK ( status ) )
    {
        view->extent = extent ;
        view->offset = self->offset + start * self->stride ;
        view->size   = extent ;
        view->stride = stride * self->stride ;
    }
}
