/*==================================================================================================================================
! . 2-D array view.
!=================================================================================================================================*/

# include "IntegerUtilities.h"
# include "IteratorND.h"
# include "Memory.h"
# include "NumericalMacros.h"
# include "RowIterator2D.h"
# include "View2D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation of type.
!---------------------------------------------------------------------------------------------------------------------------------*/
View2D *View2D_Allocate ( Status *status )
{
    View2D *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( View2D ) ;
        if ( self != NULL )
        {
            View2D_InitializeFields ( self ) ;
        }
        else Status_Set ( status, Status_OutOfMemory ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check the view versus a capacity.
! . Both maximum and minimum indices are calculated and checked.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean View2D_CheckCapacity ( const View2D *self, const Integer capacity )
{
    Boolean isOK = False ;
    if ( self != NULL )
    {
        auto Integer m = 0, n = 0, p ;
        if ( self->size > 0 )
        {
            m = self->offset ;
            n = self->offset ;
            p = ( self->extent0 - 1 ) * self->stride0 ;
            if ( p < 0 ) m += p ;
            else         n += p ;
            p = ( self->extent1 - 1 ) * self->stride1 ;
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
Boolean View2D_CheckConformability ( const View2D *self, const View2D *other, Status *status )
{
    Boolean areConformable = False ;
    if ( ( self != NULL ) && ( other != NULL ) && Status_IsOK ( status ) )
    {
        areConformable = ( self->extent0 == other->extent0 ) && ( self->extent1 == other->extent1 ) ;
        if ( ! areConformable ) Status_Set ( status, Status_NonConformableArrays ) ;
    }
    return areConformable ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Column view.
!---------------------------------------------------------------------------------------------------------------------------------*/
void View2D_ColumnView ( const View2D  *self   ,
                         const Integer  i      ,
                               View1D  *view   ,
                               Status  *status )
{
    if ( ( self != NULL ) && ( view != NULL ) && Status_IsOK ( status ) )
    {
        if ( Integer_IsInRange ( i, 0, self->extent1 ) )
        {
            view->extent  = self->extent0 ;
            view->offset  = self->offset + i * self->stride1 ;
            view->size    = self->extent0 ;
            view->stride  = self->stride0 ;
        }
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copying.
!---------------------------------------------------------------------------------------------------------------------------------*/
void View2D_CopyTo ( const View2D *self, View2D *other )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        other->extent0 = self->extent0 ;
        other->extent1 = self->extent1 ;
        other->offset  = self->offset  ;
        other->size    = self->size    ;
        other->stride0 = self->stride0 ;
        other->stride1 = self->stride1 ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void View2D_Deallocate ( View2D **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) ) Memory_Deallocate ( (*self) ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Columns.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer View2D_GetColumns ( const View2D *self ) { return ( ( (self) == NULL ) ? 0 : (self)->extent1 ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Offset.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer View2D_GetOffset ( const View2D *self ) { return ( ( (self) == NULL ) ? 0 : (self)->offset ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Rows.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer View2D_GetRows ( const View2D *self ) { return ( ( (self) == NULL ) ? 0 : (self)->extent0 ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Size.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer View2D_GetSize ( const View2D *self )
{
    Integer size = 0 ;
    if ( self != NULL ) size = self->size ;
    return size ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Stride along a dimension.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer View2D_GetStride ( const View2D *self, const Integer dimension, Status *status )
{
    Integer stride = 0 ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
             if ( dimension == 0 ) stride = self->stride0 ;
        else if ( dimension == 1 ) stride = self->stride1 ;
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
    return stride ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization of a compact view given a shape.
!---------------------------------------------------------------------------------------------------------------------------------*/
void View2D_Initialize ( View2D *self, const Integer rows, const Integer columns, Status *status )
{
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        if ( ( rows < 0 ) || ( columns < 0 ) ) Status_Set ( status, Status_InvalidArgument ) ;
        self->extent0 = Maximum ( rows   , 0 ) ;
        self->extent1 = Maximum ( columns, 0 ) ;
        self->offset  = 0 ;
        self->size    = self->extent0 * self->extent1 ;
        self->stride0 = self->extent1 ;
        self->stride1 = 1 ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the default iterator (flattened as much as possible).
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Why not 1D here as well? */
Iterator *View2D_MakeIterator ( const View2D *self, Status *status )
{
    Iterator *iterator = NULL ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto IteratorND *iteratorND ;
        if ( View2D_IsUniform ( self ) )
        {
            iteratorND = IteratorND_Allocate ( 1, status ) ;
            IteratorND_Initialize            ( iteratorND, self->offset, &(self->size), &(self->stride1), status ) ;
        }
        else
        {
            auto Integer extents[2], strides[2] ;
            extents[0] = self->extent0 ; extents[1] = self->extent1 ;
            strides[0] = self->stride0 ; strides[1] = self->stride1 ;
            iteratorND = IteratorND_Allocate ( 2, status ) ;
            IteratorND_Initialize            ( iteratorND, self->offset, extents, strides, status ) ;
        }
        iterator = Iterator_Allocate ( status ) ;
        IteratorND_MakeIterator ( iteratorND, iterator ) ;
    }
    return iterator ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a row iterator given a selection.
!---------------------------------------------------------------------------------------------------------------------------------*/
Iterator *View2D_MakeRowIterator ( const View2D *self, const Selection *rows, Status *status )
{
    Iterator *iterator = NULL ;
    if ( ( self != NULL ) && ( rows != NULL ) && Status_IsOK ( status ) )
    {
        if ( Selection_UpperBound ( rows ) <= self->extent0 )
        {
            auto Integer        extent0 ;
            auto RowIterator2D *iterator2D ;
            extent0    = Selection_Capacity ( rows ) ;
            iterator2D = RowIterator2D_Allocate ( extent0, status ) ;
            RowIterator2D_Initialize ( iterator2D, self->extent1 ,
                                                   self->offset  ,
                                                   self->stride0 ,
                                                   self->stride1 ,
                                                   rows->indices ,
                                                   status        ) ;
            iterator = Iterator_Allocate ( status ) ;
            RowIterator2D_MakeIterator ( iterator2D, iterator ) ;
        }
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
    return iterator ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Row view.
!---------------------------------------------------------------------------------------------------------------------------------*/
void View2D_RowView ( const View2D  *self   ,
                      const Integer  i      ,
                            View1D  *view   ,
                            Status  *status )
{
    if ( ( self != NULL ) && ( view != NULL ) && Status_IsOK ( status ) )
    {
        if ( Integer_IsInRange ( i, 0, self->extent0 ) )
        {
            view->extent  = self->extent1 ;
            view->offset  = self->offset + i * self->stride0 ;
            view->size    = self->extent1 ;
            view->stride  = self->stride1 ;
        }
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set the state of the view - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
void View2D_SetState (       View2D  *self    ,
                       const Integer  extent0 ,
                       const Integer  extent1 ,
                       const Integer  offset  ,
                       const Integer  size    ,
                       const Integer  stride0 ,
                       const Integer  stride1 )
{
    if ( self != NULL )
    {
        self->extent0 = extent0 ;
        self->extent1 = extent1 ;
        self->offset  = offset  ;
        self->size    = size    ;
        self->stride0 = stride0 ;
        self->stride1 = stride1 ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . View.
!---------------------------------------------------------------------------------------------------------------------------------*/
void View2D_View ( const View2D  *self    ,
                   const Integer  start0  ,
                   const Integer  start1  ,
                   const Integer  extent0 ,
                   const Integer  extent1 ,
                   const Integer  stride0 ,
                   const Integer  stride1 ,
                         View2D  *view    ,
                         Status  *status  )
{
    if ( ( self != NULL ) && ( view != NULL ) && Status_IsOK ( status ) )
    {
        view->extent0 = extent0 ;
        view->extent1 = extent1 ;
        view->offset  = self->offset + View2D_ItemIndex ( self, start0, start1 ) ;
        view->size    = extent0 * extent1 ;
        view->stride0 = self->stride0 * stride0 ;
        view->stride1 = self->stride1 * stride1 ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . View 1-D along dimension.
!---------------------------------------------------------------------------------------------------------------------------------*/
void View2D_View1D ( const View2D  *self      ,
                     const Integer  dimension ,
                     const Integer  start0    ,
                     const Integer  start1    ,
                     const Integer  extent    ,
                     const Integer  stride    ,
                           View1D  *view      ,
                           Status  *status    )
{
    if ( ( self != NULL ) && ( view != NULL ) && Status_IsOK ( status ) )
    {
        view->extent = extent ;
        view->offset = self->offset + View2D_ItemIndex ( self, start0, start1 ) ;
        view->size   = extent ;
        if ( dimension == 0 ) view->stride = self->stride0 * stride ;
        else                  view->stride = self->stride1 * stride ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a 1-D view from a multislice.
!---------------------------------------------------------------------------------------------------------------------------------*/
void View2D_View1DMultiSlice ( const View2D *self, const MultiSlice *multiSlice, View1D *view, Status *status )
{
    if ( ( self != NULL ) && ( multiSlice != NULL ) && ( view != NULL ) && Status_IsOK ( status ) )
    {
        if ( ( multiSlice->capacity == 2 ) && ( multiSlice->rank == 1 ) )
        {
            auto Integer rank = 0 ;
            view->offset = self->offset + ( multiSlice->items[0].start * self->stride0 ) + 
                                          ( multiSlice->items[1].start * self->stride1 ) ;
            if ( ! multiSlice->items[0].isScalar )
            {
                rank += 1 ;
                view->extent = multiSlice->items[0].extent ;
                view->stride = multiSlice->items[0].stride * self->stride0 ;
            }
            if ( ! multiSlice->items[1].isScalar )
            {
                rank += 1 ;
                view->extent = multiSlice->items[1].extent ;
                view->stride = multiSlice->items[1].stride * self->stride1 ;
            }
            view->size = view->extent ;
            if ( rank != 1 ) Status_Set ( status, Status_NonConformableArrays ) ;
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}


/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a 2-D view from a multislice.
!---------------------------------------------------------------------------------------------------------------------------------*/
void View2D_ViewMultiSlice ( const View2D *self, const MultiSlice *multiSlice, View2D *view, Status *status )
{
    if ( ( self != NULL ) && ( multiSlice != NULL ) && ( view != NULL ) && Status_IsOK ( status ) )
    {
        if ( multiSlice->capacity == 2 )
        {
            view->offset  = self->offset + ( multiSlice->items[0].start * self->stride0 ) + 
                                           ( multiSlice->items[1].start * self->stride1 ) ;
            view->extent0 = multiSlice->items[0].extent ;
            view->extent1 = multiSlice->items[1].extent ;
            view->stride0 = multiSlice->items[0].stride * self->stride0 ;
            view->stride1 = multiSlice->items[1].stride * self->stride1 ;
            view->size    = ( view->extent0 * view->extent1 ) ;
            if ( multiSlice->items[0].isScalar || multiSlice->items[1].isScalar ) Status_Set ( status, Status_NonConformableArrays ) ;
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}
