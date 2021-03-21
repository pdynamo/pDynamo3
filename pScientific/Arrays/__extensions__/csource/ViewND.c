/*==================================================================================================================================
! . N-D array view.
!=================================================================================================================================*/

# include "IntegerUtilities.h"
# include "IteratorND.h"
# include "Memory.h"
# include "NumericalMacros.h"
# include "ViewND.h"

# include <stdio.h>

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation of type.
!---------------------------------------------------------------------------------------------------------------------------------*/
ViewND *ViewND_Allocate ( Status *status )
{
    ViewND *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( ViewND ) ;
        if ( self != NULL )
        {
            self->extents = NULL  ;
            self->strides = NULL  ;
            self->rank    = 0 ;
            self->offset  = 0 ;
            self->size    = 0 ;
        }
        else Status_Set ( status, Status_OutOfMemory ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation with rank (>= 0).
!---------------------------------------------------------------------------------------------------------------------------------*/
ViewND *ViewND_AllocateWithRank ( const Integer rank, Status *status )
{
    ViewND *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        if ( rank >= 0 )
        {
            self = ViewND_Allocate ( status ) ;
            if ( self != NULL )
            {
                self->rank    = rank ;
                if ( rank > 0 )
                {
                    self->extents = Integer_Allocate ( rank, status ) ;
                    self->strides = Integer_Allocate ( rank, status ) ;
                    if ( ( self->extents == NULL ) || ( self->strides == NULL ) ) ViewND_Deallocate ( &self ) ;
                }
            }
            if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
        }
        else Status_Set ( status, Status_InvalidArgument ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation with shape.
!---------------------------------------------------------------------------------------------------------------------------------*/
ViewND *ViewND_AllocateWithShape ( const Integer rank, const Integer *extents, Status *status )
{
    ViewND *self = ViewND_AllocateWithRank ( rank, status ) ;
    ViewND_Initialize ( self, extents, status ) ;
    if ( ! Status_IsOK ( status ) ) ViewND_Deallocate ( &self ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check two views for conformability.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean ViewND_AreConformable ( const ViewND *self, const ViewND *other, Status *status )
{
    Boolean areConformable = False ;
    if ( ( self != NULL ) && ( other != NULL ) && Status_IsOK ( status ) )
    {
        if ( ( self->rank == other->rank ) && ( self->size == other->size ) )
        {
            auto Integer d ;
            areConformable = True ;
            for ( d = 0 ; d < self->rank ; d++ )
            {
                if ( self->extents[d] != other->extents[d] ) { areConformable = False ; break ; }
            }
        }
        if ( ! areConformable ) Status_Set ( status, Status_NonConformableArrays ) ;
    }
    return areConformable ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check the view versus a capacity.
! . Both maximum and minimum indices are calculated and checked.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean ViewND_CheckCapacity ( const ViewND *self, const Integer capacity )
{
    Boolean isOK = False ;
    if ( self != NULL )
    {
        auto Integer d, m = 0, n = 0, p ;
        if ( self->size > 0 )
        {
            m = self->offset ;
            n = self->offset ;
            for ( d = 0 ; d < self->rank ; d++ )
            {
                p = ( self->extents[d] - 1 ) * self->strides[d] ;
                if ( p < 0 ) m += p ;
                else         n += p ;
            }
        }
        isOK = ( ( m >= 0 ) && ( n < capacity ) ) ;
    }
    return isOK ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
ViewND *ViewND_Clone ( const ViewND *self, Status *status )
{
    ViewND *clone = NULL ;
    if ( self != NULL )
    {
        clone = ViewND_AllocateWithShape ( self->rank, self->extents, status ) ;
        if ( clone != NULL )
        {
            auto Integer d ;
            clone->offset = self->offset ;
            for ( d = 0 ; d < self->rank ; d++ ) clone->strides[d] = self->strides[d] ;
        }
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ViewND_Deallocate ( ViewND **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        Memory_Deallocate ( (*self)->extents ) ;
        Memory_Deallocate ( (*self)->strides ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Flattening.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer ViewND_Flatten ( const ViewND *self, Integer *extents, Integer *strides )
{
    Integer rank = 0 ;
    if ( ( self != NULL ) && ( self->rank > 0 ) )
    {
        auto Integer d, e, s ;
        extents[0] = self->extents[0] ;
        strides[0] = self->strides[0] ;
        for ( d = 1 ; d < self->rank ; d++ )
        {
            e = self->extents[d] ;
            s = self->strides[d] ;
            if ( self->strides[d-1] == ( e * s ) ) {             extents[rank] *= e ; }
            else                                   { rank += 1 ; extents[rank]  = e ; }
            strides[rank] = s ;
        }
        rank += 1 ;
    }
    return rank ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Create view from state - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
ViewND *ViewND_FromState ( const Integer  rank   ,
                           const Integer  offset ,
                           const Integer  size   ,
                                 Status  *status )
{
    ViewND *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        self = ViewND_AllocateWithRank ( rank, status ) ;
        if ( self != NULL )
        {
            self->offset = offset ;
            self->size   = size   ;
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Extent along a dimension.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer ViewND_GetExtent ( const ViewND *self, const Integer dimension, Status *status )
{
    Integer extent = 0 ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        if ( ( dimension >= 0 ) && ( dimension < self->rank ) ) extent = self->extents[dimension] ; 
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
    return extent ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the index of an item.
! . The input rank can be less than the view rank in which case zero is assumed for the remaining indices.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer ViewND_GetIndex ( const ViewND *self, const Integer rank, const Integer *indices, Status *status )
{
    Integer index = -1 ;
    if ( ( self != NULL ) && ( rank >= 0 ) && ( indices != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer d, i ;
        index = 0 ; /*self->offset ;*/
        for ( d = 0 ; d < Minimum ( rank, self->rank ) ; d++ )
        {
            i = indices[d] ;
            if ( ( i < 0 ) || ( i >= self->extents[d] ) )
            {
                Status_Set ( status, Status_IndexOutOfRange ) ;
                index = -1 ;
                break ;
            }
            index += i * self->strides[d] ;
        }
    }
    return index ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the index of an item given a scalar multislice (no checking).
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer ViewND_GetIndexMultiSlice ( const ViewND *self, const MultiSlice *multiSlice, Status *status )
{
    Integer index = -1 ;
    if ( ( self != NULL ) && ( multiSlice != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer d, i ;
        index = 0 ; /*self->offset ;*/
        for ( d = 0 ; d < Minimum ( self->rank, multiSlice->capacity ) ; d++ )
        {
            i = multiSlice->items[d].start ;
            if ( ( i < 0 ) || ( i >= self->extents[d] ) )
            {
                Status_Set ( status, Status_IndexOutOfRange ) ;
                index = -1 ;
                break ;
            }
            index += i * self->strides[d] ;
        }
    }
    return index ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Offset.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer ViewND_GetOffset ( const ViewND *self ) { return ( ( (self) == NULL ) ? 0 : (self)->offset ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Rank.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer ViewND_GetRank ( const ViewND *self )
{
    Integer rank = 0 ;
    if ( self != NULL ) rank = self->rank ;
    return rank ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Size.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer ViewND_GetSize ( const ViewND *self )
{
    Integer size = 0 ;
    if ( self != NULL ) size = self->size ;
    return size ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Stride along a dimension.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer ViewND_GetStride ( const ViewND *self, const Integer dimension, Status *status )
{
    Integer stride = 0 ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        if ( ( dimension >= 0 ) && ( dimension < self->rank ) ) stride = self->strides[dimension] ; 
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
    return stride ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization of a compact view given a shape.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ViewND_Initialize ( ViewND *self, const Integer *extents, Status *status )
{
    if ( self != NULL )
    {
        auto Integer d, extent, n ;
        for ( d = self->rank - 1, n = 1 ; d >= 0 ; d-- )
        {
            extent = extents[d] ;
            if ( extent < 0 )
            {
                extent = 0 ;
                Status_Set ( status, Status_InvalidArgument ) ;
            }
            self->extents[d] = extent ;
            self->strides[d] = n ;
            n *= extent ;
        }
        self->size = n ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for compactness (minimal spacing of items).
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean ViewND_IsCompact ( const ViewND *self )
{
    Boolean isCompact = False ;
    if ( ( self != NULL ) && ( self->rank > 0 ) ) isCompact = ViewND_IsUniform ( self ) && ( abs ( self->strides[self->rank-1] ) == 1 ) ;
    return isCompact ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for uniformity (equal spacing of items).
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean ViewND_IsUniform ( const ViewND *self )
{
    Boolean isUniform = False ;
    if ( self != NULL )
    {
        auto Integer d ;
        isUniform = True ;
        for ( d = self->rank - 1 ; d >= 1 ; d-- )
        {
            if ( self->strides[d-1] != self->extents[d] * self->strides[d] ) { isUniform = False ; break ; }
        }
    }
    return isUniform ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the default iterator (flattened as much as possible).
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Why not 1D here as well? */
Iterator *ViewND_MakeIterator ( const ViewND *self, Status *status )
{
    Iterator *iterator = NULL ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer     rank, extents[self->rank], strides[self->rank] ;
        auto IteratorND *iteratorND ; 
        rank       = ViewND_Flatten      ( self, extents, strides ) ;
        iterator   = Iterator_Allocate   (       status ) ;
        iteratorND = IteratorND_Allocate ( rank, status ) ;
        IteratorND_Initialize            ( iteratorND, self->offset, extents, strides, status ) ;
        IteratorND_MakeIterator          ( iteratorND, iterator ) ;
    }
    return iterator ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set the extent and stride along a dimension - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ViewND_SetExtentStride (       ViewND *self      ,
                              const Integer dimension ,
                              const Integer extent    ,
                              const Integer stride    ,
                                    Status *status    )
{
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        if ( ( dimension >= 0 ) && ( dimension < self->rank ) )
        {
            self->extents[dimension] = extent ;
            self->strides[dimension] = stride ;
        }
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a 1-D view from a multislice.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ViewND_View1DMultiSlice ( const ViewND *self, const MultiSlice *multiSlice, View1D *view, Status *status )
{
    if ( ( self != NULL ) && ( multiSlice != NULL ) && ( view != NULL ) && Status_IsOK ( status ) )
    {
        if ( ( self->rank == multiSlice->capacity ) && ( multiSlice->rank == 1 ) )
        {
            auto Integer d ;
            view->offset = self->offset ;
            for ( d = 0 ; d < self->rank ; d++ )
            {
                view->offset += ( multiSlice->items[d].start * self->strides[d] ) ;
                if ( ! multiSlice->items[d].isScalar )
                {
                    view->extent = multiSlice->items[d].extent ;
                    view->stride = multiSlice->items[d].stride * self->strides[d] ;
                }
            }
            view->size = view->extent ;
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a 2-D view from a multislice.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ViewND_View2DMultiSlice ( const ViewND *self, const MultiSlice *multiSlice, View2D *view, Status *status )
{
    if ( ( self != NULL ) && ( multiSlice != NULL ) && ( view != NULL ) && Status_IsOK ( status ) )
    {
        if ( ( self->rank == multiSlice->capacity ) && ( multiSlice->rank == 2 ) )
        {
            auto Integer c, d ;
            view->offset = self->offset ;
            for ( c = d = 0 ; d < self->rank ; d++ )
            {
                view->offset += ( multiSlice->items[d].start * self->strides[d] ) ;
                if ( ! multiSlice->items[d].isScalar )
                {
                    if ( c == 0 )
                    {
                        view->extent0 = multiSlice->items[d].extent ;
                        view->stride0 = multiSlice->items[d].stride * self->strides[d] ;
                    }
                    else
                    {
                        view->extent1 = multiSlice->items[d].extent ;
                        view->stride1 = multiSlice->items[d].stride * self->strides[d] ;
                    }
                    c += 1 ;
                }
            }
            view->size = view->extent0 * view->extent1 ;
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a view from a multislice.
! . No checking is done and the multislice must be fully defined (i.e. its capacity equals self's rank).
! . For cases where this is not so, the higher dimensions of the multislice should be initialized to those of self.
!---------------------------------------------------------------------------------------------------------------------------------*/
ViewND *ViewND_ViewMultiSlice ( const ViewND *self, const MultiSlice *multiSlice, Status *status )
{
    ViewND *view = NULL ;
    if ( ( self != NULL ) && ( multiSlice != NULL ) && Status_IsOK ( status ) )
    {
        if ( self->rank == multiSlice->capacity )
        {
            view = ViewND_AllocateWithRank ( multiSlice->rank, status ) ;
            if ( view != NULL )
            {
                auto Integer d ;
                view->offset = self->offset ;
                view->rank   = 0 ;
                view->size   = 1 ;
                for ( d = 0 ; d < self->rank ; d++ )
                {
                    view->offset += ( multiSlice->items[d].start * self->strides[d] ) ;
                    if ( ! multiSlice->items[d].isScalar )
                    {
                        view->extents[view->rank] = multiSlice->items[d].extent ;
                        view->strides[view->rank] = multiSlice->items[d].stride * self->strides[d] ;
                        view->rank               += 1 ;
                        view->size               *= multiSlice->items[d].extent ;
                    }
                }
            }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
    return view ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get a 2-D view of the last two dimensions of given index.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ViewND_ViewTail2D ( const ViewND *self, const Integer *indices, View2D *view, Status *status )
{
    if ( ( self != NULL ) && ( indices != NULL ) && ( view != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer n = ViewND_GetIndex ( self, self->rank-2, indices, status ) ;
        if ( n >= 0 )
        {
            auto Integer a, b ;
            a = self->rank - 2 ;
            b = self->rank - 1 ;
            view->offset  = self->offset + n ;
            view->size    = ( self->extents[a] * self->extents[b] ) ;
            view->extent0 = self->extents[a] ;
            view->extent1 = self->extents[b] ;
            view->stride0 = self->strides[a] ;
            view->stride1 = self->strides[b] ;
        }
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
}
