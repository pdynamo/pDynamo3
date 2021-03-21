/*==================================================================================================================================
! . Basic 1-D array functions.
!=================================================================================================================================*/

/* . Need _ArrayDataFormat, _ArrayDataPerLine, _ArrayDataType, _ArrayDataTypeInitializer, _UseCBLAS. */

# include <stdlib.h>
# include <stdio.h>

# include "Array_Macros.h"
# include "Memory.h"
# include "Status.h"
# include "TemplateMacros.h"
# ifndef _NoNumeric
# include "TypedMemoryBlock_Macros.h"
# endif
# include "Utilities1D.h"

# define _ArrayType Array1D
# define _BlockType _MakeToken ( _ArrayDataType,  Block     )
# define _ArrayName _MakeToken ( _ArrayDataType, _ArrayType )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local sort function.
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifndef _NoNumeric
static Integer _ItemCompare ( const void *vTerm1, const void *vTerm2 )
{
    Integer i = 0 ;
    _ArrayDataType *term1 = ( _ArrayDataType * ) vTerm1 ;
    _ArrayDataType *term2 = ( _ArrayDataType * ) vTerm2 ;
         if ( term1[0] < term2[0] ) i = -1 ;
    else if ( term1[0] > term2[0] ) i =  1 ;
    return i ;
}
# endif

/*----------------------------------------------------------------------------------------------------------------------------------
! . Absolute maximum.
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifndef _NoNumeric
_ArrayDataType _MakeToken ( _ArrayName, _AbsoluteMaximum ) ( const _ArrayName *self )
{
    _ArrayDataType value = _ArrayDataTypeInitializer ;
    if ( self != NULL ) Utilities1D_AbsoluteMaximum ( self->extent, self->data, self->stride, value ) ;
    return value ;
}
# endif

/*----------------------------------------------------------------------------------------------------------------------------------
! . Absolute maximum index.
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifndef _NoNumeric
Integer _MakeToken ( _ArrayName, _AbsoluteMaximumIndex ) ( const _ArrayName *self )
{
    Integer index = -1 ;
    if ( ( self != NULL ) && ( self->extent > 0 ) )
    {
# ifdef _UseCBLAS
        index = cblas_idamax ( self->extent, Array_DataPointer ( self ), self->stride ) ;
# else
        auto  Integer i ;
        auto _ArrayDataType maximum = abs ( Array1D_Item ( self, 0 ) ), next ;
        index = 0 ;
        for ( i = 1 ; i < self->extent ; i++ )
        {
            next = abs ( Array1D_Item ( self, i ) ) ;
            if ( next > maximum ) { maximum = next ; index = i ; }
        }
# endif
    }
    return index ;
}
# endif

/*----------------------------------------------------------------------------------------------------------------------------------
! . Adding.
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifndef _NoNumeric
void  _MakeToken ( _ArrayName, _Add ) (       _ArrayName     *self   ,
                                        const _ArrayDataType  value  ,
                                        const _ArrayName     *other  ,
                                               Status        *status )
{
    if ( ( self != NULL ) && ( other != NULL )  && Status_IsOK ( status ) )
    {
        Utilities1D_Add ( self->extent, self->data, self->stride, other->extent, other->data, other->stride, value, status ) ;
    }
}
# endif

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
_ArrayName * _MakeToken ( _ArrayName, _Allocate ) ( Status *status )
{
    _ArrayName *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( _ArrayName ) ;
        if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
        else
        {
            self->block = NULL ;
            self->data  = NULL ;
            View1D_InitializeFields ( self );
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation with extent.
!---------------------------------------------------------------------------------------------------------------------------------*/
_ArrayName * _MakeToken ( _ArrayName, _AllocateWithExtent ) ( const Integer extent, Status *status )
{
    _ArrayName *self = _MakeToken ( _ArrayName, _Allocate ) ( status ) ;
    View1D_Initialize ( ( View1D * ) self, extent, status ) ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        _BlockType *block = _MakeToken ( _BlockType, _Allocate ) ( self->size, status ) ;
        if ( block != NULL ) Array_AssignBlockWithReference ( self, block, 0 ) ;
    }
    if ( ! Status_IsOK ( status ) ) _MakeToken ( _ArrayName, _Deallocate ) ( &self ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Assign a block to the array.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _AssignBlock ) (       _ArrayName     *self          ,
                                                     _BlockType     *block         ,
                                               const  Boolean        withReference ,
                                                      Status        *status        )
{
    if ( ( self != NULL ) && ( block != NULL ) && Status_IsOK ( status ) )
    {
        if ( withReference ) { Array_AssignBlockWithReference    ( self, block, self->offset ) ; }
        else                 { Array_AssignBlockWithoutReference ( self, block, self->offset ) ; }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deep cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
_ArrayName * _MakeToken ( _ArrayName, _CloneDeep ) ( const _ArrayName *self, Status *status )
{
    _ArrayName *clone = NULL ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        clone = _MakeToken ( _ArrayName, _AllocateWithExtent ) ( self->extent, status ) ;
        _MakeToken ( _ArrayName, _CopyTo ) ( self, clone, status ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Shallow cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
_ArrayName * _MakeToken ( _ArrayName, _CloneShallow ) ( const _ArrayName *self, Status *status )
{
    _ArrayName *clone = NULL ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        clone = _MakeToken ( _ArrayName, _Allocate ) ( status ) ;
        View1D_CopyTo ( ( View1D * ) self, ( View1D * ) clone ) ;
        _MakeToken ( _ArrayName, _AssignBlock ) ( clone, self->block, True, status ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copying.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _CopyTo ) ( const _ArrayName *self, _ArrayName *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) && Status_IsOK ( status ) )
    {
        Utilities1D_CopyTo ( self->extent, self->data, self->stride, other->extent, other->data, other->stride, status ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _Deallocate ) ( _ArrayName **self )
{
    if ( (*self) != NULL )
    {
        (*self)->data = NULL ;
        _MakeToken ( _BlockType, _Deallocate ) ( &((*self)->block) ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Dot product.
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifndef _NoNumeric
_ArrayDataType _MakeToken ( _ArrayName, _Dot ) ( const _ArrayName *self   ,
                                                 const _ArrayName *other  ,
                                                        Status    *status )
{
    _ArrayDataType value = _ArrayDataTypeInitializer ;
    if ( ( self != NULL ) && ( other != NULL )  && Status_IsOK ( status ) )
    {
        Utilities1D_Dot ( self->extent, self->data, self->stride, other->extent, other->data, other->stride, value, status ) ;
    }
    return value ;
}
# endif

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
_ArrayDataType _MakeToken ( _ArrayName, _GetItem ) ( const _ArrayName *self, const Integer i, Status *status )
{
    _ArrayDataType value = _ArrayDataTypeInitializer ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        if ( ( Integer_IsInRange ( i, 0, self->extent ) ) ) value = Array1D_Item ( self, i ) ;
        else                                                Status_Set ( status, Status_IndexOutOfRange ) ;
    }
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Incrementing.
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifndef _NoNumeric
void _MakeToken ( _ArrayName, _Increment ) ( _ArrayName *self, const _ArrayDataType value )
{
    if ( self != NULL ) Utilities1D_Increment ( self->extent, self->data, self->stride, value ) ;
}
# endif

/*----------------------------------------------------------------------------------------------------------------------------------
! . Left circular shift.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _LeftCircularShift ) ( _ArrayName *self )
{
    if ( ( self != NULL ) && ( self->extent > 1 ) )
    {
        auto  Integer       i ;
        auto _ArrayDataType first = Array1D_Item ( self, 0 ) ;
        for ( i = 1 ; i < self->extent ; i++ ) Array1D_Item ( self, i-1 ) = Array1D_Item ( self, i ) ;
        Array1D_Item ( self, self->extent-1 ) = first ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Maximum.
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifndef _NoNumeric
_ArrayDataType _MakeToken ( _ArrayName, _Maximum ) ( const _ArrayName *self )
{
    _ArrayDataType value = _ArrayDataTypeInitializer ;
    if ( self != NULL ) Utilities1D_Maximum ( self->extent, self->data, self->stride, value ) ;
    return value ;
}
# endif

/*----------------------------------------------------------------------------------------------------------------------------------
! . Minimum.
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifndef _NoNumeric
_ArrayDataType _MakeToken ( _ArrayName, _Minimum ) ( const _ArrayName *self )
{
    _ArrayDataType value = _ArrayDataTypeInitializer ;
    if ( self != NULL ) Utilities1D_Minimum ( self->extent, self->data, self->stride, value ) ;
    return value ;
}
# endif

/*----------------------------------------------------------------------------------------------------------------------------------
! . Multiplying.
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifndef _NoNumeric
void _MakeToken ( _ArrayName, _Multiply ) ( _ArrayName *self, const _ArrayName *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) && Status_IsOK ( status ) )
    {
        Utilities1D_Multiply ( self->extent, self->data, self->stride, other->extent, other->data, other->stride, status ) ;
    }
}
# endif

/*----------------------------------------------------------------------------------------------------------------------------------
! . Pointer to data.
!---------------------------------------------------------------------------------------------------------------------------------*/
_ArrayDataType * _MakeToken ( _ArrayName, _PointerToData ) ( const _ArrayName *self ) { return Array_DataPointer ( self ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Printing.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _Print ) ( const _ArrayName  *self )
{
    if ( ( self == NULL ) || ( self->extent <= 0 ) ) printf ( "Empty 1-D array.\n" ) ;
    else
    {
        auto Integer i, i0, i1, l, n ;
        auto div_t q = div ( self->extent, _ArrayDataPerLine ) ;
        n = q.quot ; if ( q.rem > 0 ) n += 1 ;
        for ( l = 0 ; l < n ; l++ )
        {
            i0 = ( l * _ArrayDataPerLine ) ;
            i1 = Minimum ( i0 + _ArrayDataPerLine, self->extent ) ;
            printf ( "%4d--%4d", i0, i1 ) ;
            for ( i = i0 ; i < i1 ; i++ ) printf ( _ArrayDataFormat, Array1D_Item ( self, i ) ) ;
            printf ( "\n" ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Resizing.
! . This is only done if the array owns the block and their are no other views.
! . Care should be taken not to resize when there are unreferenced views.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _Resize ) (       _ArrayName *self   ,
                                          const  Integer    extent ,
                                                 Status    *status )
{
    if ( ( self != NULL ) && ( extent != self->extent ) && Status_IsOK ( status ) )
    {
        auto _BlockType *block = self->block ;
        if ( ( block != NULL ) && ( block->references <= 2 ) && ( self->size == block->capacity ) )
        {
            auto Integer n = Maximum ( extent, 0 ) ;
            _MakeToken ( _BlockType, _Resize ) ( block, n, status ) ;
            if ( Status_IsOK ( status ) )
            {
                self->extent = n ;
                self->size   = n ;
                self->data   = Array_BlockDataPointer ( self, self->offset ) ;
            }
            else
            {
                _MakeToken ( _BlockType, _Deallocate ) ( &block ) ;
                self->block = NULL ;
                self->data  = NULL ;
            }
        }
        else Status_Set ( status, Status_InvalidArrayOperation ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Resizing with initialization of any excess data.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _ResizeWithInitializer ) (       _ArrayName     *self        ,
                                                         const  Integer        extent      ,
                                                         const _ArrayDataType  initializer ,
                                                                Status        *status      )
{
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer oldE = self->extent ;
        _MakeToken ( _ArrayName, _Resize ) ( self, extent, status ) ;
        if ( ( oldE < extent ) && Status_IsOK ( status ) )
        {
            auto _ArrayName view ;
            _MakeToken ( _ArrayName, _View ) ( self, oldE, ( extent - oldE ), 1, False, &view, status ) ;
            _MakeToken ( _ArrayName, _Set  ) ( self, initializer ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . In-place reverse.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _Reverse ) ( _ArrayName *self )
{
    if ( ( self != NULL ) && ( self->extent > 1 ) )
    {
        auto  Integer       i, j ;
        auto _ArrayDataType t ;
        for ( i = 0, j = self->extent - 1 ; i < ( self->extent / 2 ) ; i++, j-- )
        {
            t = Array1D_Item ( self, i ) ;
            Array1D_Item ( self, i ) = Array1D_Item ( self, j ) ;
            Array1D_Item ( self, j ) = t ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Right circular shift.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _RightCircularShift ) ( _ArrayName *self )
{
    if ( ( self != NULL ) && ( self->extent > 1 ) )
    {
        auto  Integer       i ;
        auto _ArrayDataType last = Array1D_Item ( self, self->extent-1 ) ;
        for ( i = self->extent-1 ; i > 0 ; i-- ) Array1D_Item ( self, i ) = Array1D_Item ( self, i-1 ) ;
        Array1D_Item ( self, 0 ) = last ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Scaling.
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifndef _NoNumeric
void _MakeToken ( _ArrayName, _Scale ) ( _ArrayName *self, const _ArrayDataType value )
{
    if ( self != NULL ) Utilities1D_Scale ( self->extent, self->data, self->stride, value ) ;
}
# endif

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set all items.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _Set ) ( _ArrayName *self, const _ArrayDataType value )
{
    if ( self != NULL ) Utilities1D_Set ( self->extent, self->data, self->stride, value ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _SetItem ) ( _ArrayName *self, const Integer i, const _ArrayDataType value, Status *status )
{
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        if ( Integer_IsInRange ( i, 0, self->extent ) ) Array1D_Item ( self, i ) = value ;
        else                                            Status_Set ( status, Status_IndexOutOfRange ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Sort - ascending in-place.
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifndef _NoNumeric
void _MakeToken ( _ArrayName, _Sort ) ( _ArrayName *self )
{
    if ( ( self != NULL ) && ( self->extent > 1 ) )
    {
        qsort ( ( void * ) Array_DataPointer ( self )    ,
                ( size_t ) self->extent                  ,
                self->stride * sizeof ( _ArrayDataType ) ,
                ( void * ) _ItemCompare                  ) ;
    }
}
# endif

/*----------------------------------------------------------------------------------------------------------------------------------
! . Index sort.
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifndef _NoNumeric
void _MakeToken ( _ArrayName, _SortIndex ) ( const _ArrayName *self, IntegerBlock *indices, Status *status )
{
    if ( ( self != NULL ) && ( indices != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer c = Block_Capacity ( indices ) , n = self->extent ;
        if ( c >= n )
        {
            auto Integer j, k, ki, l, m, tmp ;
            for ( k = 0 ; k < n ; k++ ) Block_Item ( indices, k ) =  k ;
            for ( k = n ; k < c ; k++ ) Block_Item ( indices, k ) = -1 ;
            m = n - 1 ;
            k = m / 2;
            k ++ ;
            do
            {
                k -- ;
                ki = Block_Item ( indices, k ) ;
                l  = k ;
                while ( l <= ( m / 2 ) )
                {
                    j = 2 * l ;
                    if ( ( j < m ) && ( Array1D_Item ( self, Block_Item ( indices, j ) ) < Array1D_Item ( self, Block_Item ( indices, j+1 ) ) ) ) j++ ;
                    if ( Array1D_Item ( self, ki ) >= Array1D_Item ( self, Block_Item ( indices, j ) ) ) break ;
                    Block_Item ( indices, l ) = Block_Item ( indices, j );
                    l = j;
                }
                Block_Item ( indices, l ) = ki;
            }
            while ( k > 0 ) ;
            while ( m > 0 )
            {
                tmp = Block_Item ( indices, 0 ) ; Block_Item ( indices, 0 ) = Block_Item ( indices, m ) ; Block_Item ( indices, m ) = tmp ;
                m -- ;
                ki = Block_Item ( indices, 0 ) ;
                l  = k ;
                while ( l <= ( m / 2 ) )
                {
                    j = 2 * l ;
                    if ( j < m && ( Array1D_Item ( self, Block_Item ( indices, j ) ) < Array1D_Item ( self, Block_Item ( indices, j+1 ) ) ) ) j++ ;
                    if ( Array1D_Item ( self, ki ) >= Array1D_Item ( self, Block_Item ( indices, j ) ) ) break ;
                    Block_Item ( indices, l ) = Block_Item ( indices, j );
                    l = j;
                }
                Block_Item ( indices, l ) = ki;
            }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}
# endif

/*----------------------------------------------------------------------------------------------------------------------------------
! . Sort unique.
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifndef _NoNumeric
Integer _MakeToken ( _ArrayName, _SortUnique ) ( _ArrayName *self )
{
    Integer unique = 0 ;
    _MakeToken ( _ArrayName, _Sort ) ( self ) ;
    if ( ( self != NULL ) && ( self->extent > 1 ) )
    {
        auto Integer i, n ;
        for ( i = 1, n = 1 ; i < self->extent ; i++ )
        {
            if ( Array1D_Item ( self, i-1 ) != Array1D_Item ( self, i ) )
            {
                if ( n < i ) Array1D_Item ( self, n ) = Array1D_Item ( self, i ) ;
                n++ ;
            }
        }
        unique = n ;
    }
    return unique ;
}
# endif

/*----------------------------------------------------------------------------------------------------------------------------------
! . Summing.
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifndef _NoNumeric
_ArrayDataType _MakeToken ( _ArrayName, _Sum ) ( const _ArrayName *self )
{
    _ArrayDataType value = _ArrayDataTypeInitializer ;
    if ( self != NULL ) Utilities1D_Sum ( self->extent, self->data, self->stride, value ) ;
    return value ;
}
# endif

/*----------------------------------------------------------------------------------------------------------------------------------
! . View.
! . Extent, start and stride need to be checked before entry.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _View ) ( const _ArrayName *self          ,
                                        const  Integer    start         ,
                                        const  Integer    extent        ,
                                        const  Integer    stride        ,
                                        const  Boolean    withReference ,
                                              _ArrayName *view          ,
                                               Status    *status        )
{
    if ( ( self != NULL ) && ( view != NULL ) && Status_IsOK ( status ) )
    {
        View1D_View ( ( View1D * ) self, start, extent, stride, ( View1D * ) view, status ) ; 
        if ( withReference ) { Array_AssignBlockWithReference    ( view, self->block, view->offset ) ; }
        else                 { Array_AssignBlockWithoutReference ( view, self->block, view->offset ) ; }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . View of raw.
! . Offset, extent and stride need to be checked before entry.
! . For limited, local use only as no block.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _ViewOfRaw ) (       _ArrayName     *self   ,
                                             const  Integer        offset ,
                                             const  Integer        extent ,
                                             const  Integer        stride ,
                                                   _ArrayDataType *data   )
{
    if ( ( self != NULL ) && ( data != NULL ) )
    {
        self->block  = NULL   ;
        self->data   = data   ;
        self->extent = extent ;
        self->offset = offset ;
        self->size   = extent ;
        self->stride = stride ;
    }
}

# undef _ArrayName
# undef _ArrayType
# undef _BlockType
