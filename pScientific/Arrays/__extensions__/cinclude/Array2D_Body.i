/*==================================================================================================================================
! . Basic 2-D array functions.
!=================================================================================================================================*/

/* . Need _ArrayDataFormat, _ArrayDataPerLine, _ArrayDataType, _ArrayDataTypeInitializer, _UseCBLAS. */

# include <stdlib.h>
# include <stdio.h>

# include "Array_Macros.h"
# include "Integer.h"
# include "Memory.h"
# include "Status.h"
# include "TemplateMacros.h"
# include "Utilities1D.h"

# define _Array1DType Array1D
# define _ArrayType   Array2D
# define _Array1DName _MakeToken ( _ArrayDataType, _Array1DType )
# define _ArrayName   _MakeToken ( _ArrayDataType, _ArrayType   )
# define _BlockType   _MakeToken ( _ArrayDataType,  Block       )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Absolute maximum.
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifndef _NoNumeric
_ArrayDataType _MakeToken ( _ArrayName, _AbsoluteMaximum ) ( const _ArrayName *self )
{
    _ArrayDataType value = _ArrayDataTypeInitializer ;
    if ( self != NULL )
    {
        if ( View2D_IsUniform ( self ) ) { Utilities1D_AbsoluteMaximum ( self->size, self->data, self->stride1, value ) ; }
        else
        {
            auto Integer r ;
            for ( r = 0 ; r < View2D_Rows ( self ) ; r++ )
            {
                Utilities1D_AbsoluteMaximum ( self->extent1, Array2D_RowPointer ( self, r ), self->stride1, value ) ;
            }
        }
    }
    return value ;
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
        if ( View2D_AreConformable ( self, other ) )
        {
# ifdef _UseCBLAS
            if ( View2D_IsUniform ( self ) && View2D_IsUniform ( other ) )
            {
                cblas_daxpy ( View2D_Size ( self )        ,
                              value                       ,
                              Array_DataPointer ( other ) ,
                              other->stride1              ,
                              Array_DataPointer ( self  ) ,
                              self->stride1               ) ;
            }
            else
            {
                auto Integer i ;
                for ( i = 0 ; i < View2D_Rows ( self ) ; i++ )
                {
                    cblas_daxpy ( View2D_Columns ( self )         ,
                                  value                           ,
                                  Array2D_RowPointer ( other, i ) ,
                                  other->stride1                  ,
                                  Array2D_RowPointer ( self , i ) ,
                                  self->stride1                   ) ;
                }
            }
# else
            auto Integer i, j ;
            for ( i = 0 ; i < View2D_Rows ( self ) ; i++ )
            {
                for ( j = 0 ; j < View2D_Columns ( self ) ; j++ ) Array2D_Item ( other, i, j ) += ( value * Array2D_Item ( self, i, j ) ) ;
            }
# endif
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
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
            View2D_InitializeFields ( self );
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation with extent.
!---------------------------------------------------------------------------------------------------------------------------------*/
_ArrayName * _MakeToken ( _ArrayName, _AllocateWithExtents ) ( const Integer rows, const Integer columns, Status *status )
{
    _ArrayName *self = _MakeToken ( _ArrayName, _Allocate ) ( status ) ;
    View2D_Initialize ( ( View2D * ) self, rows, columns, status ) ;
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
        clone = _MakeToken ( _ArrayName, _AllocateWithExtents ) ( self->extent0, self->extent1, status ) ;
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
        View2D_CopyTo ( ( View2D * ) self, ( View2D * ) clone ) ;
        _MakeToken ( _ArrayName, _AssignBlock ) ( clone, self->block, True, status ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Column view.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _ColumnView ) ( const _ArrayName   *self          ,
                                              const  Integer      i             ,
                                              const  Boolean      withReference ,
                                                    _Array1DName *view          ,
                                                     Status      *status        )
{
    if ( ( self != NULL ) && ( view != NULL ) && Status_IsOK ( status ) )
    {
        View2D_ColumnView ( ( View2D * ) self, i, ( View1D * ) view, status ) ;
        if ( withReference ) { Array_AssignBlockWithReference    ( view, self->block, view->offset ) ; }
        else                 { Array_AssignBlockWithoutReference ( view, self->block, view->offset ) ; }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copying.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _CopyTo ) ( const _ArrayName *self, _ArrayName *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) && Status_IsOK ( status ) )
    {
        if ( View2D_AreConformable ( self, other ) )
        {
# ifdef _UseCBLAS
            if ( View2D_IsUniform ( self ) && View2D_IsUniform ( other ) )
            {
                cblas_dcopy ( View2D_Size ( self )        ,
                              Array_DataPointer ( self  ) ,
                              self->stride1               ,
                              Array_DataPointer ( other ) ,
                              other->stride1              ) ;
            }
            else
            {
                auto Integer i ;
                for ( i = 0 ; i < View2D_Rows ( self ) ; i++ )
                {
                    cblas_dcopy ( View2D_Columns ( self )         ,
                                  Array2D_RowPointer ( self , i ) ,
                                  self->stride1                   ,
                                  Array2D_RowPointer ( other, i ) ,
                                  other->stride1                  ) ;
                }
            }
# else
            auto Integer i, j ;
            for ( i = 0 ; i < View2D_Rows ( self ) ; i++ )
            {
                for ( j = 0 ; j < View2D_Columns ( self ) ; j++ ) Array2D_Item ( other, i, j ) = Array2D_Item ( self, i, j ) ;
            }
# endif
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
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
! . Get an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
_ArrayDataType _MakeToken ( _ArrayName, _GetItem ) ( const _ArrayName *self, const Integer i, const Integer j, Status *status )
{
    _ArrayDataType value = _ArrayDataTypeInitializer ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        if ( Integer_IsInRange ( i, 0, self->extent0 ) &&
             Integer_IsInRange ( j, 0, self->extent1 ) ) value = Array2D_Item ( self, i, j ) ;
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get an item from a multislice.
!---------------------------------------------------------------------------------------------------------------------------------*/
_ArrayDataType _MakeToken ( _ArrayName, _GetItemMultiSlice ) ( const _ArrayName  *self       ,
                                                               const  MultiSlice *multiSlice ,
                                                                      Status     *status     )
{
    _ArrayDataType value = _ArrayDataTypeInitializer ;
    if ( ( multiSlice != NULL ) && ( multiSlice->capacity >= 2 ) )
    {
        value = _MakeToken ( _ArrayName, _GetItem ) ( self, multiSlice->items[0].start, multiSlice->items[1].start, status ) ;
    }
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Pointer to data.
!---------------------------------------------------------------------------------------------------------------------------------*/
_ArrayDataType * _MakeToken ( _ArrayName, _PointerToData ) ( const _ArrayName *self ) { return Array_DataPointer ( self ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Printing.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _Print ) ( const _ArrayName  *self )
{
    if ( ( self == NULL ) || ( self->size <= 0 ) ) printf ( "Empty 2-D array.\n" ) ;
    else
    {
        auto Integer i, i0, i1, l, n, r ;
        auto div_t q = div ( self->extent1, _ArrayDataPerLine ) ;
        n = q.quot ; if ( q.rem > 0 ) n += 1 ;
        for ( l = 0 ; l < n ; l++ )
        {
            i0 = ( l * _ArrayDataPerLine ) ;
            i1 = Minimum ( i0 + _ArrayDataPerLine, self->extent1 ) ;
            printf ( "Columns %4d--%4d:\n", i0, i1 ) ;
            for ( r = 0 ; r < self->extent0 ; r++ )
            {
                printf ( "%4d  ", r ) ;
                for ( i = i0 ; i < i1 ; i++ ) printf ( _ArrayDataFormat, Array2D_Item ( self, r, i ) ) ;
                printf ( "\n" ) ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Prune by row.
!---------------------------------------------------------------------------------------------------------------------------------*/
_ArrayName * _MakeToken ( _ArrayName, _PruneByRow ) ( const _ArrayName *self   ,
                                                      const  Selection *toKeep ,
                                                             Status    *status )
{
    _ArrayName *new = NULL ;
    if ( ( self != NULL ) && ( toKeep != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer c, i, m, n ;
        c = Selection_Capacity ( toKeep ) ;
        for ( i = 0, n = 0 ; i < c ; i++ )
        {
            if ( Integer_IsInRange ( Selection_Item ( toKeep, i ), 0, View2D_Rows ( self ) ) ) n += 1 ;
        }
        /* . Create the matrix. */
        new = _MakeToken ( _ArrayName, _AllocateWithExtents) ( n, View2D_Columns ( self ), status ) ;
	if ( new != NULL )
	{
# ifndef _UseCBLAS
            auto Integer j ;
# endif
            for ( i = 0, n = 0 ; i < c ; i++ )
            {
        	m = Selection_Item ( toKeep, i ) ;
        	if ( Integer_IsInRange ( m, 0, self->extent0 ) )
        	{
# ifdef _UseCBLAS
                    cblas_dcopy ( self->extent1, Array2D_RowPointer ( self, m ), self->stride1, Array2D_RowPointer ( new, n ), new->stride1 ) ;
# else
                    for ( j = 0 ; j < self->extent1 ; j++ ) Array2D_Item ( new, n, j ) = Array2D_Item ( self, m, j ) ;
# endif
                    n += 1 ;
        	}
            }
	}
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Resizing.
! . This is only done if the array owns the block and their are no other views.
! . Care should be taken not to resize when there are unreferenced views.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _Resize ) (       _ArrayName *self    ,
                                          const  Integer    extent0 ,
                                                 Status    *status  )
{
    if ( ( self != NULL ) && ( extent0 != self->extent0 ) && Status_IsOK ( status ) )
    {
        auto _BlockType *block = self->block ;
        if ( ( block != NULL ) && ( block->references <= 2 ) && ( self->size == block->capacity ) )
        {
            auto Integer n0 = Maximum ( extent0, 0 ), nS = n0 * self->extent1 ;
            _MakeToken ( _BlockType, _Resize ) ( block, nS, status ) ;
            if ( Status_IsOK ( status ) )
            {
                self->extent0 = n0 ;
                self->size    = nS ;
                self->data    = Array_BlockDataPointer ( self, self->offset ) ;
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
                                                         const  Integer        extent0     ,
                                                         const _ArrayDataType  initializer ,
                                                                Status        *status      )
{
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer oldE = self->extent0 ;
        _MakeToken ( _ArrayName, _Resize ) ( self, extent0, status ) ;
        if ( ( oldE < extent0 ) && Status_IsOK ( status ) )
        {
            auto _ArrayName view ;
            _MakeToken ( _ArrayName, _View ) ( self, oldE, 0, ( extent0 - oldE ), self->extent1, 1, 1, False, &view, status ) ;
            _MakeToken ( _ArrayName, _Set  ) ( self, initializer ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Row view.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _RowView ) ( const _ArrayName   *self          ,
                                           const  Integer      i             ,
                                           const  Boolean      withReference ,
                                                 _Array1DName *view          ,
                                                  Status      *status        )
{
    if ( ( self != NULL ) && ( view != NULL ) && Status_IsOK ( status ) )
    {
        View2D_RowView ( ( View2D * ) self, i, ( View1D * ) view, status ) ;
        if ( withReference ) { Array_AssignBlockWithReference    ( view, self->block, view->offset ) ; }
        else                 { Array_AssignBlockWithoutReference ( view, self->block, view->offset ) ; }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Scaling.
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifndef _NoNumeric
void _MakeToken ( _ArrayName, _Scale ) ( _ArrayName *self, const _ArrayDataType value )
{
    if ( self != NULL )
    {
        if ( View2D_IsUniform ( self ) ) { Utilities1D_Scale ( self->size, self->data, self->stride1, value ) ; }
        else
        {
            auto Integer r ;
            for ( r = 0 ; r < View2D_Rows ( self ) ; r++ )
            {
                Utilities1D_Scale ( self->extent1, Array2D_RowPointer ( self, r ), self->stride1, value ) ;
            }
        }
    }
}
# endif

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set all items.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _Set ) ( _ArrayName *self, const _ArrayDataType value )
{
    if ( self != NULL )
    {
        if ( View2D_IsUniform ( self ) ) { Utilities1D_Set ( self->size, self->data, self->stride1, value ) ; }
        else
        {
            auto Integer r ;
            for ( r = 0 ; r < View2D_Rows ( self ) ; r++ )
            {
                Utilities1D_Set ( self->extent1, Array2D_RowPointer ( self, r ), self->stride1, value ) ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _SetItem ) ( _ArrayName *self, const Integer i, const Integer j, const _ArrayDataType value, Status *status )
{
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        if ( Integer_IsInRange ( i, 0, self->extent0 ) &&
             Integer_IsInRange ( j, 0, self->extent1 ) ) Array2D_Item ( self, i, j ) = value ;
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set an item from a multislice.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _SetItemMultiSlice ) (       _ArrayName     *self       ,
                                                               const  MultiSlice    *multiSlice ,
                                                               const _ArrayDataType  value      ,
                                                                      Status        *status     )
{
    if ( ( multiSlice != NULL ) && ( multiSlice->capacity >= 2 ) )
    {
        _MakeToken ( _ArrayName, _SetItem ) ( self, multiSlice->items[0].start, multiSlice->items[1].start, value, status ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . View.
! . Extents, starts and strides need to be checked before entry.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _View ) ( const _ArrayName *self          ,
                                        const  Integer    start0        ,
                                        const  Integer    start1        ,
                                        const  Integer    extent0       ,
                                        const  Integer    extent1       ,
                                        const  Integer    stride0       ,
                                        const  Integer    stride1       ,
                                        const  Boolean    withReference ,
                                              _ArrayName *view          ,
                                               Status    *status        )
{
    if ( ( self != NULL ) && ( view != NULL ) && Status_IsOK ( status ) )
    {
        View2D_View ( ( View2D * ) self, start0, start1, extent0, extent1, stride0, stride1, ( View2D * ) view, status ) ;
        if ( withReference ) { Array_AssignBlockWithReference    ( view, self->block, view->offset ) ; }
        else                 { Array_AssignBlockWithoutReference ( view, self->block, view->offset ) ; }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . View 1-D.
! . Offset, extents and strides need to be checked before entry.
! . One of the extents needs to be 1.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _View1D ) ( const _ArrayName   *self          ,
                                          const  Integer      dimension     , /* 0 or 1 */
                                          const  Integer      start0        ,
                                          const  Integer      start1        ,
                                          const  Integer      extent        ,
                                          const  Integer      stride        ,
                                          const  Boolean      withReference ,
                                                _Array1DName *view          ,
                                                 Status      *status        )
{
    if ( ( self != NULL ) && ( view != NULL ) && Status_IsOK ( status ) )
    {
        View2D_View1D ( ( View2D * ) self, dimension, start0, start1, extent, stride, ( View1D * ) view, status ) ; 
        if ( withReference ) { Array_AssignBlockWithReference    ( view, self->block, view->offset ) ; }
        else                 { Array_AssignBlockWithoutReference ( view, self->block, view->offset ) ; }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . View of raw.
! . Offset, extents and strides need to be checked before entry.
! . For limited, local use only as no block.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _ViewOfRaw ) (       _ArrayName     *self    ,
                                             const  Integer        offset  ,
                                             const  Integer        extent0 ,
                                             const  Integer        extent1 ,
                                             const  Integer        stride0 ,
                                             const  Integer        stride1 ,
                                                   _ArrayDataType *data    )
{
    if ( ( self != NULL ) && ( data != NULL ) )
    {
        self->block   = NULL              ;
        self->data    = data              ;
        self->extent0 = extent0           ;
        self->extent1 = extent1           ;
        self->offset  = offset            ;
        self->size    = extent0 * extent1 ;
        self->stride0 = stride0           ;
        self->stride1 = stride1           ;
    }
}

# undef _Array1DName
# undef _Array1DType
# undef _ArrayName
# undef _ArrayType
# undef _BlockType
