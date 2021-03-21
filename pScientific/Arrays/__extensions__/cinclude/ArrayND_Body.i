/*==================================================================================================================================
! . N-D arrays.
!=================================================================================================================================*/

/* . Need _ArrayDataType, _ArrayDataTypeInitializer. */

# include <stdlib.h>
# include <stdio.h>

# include "Array_Macros.h"
# include "Integer.h"
# include "Memory.h"
# include "NumericalMacros.h"
# include "Status.h"
# include "TemplateMacros.h"

# define _Array1DType   Array1D
# define _Array2DType   Array2D
# define _ArrayType     ArrayND
# define _Array1DName  _MakeToken ( _ArrayDataType, _Array1DType )
# define _Array2DName  _MakeToken ( _ArrayDataType, _Array2DType )
# define _ArrayName    _MakeToken ( _ArrayDataType, _ArrayType   )
# define _BlockType    _MakeToken ( _ArrayDataType,  Block       )
# define _IteratorName _MakeToken ( _ArrayDataType, Iterator     )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation (type only).
!---------------------------------------------------------------------------------------------------------------------------------*/
_ArrayName * _MakeToken ( _ArrayName, _Allocate ) ( Status *status )
{
    _ArrayName *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( _ArrayName ) ;
        if ( self != NULL )
        {
            self->view  = NULL ;
            self->block = NULL ;
            self->data  = NULL ;
        }
        else Status_Set ( status, Status_OutOfMemory ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation with a rank (view only).
!---------------------------------------------------------------------------------------------------------------------------------*/
_ArrayName * _MakeToken ( _ArrayName, _AllocateWithRank ) ( const Integer rank, Status *status )
{
    _ArrayName *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        auto ViewND *view = NULL ;
        view = ViewND_AllocateWithRank ( rank, status ) ;
        self = _MakeToken ( _ArrayName, _FromViewBlock ) ( view, NULL, False, status ) ;
        if ( self == NULL ) ViewND_Deallocate ( &view  ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation with a shape (view and block).
!---------------------------------------------------------------------------------------------------------------------------------*/
_ArrayName * _MakeToken ( _ArrayName, _AllocateWithShape ) ( const Integer rank, const Integer *extents, Status *status )
{
    _ArrayName *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        auto _BlockType *block = NULL ;
        auto ViewND     *view  = NULL ;
        view  = ViewND_AllocateWithShape                  ( rank, extents       , status ) ;
        block = _MakeToken ( _BlockType, _Allocate      ) ( ViewND_Size ( view ), status ) ;
        self  = _MakeToken ( _ArrayName, _FromViewBlock ) ( view, block, True   , status ) ;
        if ( self == NULL )
        {
            ViewND_Deallocate                      ( &view  ) ;
            _MakeToken ( _BlockType, _Deallocate ) ( &block ) ;
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deep cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
_ArrayName * _MakeToken ( _ArrayName, _CloneDeep ) ( const _ArrayName *self, Status *status )
{
    _ArrayName *clone = NULL ;
    if ( ( self != NULL ) && ( self->view != NULL ) && Status_IsOK ( status ) )
    {
        clone = _MakeToken ( _ArrayName, _AllocateWithShape ) ( self->view->rank, self->view->extents, status ) ;
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
        auto ViewND *view = ViewND_Clone ( self->view, status ) ;
        clone = _MakeToken ( _ArrayName, _FromViewBlock ) ( view, self->block, True, status ) ;
        if ( clone == NULL ) ViewND_Deallocate ( &view ) ;
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
        if ( ViewND_AreConformable ( self->view, other->view, status ) )
        {
            auto Iterator *oIterator, *sIterator ;
            oIterator = ViewND_MakeIterator ( other->view, status ) ;
            sIterator = ViewND_MakeIterator ( self->view , status ) ;
            _MakeToken ( _IteratorName, _CopyTo ) ( sIterator, self->data, oIterator, other->data, status ) ;
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
        ViewND_Deallocate ( &((*self)->view) ) ;
        _MakeToken ( _BlockType, _Deallocate ) ( &((*self)->block) ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Constructor given a view and, optionally, a block.
! . No checking for the moment - add capacity check between view and block?
!---------------------------------------------------------------------------------------------------------------------------------*/
_ArrayName * _MakeToken ( _ArrayName, _FromViewBlock ) (       ViewND     *view          ,
                                                               _BlockType *block         ,
                                                         const Boolean     withReference ,
                                                               Status     *status        )
{
    _ArrayName *self = NULL ;
    if ( ( view != NULL ) && Status_IsOK ( status ) )
    {
        self = _MakeToken ( _ArrayName, _Allocate ) ( status ) ;
        if ( self != NULL )
        {
            self->view = view ;
            if ( block != NULL )
            {
                auto Integer offset = self->view->offset ;
                if ( withReference ) { Array_AssignBlockWithReference    ( self, block, offset ) ; }
                else                 { Array_AssignBlockWithoutReference ( self, block, offset ) ; }
            }
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
_ArrayDataType _MakeToken ( _ArrayName, _GetItem ) ( const _ArrayName *self, const Integer *indices, Status *status )
{
    _ArrayDataType value = _ArrayDataTypeInitializer ;
    if ( ( self != NULL ) && ( indices != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer n ;
        n = ViewND_GetIndex ( self->view, self->view->rank, indices, status ) ;
        if ( n >= 0 ) value = self->data[n] ;
    }
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get an item using a scalar multislice.
!---------------------------------------------------------------------------------------------------------------------------------*/
_ArrayDataType _MakeToken ( _ArrayName, _GetItemMultiSlice ) ( const _ArrayName *self, const MultiSlice *multiSlice, Status *status )
{
    _ArrayDataType value = _ArrayDataTypeInitializer ;
    if ( ( self != NULL ) && ( multiSlice != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer n ;
        n = ViewND_GetIndexMultiSlice ( self->view, multiSlice, status ) ;
        if ( n >= 0 ) value = self->data[n] ;
    }
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Resizing.
! . This is only done if the array owns the block and their are no other views.
! . Care should be taken not to resize when there are unreferenced views.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _Resize ) (       _ArrayName  *self    ,
                                          const  Integer     extent0 ,
                                                 Status     *status  )
{
    if ( ( self != NULL ) && ( self->view != NULL ) && Status_IsOK ( status ) )
    {
        auto _BlockType *block = self->block ;
        auto ViewND     *view  = self->view  ;
        if ( ( view->rank        >  0               ) &&
             ( view->extents[0]  >  0               ) &&
             ( block             != NULL            ) &&
             ( block->references <= 2               ) &&
             ( view->size        == block->capacity ) )
        {
            if ( extent0 != view->extents[0] )
            {
                auto Integer newSize = extent0 * ( view->size / view->extents[0] ) ;
                _MakeToken ( _BlockType, _Resize ) ( block, newSize, status ) ;
                if ( Status_IsOK ( status ) )
                {
                    view->extents[0] = extent0 ;
                    view->size       = newSize ;
                    self->data       = Array_BlockDataPointer ( self, view->offset ) ;
                }
                else
                {
                    _MakeToken ( _BlockType, _Deallocate ) ( &block ) ;
                    self->block = NULL ;
                    self->data  = NULL ;
                }
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
    if ( ( self != NULL ) && ( self->view != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer oldE = self->view->extents[0] ;
        _MakeToken ( _ArrayName, _Resize ) ( self, extent0, status ) ;
        if ( ( oldE < extent0 ) && Status_IsOK ( status ) )
        {
            /* . A fudge until ArrayNDs have same structure as 1D and 2D - works because resizing only valid on compact arrays. */
            auto Integer i, start, stop ;
            start = ViewND_Index1D( self->view, oldE ) ;
            stop  = self->view->size ;
            for ( i = start ; i < stop ; i++ ) self->data[i] = initializer ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set all items.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _Set ) ( _ArrayName *self, const _ArrayDataType value, Status *status )
{
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto Iterator *iterator = ViewND_MakeIterator ( self->view, status ) ;
        _MakeToken ( _IteratorName, _Set ) ( iterator, self->data, value, status ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _SetItem ) ( _ArrayName *self, const Integer *indices, const _ArrayDataType value, Status *status )
{
    if ( ( self != NULL ) && ( self->view != NULL ) && ( indices != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer n ;
        n = ViewND_GetIndex ( self->view, self->view->rank, indices, status ) ;
        if ( n >= 0 ) self->data[n] = value ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set an item using a scalar multislice.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _ArrayName, _SetItemMultiSlice ) ( _ArrayName *self, const MultiSlice *multiSlice, const _ArrayDataType value, Status *status )
{
    if ( ( self != NULL ) && ( multiSlice != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer n ;
        n = ViewND_GetIndexMultiSlice ( self->view, multiSlice, status ) ;
        if ( n >= 0 ) self->data[n] = value ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Tail view 2-D.
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifndef _NoViewTail2D
void _MakeToken ( _ArrayName, _ViewTail2D ) ( const _ArrayName   *self          ,
                                                    Integer      *indices       ,
                                              const Boolean       withReference ,
                                                    _Array2DName *tail          ,
                                                    Status       *status        )
{
    if ( ( self != NULL ) && ( indices != NULL ) && ( tail != NULL ) && Status_IsOK ( status ) )
    {
        ViewND_ViewTail2D ( self->view, indices, ( View2D * ) tail, status ) ;
        if ( Status_IsOK ( status ) )
        {
            if ( withReference ) { Array_AssignBlockWithReference    ( tail, self->block, tail->offset ) ; }
            else                 { Array_AssignBlockWithoutReference ( tail, self->block, tail->offset ) ; }
        }
    }
}
# endif

# undef _Array1DName
# undef _Array1DType
# undef _Array2DName
# undef _Array2DType
# undef _ArrayName
# undef _ArrayType
# undef _BlockType
# undef _IteratorName
