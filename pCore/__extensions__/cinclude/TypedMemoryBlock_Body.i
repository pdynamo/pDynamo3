/*==================================================================================================================================
! . Basic typed memory block functions.
!=================================================================================================================================*/

/* . Need _CoreDataType, _CoreDataTypeInitializer. */

# ifdef _TMBDEBUG
# include <stdio.h>
# endif

# include "Integer.h"
# include "Memory.h"
# include "Status.h"
# include "TemplateMacros.h"

# define _BlockType Block
# define _BlockName _MakeToken ( _CoreDataType, _BlockType )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
_BlockName * _MakeToken ( _BlockName, _Allocate ) ( const Integer capacity, Status *status )
{
    _BlockName *self = NULL ;
    if ( Integer_IsValidCapacity ( capacity ) && Status_IsOK ( status ) )
    {
        _CoreDataType *items = NULL ;
        self  = Memory_AllocateType         ( _BlockName              ) ;
        items = Memory_AllocateArrayOfTypes ( capacity, _CoreDataType ) ;
        if ( ( self == NULL ) || ( ( capacity > 0 ) && ( items == NULL ) ) )
        {
            Memory_Deallocate ( items ) ;
            Memory_Deallocate ( self  ) ;
            Status_Set ( status, Status_OutOfMemory ) ;
        }
        else
        {
            self->capacity   = capacity ;
            self->references = 0        ;
            self->items      = items    ;
/*            _MakeToken ( _BlockName, _Reference ) ( self ) ;*/
# ifdef _TMBDEBUG
            printf ( "** Allocating %d (%zu) **\n", self->capacity, sizeof ( _CoreDataType ) ) ;
            fflush ( stdout ) ;
# endif
        }
    }
    else Status_Set ( status, Status_InvalidArgument ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
_BlockName * _MakeToken ( _BlockName, _Clone ) ( const _BlockName *self, Status *status )
{
    _BlockName *clone = NULL ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        clone = _MakeToken ( _BlockName, _Allocate ) ( self->capacity, status ) ;
        _MakeToken ( _BlockName, _CopyTo ) ( self, clone, status ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copying.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _BlockName, _CopyTo ) ( const _BlockName *self, _BlockName *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) && Status_IsOK ( status ) )
    {
        if ( self->capacity == other->capacity )
        {
            auto Integer i ;
            for ( i = 0 ; i < self->capacity ; i++ ) other->items[i] = self->items[i] ;
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _BlockName, _Deallocate ) ( _BlockName **self )
{
    if ( (*self) != NULL )
    {
        _MakeToken ( _BlockName, _Dereference ) ( (*self) ) ;
        if ( (*self)->references <= 0 )
        {
# ifdef _TMBDEBUG
            printf ( "** Deallocating %d (%zu) **\n", (*self)->capacity, sizeof ( _CoreDataType ) ) ;
            fflush ( stdout ) ;
# endif
            Memory_Deallocate ( (*self)->items ) ;
            Memory_Deallocate ( (*self)        ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Dereferencing.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _BlockName, _Dereference ) ( _BlockName *self )
{
    if ( self != NULL )
    {
        self->references -= 1 ;
# ifdef _TMBDEBUGR
        printf ( "** Dereferencing %d for %d\n", self->references, self->capacity ) ;
        fflush ( stdout ) ;
# endif
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
_CoreDataType _MakeToken ( _BlockName, _GetItem ) ( const _BlockName *self, const Integer index, Status *status )
{
    _CoreDataType value = _CoreDataTypeInitializer ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        if ( Integer_IsInRange ( index, 0, self->capacity ) ) value = self->items[index] ;
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
    return value ;
}

# ifndef _ExcludeNegate
/*----------------------------------------------------------------------------------------------------------------------------------
! . Negation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _BlockName, _Negate ) ( _BlockName *self )
{
    if ( self != NULL )
    {
        auto Integer       i     ;
        auto _CoreDataType value ;
        for ( i = 0 ; i < self->capacity ; i++ )
        {
            value = self->items[i] ;
            self->items[i] = _UnaryNegationOperator value ;
        }
    }
}
# endif

/*----------------------------------------------------------------------------------------------------------------------------------
! . Referencing.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _BlockName, _Reference ) ( _BlockName *self )
{
    if ( self != NULL )
    {
        self->references += 1 ;
# ifdef _TMBDEBUGR
        printf ( "** Referencing %d for %d\n", self->references, self->capacity ) ;
        fflush ( stdout ) ;
# endif
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Resizing.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _BlockName, _Resize ) (       _BlockName  *self     ,
                                          const  Integer     capacity ,
                                                 Status     *status   )
{
    if ( ( self != NULL ) && ( capacity != self->capacity ) )
    {
        if ( Integer_IsValidCapacity ( capacity ) )
        {
            self->items = Memory_ReallocateArrayOfTypes ( self->items, capacity, _CoreDataType ) ;
            if ( ( self->items == NULL ) && ( capacity > 0 ) ) { self->capacity = 0 ; Status_Set ( status, Status_OutOfMemory ) ; }
            else self->capacity = capacity ;
        }
        else Status_Set ( status, Status_InvalidArgument ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set all items.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _BlockName, _Set ) ( _BlockName *self, const _CoreDataType value )
{
    if ( self != NULL )
    {
        auto Integer  i ;
        for ( i = 0 ; i < self->capacity ; i++ ) self->items[i] = value ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _BlockName, _SetItem ) ( _BlockName *self, const Integer  index, const _CoreDataType value, Status *status )
{
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        if ( Integer_IsInRange ( index, 0, self->capacity ) ) self->items[index] = value ;
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
}

# undef _BlockName
# undef _BlockType
