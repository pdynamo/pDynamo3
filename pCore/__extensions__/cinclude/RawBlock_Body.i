/*==================================================================================================================================
! . Basic raw array functions.
!=================================================================================================================================*/

/* . Need _CoreDataType, _CoreDataTypeInitializer. */

# include "Integer.h"
# include "Memory.h"
# include "Size.h"
# include "Status.h"
# include "TemplateMacros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
_CoreDataType * _MakeToken ( _CoreDataType, _Allocate ) ( const Integer capacity, Status *status )
{
    _CoreDataType *self = NULL ;
    if ( Integer_IsValidCapacity ( capacity ) && Status_IsOK ( status ) )
    {
        self = Memory_AllocateArrayOfTypes ( capacity, _CoreDataType ) ;
        if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    }
    else Status_Set ( status, Status_InvalidArgument ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
_CoreDataType * _MakeToken ( _CoreDataType, _Clone ) ( const _CoreDataType *self, const Integer capacity, Status *status )
{
    _CoreDataType *clone = NULL ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        clone = _MakeToken ( _CoreDataType, _Allocate ) ( capacity, status ) ;
        _MakeToken ( _CoreDataType, _CopyTo ) ( self, capacity, clone, status ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copying.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _CoreDataType, _CopyTo ) ( const _CoreDataType *self, const Integer capacity, _CoreDataType *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer i ;
        for ( i = 0 ; i < capacity ; i++ ) other[i] = self[i] ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Counting.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer _MakeToken ( _CoreDataType, _Count ) ( const _CoreDataType  *self, const Integer capacity, const _CoreDataType value )
{
    auto Integer n = 0 ;
    if ( ( self != NULL ) && ( capacity > 0 ) )
    {
        auto Integer i ;
        for ( i = 0 ; i < capacity ; i++ ) { if ( self[i] == value ) n++ ; }
    }
    return n ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _CoreDataType, _Deallocate ) ( _CoreDataType **self )
{
    if ( (*self) != NULL ) Memory_Deallocate ( (*self) ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
_CoreDataType _MakeToken ( _CoreDataType, _GetItem ) ( const _CoreDataType *self, const Integer capacity, const Integer  index, Status *status )
{
    _CoreDataType value = _CoreDataTypeInitializer ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        if ( Integer_IsInRange ( index, 0, capacity ) ) value = self[index] ;
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set all items.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _CoreDataType, _Set ) ( _CoreDataType *self, const Integer capacity, const _CoreDataType value )
{
    if ( self != NULL )
    {
        auto Integer i ;
        for ( i = 0 ; i < capacity ; i++ ) self[i] = value ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
void _MakeToken ( _CoreDataType, _SetItem ) ( _CoreDataType  *self, const Integer capacity, const Integer  index, const _CoreDataType value, Status *status )
{
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        if ( Integer_IsInRange ( index, 0, capacity ) ) self[index] = value ;
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
}

# ifndef _ExcludeSort
/*----------------------------------------------------------------------------------------------------------------------------------
! . Sort items in ascending order.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Integer _MakeToken ( _CoreDataType, _Compare ) ( const void *vSelf, const void *vOther )
{
    Integer result ;
    _CoreDataType *self, *other ;
    self  = ( _CoreDataType * ) vSelf  ;
    other = ( _CoreDataType * ) vOther ;
         if ( (*self) < (*other) ) result = -1 ;
    else if ( (*self) > (*other) ) result =  1 ;
    else result = 0 ;
    return result ;
}

void _MakeToken ( _CoreDataType, _Sort ) ( _CoreDataType *self, const Integer capacity )
{
    if ( ( self != NULL ) && ( capacity > 1 ) )
    {
        qsort ( ( void * ) self, ( Size ) capacity, SizeOf ( _CoreDataType ), ( void * ) _MakeToken ( _CoreDataType, _Compare ) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Sort and remove duplicates.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer  _MakeToken ( _CoreDataType, _SortUnique ) ( _CoreDataType  *self, const Integer capacity )
{
    auto Integer n = 0 ;
    if ( ( self != NULL ) && ( capacity > 0 ) )
    {
        if ( capacity > 1 )
        {
            auto Integer       i ;
            auto _CoreDataType u, v ;
            _MakeToken ( _CoreDataType, _Sort ) ( self, capacity ) ;
            u = self[0] ;
            for ( i = n = 1 ; i < capacity ; i++ )
            {
                v = self[i] ;
                if ( v != u ) { self[n] = v ; u = v ; n++ ; }
            }
        }
        else n = 1 ;
    }
    return n ;
}

# endif
