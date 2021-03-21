/*==================================================================================================================================
! . A container for selections.
!=================================================================================================================================*/

# include "BooleanUtilities.h"
# include "Memory.h"
# include "NumericalMacros.h"
# include "SelectionContainer.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void SelectionContainer_Initialize ( SelectionContainer *self ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
SelectionContainer *SelectionContainer_Allocate ( const Integer capacity, Status *status )
{
    SelectionContainer *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        if ( Integer_IsValidCapacity ( capacity ) )
        {
            self = Memory_AllocateType ( SelectionContainer ) ;
            SelectionContainer_Initialize ( self ) ;
            if ( ( self != NULL ) && ( capacity > 0 ) )
            {
	        self->items = Memory_AllocateArrayOfReferences ( capacity, Selection ) ;
                if ( self->items == NULL ) SelectionContainer_Deallocate ( &self ) ;
                else
                {
                    auto Integer i ;
                    for ( i = 0 ; i < capacity ; i++ ) self->items[i] = NULL ;
                    self->capacity = capacity ;
                }
            }
            if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
        }
        else Status_Set ( status, Status_InvalidArgument ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Capacity.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer SelectionContainer_Capacity ( const SelectionContainer *self ) { return ( ( self == NULL ) ? 0 : self->capacity ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
SelectionContainer *SelectionContainer_Clone ( const SelectionContainer *self, Status *status )
{
    SelectionContainer *clone = NULL ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        clone = SelectionContainer_Allocate ( self->capacity, status ) ;
        if ( clone != NULL )
        {
            auto Boolean isOK = True ;
            auto Integer i ;
            auto Selection *item ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                item = Selection_Clone ( self->items[i], status ) ;
                if ( item == NULL ) { isOK = False ; break ; }
                else clone->items[i] = item ;
            }
            if ( ! isOK ) SelectionContainer_Deallocate ( &clone ) ;
        }
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SelectionContainer_Deallocate ( SelectionContainer **self )
{
    if ( (*self) != NULL )
    {
        auto Integer i ;
        for ( i = 0 ; i < (*self)->capacity ; i++ ) Selection_Deallocate ( &((*self)->items[i]) ) ;
        Memory_Deallocate ( (*self)->items ) ;
        Memory_Deallocate ( (*self)        ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Constructor of single-item selections given a capacity.
!---------------------------------------------------------------------------------------------------------------------------------*/
SelectionContainer *SelectionContainer_FromCapacity ( const Integer capacity, Status *status )
{
    SelectionContainer *self = NULL ;
    if ( ( capacity > 0 ) && Status_IsOK ( status ) )
    {
        self = SelectionContainer_Allocate ( capacity, status ) ;
        if ( self != NULL )
        {
            auto Boolean    isOK = True ;
            auto Integer    i ;
            auto Selection *item ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                item = Selection_FromIntegers ( 1, &i, status ) ;
                if ( item == NULL ) { isOK = False ; break ; }
                else self->items[i] = item ;
            }
            if ( ! isOK ) SelectionContainer_Deallocate ( &self ) ;
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . In-place fusion of items with indices in toFuse.
! . The boolean array toFuse can be created with a call to MakeMembershipFlags.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean SelectionContainer_FuseItems ( SelectionContainer *self, BooleanBlock *toFuse, Status *status )
{
    auto Boolean isOK = True ;
    if ( ( self != NULL ) && ( toFuse != NULL ) && Status_IsOK ( status ) )
    {
        auto Selection **others = Memory_AllocateArrayOfReferences ( self->capacity, Selection ) ;
        if ( others != NULL )
        {
            auto Integer i, n ;
            for ( i = n = 0 ; i < self->capacity ; i++ )
            {
                if ( Block_Item ( toFuse, i ) ) { others[n] = self->items[i] ; n += 1 ; }
            }
            if ( n > 1 )
            {
                auto Selection *item, *new = Selection_Union ( n, others, status ) ;
                if ( new != NULL )
                {
                    for ( i = n = 0 ; i < self->capacity ; i++ )
                    {
                        item = self->items[i] ;
                        if ( Block_Item ( toFuse, i ) ) Selection_Deallocate ( &item ) ;
                        else { self->items[n] = item ; n++ ; }
                    }
                    self->items[n] = new ;
                    self->capacity = n+1 ;
                }
                else isOK = False ;
            }
            Memory_Deallocate ( others ) ;
        }
        else isOK = False ;
    }
    return isOK ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void SelectionContainer_Initialize ( SelectionContainer *self )
{
    if ( self != NULL )
    {
        self->capacity = 0    ;
        self->items    = NULL ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make an array of flags indicating whether an item contains indices in the members selection.
! . Two tests are permitted:
!   - the AND test flags an item in which all members belong to the selection                (andTest = True )
!   - the OR test flags an item in which at least one of the members belong to the selection (andTest = False)
!---------------------------------------------------------------------------------------------------------------------------------*/
BooleanBlock *SelectionContainer_MakeMembershipFlags (       SelectionContainer *self    ,
                                                             Selection          *members ,
                                                       const Boolean             andTest ,
                                                             Status             *status  )
{
    auto BooleanBlock *flags = NULL ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        flags = BooleanBlock_Allocate ( self->capacity, status ) ;
        if ( flags != NULL )
        {
            BooleanBlock_Set ( flags, False ) ;
            if ( members != NULL )
            {
                auto BooleanBlock *memberFlags = Selection_MakeFlags ( members, SelectionContainer_UpperBound ( self ), status ) ;
                if ( memberFlags != NULL )
                {
                    auto Integer    i, n, s ;
                    auto Selection *item ;
                    for ( i = 0 ; i < self->capacity ; i++ )
                    {
                        item = self->items[i] ;
                        for ( n = s = 0 ;  s < item->capacity ; s++ )
                        {
                            if ( Block_Item ( memberFlags, Selection_Item ( item, s ) ) ) n += 1 ;
                        }
                        if ( andTest ) { Block_Item ( flags, i ) = ( n >= item->capacity ) ; }
                        else           { Block_Item ( flags, i ) = ( n >  0              ) ; }
                    }
                }
            }
        }
    }
    return flags ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . In-place removal of items with indices in toRemove.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean SelectionContainer_RemoveItems ( SelectionContainer *self, Selection *toRemove, Status *status )
{
    auto Boolean isOK = True ;
    if ( ( self != NULL ) && ( toRemove != NULL ) && Status_IsOK ( status ) )
    {
        auto BooleanBlock *flags ;
        flags = Selection_MakeFlags ( toRemove, SelectionContainer_UpperBound ( self ), status ) ;
        if ( flags != NULL )
        {
            auto Boolean    isFlagged ;
            auto Integer    i, n = 0, s ;
            auto Selection *item ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                isFlagged = False ;
                item      = self->items[i] ;
                for ( s = 0 ; s < item->capacity ; s++ )
                {
                    if ( Block_Item ( flags, Selection_Item ( item, s ) ) ) { isFlagged = True ; break ; }
                }
                if ( isFlagged ) Selection_Deallocate ( &item ) ;
                else { self->items[n] = item ; n++ ; }
            }
            self->capacity =  n ;
        }
        else isOK = False ;
    }
    return isOK ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the union of all items with indices in toUnion.
!---------------------------------------------------------------------------------------------------------------------------------*/
Selection *SelectionContainer_UnionOfItems ( const SelectionContainer *self, Selection *toUnion, Status *status )
{
    Selection *new = NULL ;
    if ( ( self != NULL ) && ( toUnion != NULL ) && Status_IsOK ( status ) )
    {
        auto BooleanBlock  *flags  = Selection_MakeFlags ( toUnion, SelectionContainer_UpperBound ( self ), status ) ;
        auto Selection    **others = Memory_AllocateArrayOfReferences ( self->capacity, Selection ) ;
        if ( ( flags != NULL ) && ( others != NULL ) )
        {
            auto Integer    i, n, s ;
            auto Selection *item ;
            for ( i = n = 0 ; i < self->capacity ; i++ )
            {
                item = self->items[i] ;
                for ( s = 0 ; s < item->capacity ; s++ )
                {
                    if ( Block_Item ( flags, Selection_Item ( item, s ) ) ) { others[n] = item ; n += 1 ; break ; }
                }
            }
            new = Selection_Union ( n, others, status ) ;
            Memory_Deallocate ( others ) ;
        }
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Determine the upper bound of the container.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer SelectionContainer_UpperBound ( const SelectionContainer *self )
{
    Integer upperBound = 0 ;
    if ( self != NULL )
    {
        auto Integer i ;
        for ( i = 0 ; i < self->capacity ; i++ ) { upperBound = Maximum ( Selection_UpperBound ( self->items[i] ), upperBound ) ; }
    }
    return upperBound ;
}
