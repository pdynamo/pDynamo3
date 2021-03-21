/*==================================================================================================================================
! . Procedures for selections (immutable ordered cardinal sets).
!=================================================================================================================================*/

# include <stdlib.h>

# include "BooleanUtilities.h"
# include "IntegerUtilities.h"
# include "Memory.h"
# include "NumericalMacros.h"
# include "Selection.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void Selection_Initialize ( Selection *self ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
Selection *Selection_Allocate ( const Integer capacity, Status *status )
{
    Selection *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        if ( Integer_IsValidCapacity ( capacity ) )
        {
            self = Memory_AllocateType ( Selection ) ;
            Selection_Initialize ( self ) ;
            if ( ( self != NULL ) && ( capacity > 0 ) )
            {
                self->indices = Integer_Allocate ( capacity, status ) ;
                Integer_Set ( self->indices, capacity, -1 ) ;
                self->capacity = capacity ;
                if ( self->indices == NULL ) Selection_Deallocate ( &self ) ;
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
Integer Selection_Capacity ( const Selection *self ) { return ( ( self == NULL ) ? 0 : self->capacity ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Clear the flags representation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Selection_ClearFlags ( Selection *self ) { if ( self != NULL ) BooleanBlock_Deallocate ( &(self->flags) ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Clear the position representation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Selection_ClearPositions ( Selection *self ) { if ( self != NULL ) IntegerBlock_Deallocate ( &(self->positions) ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Clear representations.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Selection_ClearRepresentations ( Selection *self )
{
    Selection_ClearFlags     ( self ) ;
    Selection_ClearPositions ( self ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
Selection *Selection_Clone ( const Selection *self, Status *status )
{
    Selection *clone = NULL ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        clone = Selection_Allocate ( self->capacity, status ) ;
        if ( clone != NULL ) Integer_CopyTo ( self->indices, self->capacity, clone->indices, NULL ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Complement.
!---------------------------------------------------------------------------------------------------------------------------------*/
Selection *Selection_Complement ( Selection *self, const Integer upperBound, Status *status )
{
    Selection *new = NULL ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto BooleanBlock *flags ;
        auto Integer       m, n ;
        m     = Maximum ( upperBound, Selection_UpperBound ( self ) ) ;
        n     = m - self->capacity ;
        flags = Selection_MakeFlags ( self, m, status ) ;
        new   = Selection_Allocate  (       n, status ) ;
        if ( ( flags != NULL ) && ( new != NULL ) && ( n > 0 ) )
        {
            auto Integer i ;
            for ( i = n = 0 ; i < m ; i++ )
            {
                if ( ! Block_Item ( flags, i ) ) { Selection_Item ( new, n ) = i ; n++ ; }
            }
        }
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Selection_Deallocate ( Selection **self )
{
    if ( (*self) != NULL )
    {
        Selection_ClearRepresentations ( (*self) ) ;
        Integer_Deallocate ( &((*self)->indices) ) ;
        Memory_Deallocate  (   (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Difference.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . "self" is returned if it is unchanged. */
Selection *Selection_Difference ( Selection *self, const Integer number, Selection **others, Status *status )
{
    Selection *new = NULL ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        if ( ( self->capacity == 0 ) || ( number <= 0 ) ) new = self ;
        else
        {
            auto BooleanBlock *flags = Selection_MakeFlags ( self, 0, status ) ;
            if ( flags != NULL )
            {
                auto Boolean   *f ;
                auto Integer    c, i, n, s, u ;
                auto Selection *other ;
                f = Block_Items    ( flags ) ;
                u = Block_Capacity ( flags ) ;
                for ( s = 0 ; s < number ; s++ )
                {
                    other = others[s] ;
                    if ( ( other != NULL ) && ( other->capacity > 0 ) )
                    {
                        for ( c = 0 ; c < other->capacity ; c++ )
                        {
                            i = Selection_Item ( other, c ) ;
                            if ( i < u ) f[i] = False ;
                        }
                    }
                }
                n = Boolean_Count ( f, u, True ) ;
                if ( n < self->capacity ) new = Selection_FromBooleans ( u, f, status ) ;
                else new = self ;
                Selection_ClearFlags ( self ) ;
            }
        }
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Constructor from booleans.
!---------------------------------------------------------------------------------------------------------------------------------*/
Selection *Selection_FromBooleans ( const Integer capacity, const Boolean *flags, Status *status )
{
    Selection *self = NULL ;
    if ( ( capacity > 0 ) && ( flags != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer n = Boolean_Count ( flags, capacity, True ) ;
        self = Selection_Allocate ( n, status ) ;
        if ( self != NULL )
        {
            auto Integer i ;
            for ( i = n = 0 ; i < capacity ; i++ )
            {
                if ( flags[i] ) { Selection_Item ( self, n ) = i ; n++ ; }
            }
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Constructor from cardinals.
!---------------------------------------------------------------------------------------------------------------------------------*/
Selection *Selection_FromIntegers ( const Integer capacity, const Integer *indices, Status *status )
{
    Selection *self = NULL ;
    if ( ( ( capacity == 0 ) || ( ( capacity > 0 ) && ( indices != NULL ) ) ) && Status_IsOK ( status ) )
    {
        self = Selection_Allocate ( capacity, status ) ;
        if ( ( self != NULL ) && ( capacity > 0 ) )
        {
            Integer_CopyTo ( indices, capacity, self->indices, NULL ) ;
            self->capacity = Integer_SortUnique ( self->indices, capacity ) ;
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for emptiness.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean Selection_IsEmpty ( const Selection *self, const Boolean nullIsFull ) { return ( ( self == NULL ) ? ( ! nullIsFull ) : ( self->capacity == 0 ) ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Membership.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean Selection_HasItem ( Selection *self, const Integer value, Status *status )
{
    Boolean flag = False ;
    if ( ( self != NULL ) && ( Integer_IsInRange ( value, 0, Selection_UpperBound ( self ) ) ) && Status_IsOK ( status ) )
    {
        auto BooleanBlock *flags = Selection_MakeFlags ( self, 0, status ) ;
        if ( flags != NULL ) flag = Block_Item ( flags, value ) ;
    }
    return flag ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void Selection_Initialize ( Selection *self )
{
    if ( self != NULL )
    {
        self->capacity  = 0     ;
        self->indices   = NULL  ;
        self->flags     = NULL  ;
        self->positions = NULL  ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Intersection.
!---------------------------------------------------------------------------------------------------------------------------------*/
Selection *Selection_Intersection ( const Integer number, Selection **others, Status *status )
{
    Selection *self = NULL ;
    if ( ( number > 0 ) && Status_IsOK ( status ) )
    {
        Boolean  isOK = True ;
        Integer *frequencies = NULL, i, m = 0, n, s, t ;
        n = Selection_UpperBound ( others[0] ) ;
        for ( i = 1 ; i < number ; i++ ) { n = Minimum ( n, Selection_UpperBound ( others[i] ) ) ; }
        if ( n > 0 )
        {
            frequencies = Integer_Allocate ( n, status ) ;
            isOK        = ( frequencies != NULL ) ;
            if ( isOK )
            {
                auto Selection *other ;
                Integer_Set ( frequencies, n, 0 ) ;
                for ( i = 0 ; i < number ; i++ )
                {
                    other = others[i] ;
                    if ( other != NULL )
                    {
                        for ( s = 0 ; s < other->capacity ; s++ )
                        {
                            t = Selection_Item ( other, s ) ;
                            if ( t < n ) frequencies[t] += 1 ;
                        }
                    }
                }
                m = Integer_Count ( frequencies, n, number ) ;
            }
        }
        if ( isOK )
        {
            self = Selection_Allocate ( m, status ) ;
            if ( self != NULL )
            {
                for ( m = s = 0 ; s < n ; s++ )
                {
                    if ( frequencies[s] == number ) { Selection_Item ( self, m ) = s ; m++ ; }
                }
            }
        }
        Integer_Deallocate ( &frequencies ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the flags representation.
!---------------------------------------------------------------------------------------------------------------------------------*/
BooleanBlock *Selection_MakeFlags ( Selection *self, const Integer upperBound, Status *status )
{
    auto BooleanBlock *flags = NULL ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer n = Maximum ( Selection_UpperBound ( self ), upperBound ) ;
        flags = self->flags ;
        if ( ( flags == NULL ) || ( flags->capacity < n ) )
        {
            Selection_ClearFlags ( self ) ;
            flags = BooleanBlock_Allocate ( n, status ) ;
            if ( flags != NULL )
            {
                auto Integer i ;
                BooleanBlock_Set ( flags, False ) ;
                for ( i = 0 ; i < self->capacity ; i++ ) { Block_Item ( flags, Selection_Item ( self, i ) ) = True ; }
                self->flags = flags ;
            }
        }
    }
    return flags ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the positions representation.
!---------------------------------------------------------------------------------------------------------------------------------*/
IntegerBlock *Selection_MakePositions ( Selection *self, const Integer upperBound, Status *status )
{
    auto IntegerBlock *positions = NULL ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer n = Maximum ( Selection_UpperBound ( self ), upperBound ) ;
        positions = self->positions ;
        if ( ( positions == NULL ) || ( positions->capacity < n ) )
        {
            Selection_ClearPositions ( self ) ;
            positions = IntegerBlock_Allocate ( n, status ) ;
            if ( positions != NULL )
            {
                auto Integer i ;
                IntegerBlock_Set ( positions, -1 ) ;
                for ( i = 0 ; i < self->capacity ; i++ ) { Block_Item ( positions, Selection_Item ( self, i ) ) = i ; }
                self->positions = positions ;
            }
        }
    }
    return positions ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Merging.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Straightforward if needed - just append all indices with renumbering by increment. */

/*----------------------------------------------------------------------------------------------------------------------------------
! . Position.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer Selection_PositionOfItem ( Selection *self, const Integer value, Status *status )
{
    auto Integer position = -1 ;
    if ( ( self != NULL ) && ( Integer_IsInRange ( value, 0, Selection_UpperBound ( self ) ) ) && Status_IsOK ( status ) )
    {
        auto IntegerBlock *positions = Selection_MakePositions ( self, 0, status ) ;
        if ( positions != NULL ) position = Block_Item ( positions, value ) ;
    }
    return position ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Pruning.
!---------------------------------------------------------------------------------------------------------------------------------*/
Selection *Selection_Prune ( Selection *self, Selection *toKeep, Status *status )
{
    Selection *new = NULL ;
    if ( ( self != NULL ) && ( toKeep != NULL ) && Status_IsOK ( status ) )
    {
        IntegerBlock *positions ;
        Selection    *others[2] ;
        others[0] = self ; others[1] = toKeep ;
        new       = Selection_Intersection  ( 2, others, status ) ;
        positions = Selection_MakePositions ( toKeep, 0, status ) ; 
        if ( ( new != NULL ) && ( new->capacity > 0 ) && ( positions != NULL ) )
        {
            auto Integer i, s ;
            for ( i = 0 ; i < new->capacity ; i++ )
            {
                s = Selection_Item ( new, i ) ;
                Selection_Item ( new, i ) = Block_Item ( positions, s ) ;
            }
        }
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Symmetric difference.
!---------------------------------------------------------------------------------------------------------------------------------*/
Selection *Selection_SymmetricDifference ( const Integer number, Selection **others, Status *status )
{
    Selection *self = NULL ;
    if ( ( number > 0 ) && Status_IsOK ( status ) )
    {
        Boolean  isOK = True ;
        Integer *frequencies = NULL, i, m = 0, n = 0, s ;
        for ( i = 0 ; i < number ; i++ ) { n = Maximum ( n, Selection_UpperBound ( others[i] ) ) ; }
        if ( n > 0 )
        {
            frequencies = Integer_Allocate ( n, status ) ;
            isOK        = ( frequencies != NULL ) ;
            if ( isOK )
            {
                auto Selection *other ;
                Integer_Set ( frequencies, n, 0 ) ;
                for ( i = 0 ; i < number ; i++ )
                {
                    other = others[i] ;
                    if ( other != NULL )
                    {
                        for ( s = 0 ; s < other->capacity ; s++ ) frequencies[Selection_Item ( other, s )] += 1 ;
                    }
                }
                m = Integer_Count ( frequencies, n, 1 ) ;
            }
        }
        if ( isOK )
        {
            self = Selection_Allocate ( m, status ) ;
            if ( self != NULL )
            {
                for ( m = s = 0 ; s < n ; s++ )
                {
                    if ( frequencies[s] == 1 ) { Selection_Item ( self, m ) = s ; m++ ; }
                }
            }
        }
        Integer_Deallocate ( &frequencies ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Union.
!---------------------------------------------------------------------------------------------------------------------------------*/
Selection *Selection_Union ( const Integer number, Selection **others, Status *status )
{
    Selection *self = NULL ;
    if ( ( number > 0 ) && Status_IsOK ( status ) )
    {
        Boolean  *flags = NULL, isOK = True ;
        Integer  i, n = 0, s ;
        for ( i = 0 ; i < number ; i++ ) { n = Maximum ( n, Selection_UpperBound ( others[i] ) ) ; }
        if ( n > 0 )
        {
            flags = Boolean_Allocate ( n, status ) ;
            isOK  = ( flags != NULL ) ;
            if ( isOK )
            {
                auto Selection *other ;
                Boolean_Set ( flags, n, False ) ;
                for ( i = 0 ; i < number ; i++ )
                {
                    other = others[i] ;
                    if ( other != NULL )
                    {
                        for ( s = 0 ; s < other->capacity ; s++ ) { flags[Selection_Item ( other, s )] = True ; }
                    }
                }
            }
        }
        if ( isOK ) self = Selection_FromBooleans ( n, flags, status ) ;
        Boolean_Deallocate ( &flags ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Return the upper bound for the container.
! . This is the value of the largest index plus one and will be the minimal size of flag and position representations.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer Selection_UpperBound ( const Selection *self )
{
    return ( ( ( self != NULL ) && ( self->capacity > 0 ) ) ? self->indices[self->capacity-1] + 1 : 0 ) ;
}
