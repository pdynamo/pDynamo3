/*==================================================================================================================================
! . Iterator.
!=================================================================================================================================*/

# include "Iterator.h"
# include "Memory.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
Iterator *Iterator_Allocate ( Status *status )
{
    Iterator *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( Iterator ) ;
        if ( self != NULL )
        {
/*            self->readOnly      = True ; */
            self->isRegular     = False ;
            self->extent        = 0     ;
            self->numberOfLoops = 0     ;
            self->size          = 0     ;
            self->type          = NULL  ;
            self->vSelf         = NULL  ;
        }
        else Status_Set ( status, Status_OutOfMemory ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
Iterator *Iterator_Clone ( const Iterator *self, Status *status )
{
    Iterator *clone = NULL ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        const IteratorType *type = self->type ;
        clone = Iterator_Allocate ( status ) ;
        if ( clone != NULL )
        {
/*            clone->readOnly      = self->readOnly      ; */
            clone->isRegular     = self->isRegular     ;
            clone->extent        = self->extent        ;
            clone->numberOfLoops = self->numberOfLoops ;
            clone->size          = self->size          ;
            clone->type          = type                ;
            if ( ( type->Clone != NULL ) && ( self->vSelf != NULL ) )
            {
                clone->vSelf = type->Clone ( self->vSelf, status ) ;
                if ( clone->vSelf == NULL )
                {
                   Iterator_Deallocate ( &clone ) ;
                   Status_Set ( status, Status_OutOfMemory ) ;
                }
            }
        }
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Return the current index - no changes to the iterator are made.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer Iterator_CurrentIndex ( const Iterator *self )
{
    Integer c = -1 ;
    if ( self != NULL )
    {
        auto const IteratorType   *type                = self->type         ;
        auto       Integer      ( *current ) ( void* ) = type->CurrentIndex ;
        c = current ( self->vSelf ) ;
    }
    return c ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the data offset.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer Iterator_DataOffSet ( const Iterator *self )
{
    Integer offset = 0 ;
    if ( ( self                   != NULL ) &&
         ( self->type             != NULL ) &&
         ( self->type->DataOffSet != NULL ) ) offset = self->type->DataOffSet ( self->vSelf ) ;
    return offset ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Iterator_Deallocate ( Iterator **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        const IteratorType *type = (*self)->type  ;
              void        *vSelf = (*self)->vSelf ;
        if ( type->Deallocate != NULL ) type->Deallocate ( &vSelf ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Dump the iterator.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _NStart 5 /* 6 */
Integer *Iterator_Dump ( const Iterator *self   ,
                               Integer  *n      ,
                               Status   *status )
{
    Integer *state = NULL ;
    if ( n != NULL ) (*n) = 0 ;
    if ( ( self             != NULL ) &&
         ( self->type       != NULL ) &&
         ( self->type->Dump != NULL ) &&
         Status_IsOK ( status ) )
    {
        state = self->type->Dump ( self->vSelf, _NStart, n, status ) ;
        if ( state != NULL )
        {
            state[0] = Iterator_TypeToInteger ( self->type, status ) ;
/*            state[1] = ( self->readOnly ? 1 : 0 ) ; */
            state[1] = ( self->isRegular ? 1 : 0 ) ;
            state[2] = self->extent        ;
            state[3] = self->numberOfLoops ;
            state[4] = self->size          ;
        }
    }
    return state ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the size.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer Iterator_GetSize ( const Iterator *self ) { return ( (self) == NULL ? 0 : (self)->size ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Load the iterator.
!---------------------------------------------------------------------------------------------------------------------------------*/
Iterator *Iterator_Load ( const Integer   n      ,
                          const Integer  *state  ,
                                Status   *status )
{
    Iterator *self = NULL ;
    if ( ( n > 0 ) && ( state != NULL ) && Status_IsOK ( status ) )
    {
        const IteratorType *type = Iterator_TypeFromInteger ( state[0], status ) ;
        self = Iterator_Allocate ( status ) ;
        if ( ( self != NULL ) && ( type != NULL ) )
        {
/*            self->readOnly               = ( state[1] == 0 ? False : True ) ; */
            self->isRegular     = ( state[1] == 0 ? False : True ) ;
            self->extent        = state[2] ;
            self->numberOfLoops = state[3] ;
            self->size          = state[4] ;
            self->type          = type ;
            self->vSelf = self->type->Load ( _NStart, n, state, status ) ;
        }
        if ( ! Status_IsOK ( status ) ) Iterator_Deallocate ( &self ) ;
    }
    return self ;
}
# undef _NStart

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the next index.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer Iterator_NextIndex ( Iterator *self )
{
    Integer n = -1 ;
    if ( self != NULL )
    {
        auto const IteratorType   *type             = self->type      ;
        auto       Integer      ( *next ) ( void* ) = type->NextIndex ;
        n = next ( self->vSelf ) ;
    }
    return n ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Reset the iterator.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Iterator_Reset ( Iterator *self )
{
    if ( self != NULL )
    {
        auto const IteratorType *type = self->type ;
        type->Reset ( self->vSelf ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Type <-> integer functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
const IteratorType *Iterator_TypeFromInteger ( const Integer  n      ,
                                                     Status  *status )
{
    const IteratorType *type = NULL ;
    if ( ( n >= 0 ) && Status_IsOK ( status ) )
    {
             if ( n == 0 ) type = &IteratorType_Regular1D ;
        else if ( n == 1 ) type = &IteratorType_RegularND ;
        else if ( n == 2 ) type = &IteratorType_Row2D     ;
        else Status_Set ( status, Status_InvalidArgument ) ;
    }
    return type ;
}

Integer Iterator_TypeToInteger ( const IteratorType *type   ,
                                       Status       *status )
{
    Integer n = -1 ;
    if ( ( type != NULL ) && Status_IsOK ( status ) )
    {
             if ( type == &IteratorType_Regular1D ) n = 0 ;
        else if ( type == &IteratorType_RegularND ) n = 1 ;
        else if ( type == &IteratorType_Row2D     ) n = 2 ;
        Status_Set ( status, Status_InvalidArgument ) ;
    }
    return n ;
}
