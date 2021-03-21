/*==================================================================================================================================
! . Regular 1-D iterators.
!=================================================================================================================================*/

# include "IntegerUtilities.h"
# include "Memory.h"
# include "NumericalMacros.h"
# include "Iterator1D.h"

/*==================================================================================================================================
! . Type definitions.
!=================================================================================================================================*/
const IteratorType IteratorType_Regular1D = { &Iterator1D_Clone          ,
                                              &Iterator1D_CurrentIndex   ,
                                              &Iterator1D_DataOffSet     ,
                                              &Iterator1D_Deallocate     ,
                                              &Iterator1D_Dump           ,
                                              &Iterator1D_Load           ,
                                              &Iterator1D_NextIndex      ,
                                              &Iterator1D_NextInnerLoop  ,
                                              &Iterator1D_Reset          } ;

/*==================================================================================================================================
! . Regular 1-D iterator.
!=================================================================================================================================*/
/*--------------------------------------------------------------------------------------------------------------------------------*/
Iterator1D *Iterator1D_Allocate ( Status *status )
{
    Iterator1D *self = NULL ;
    self = Memory_AllocateType ( Iterator1D ) ;
    if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    return self ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
void *Iterator1D_Clone ( void *vSelf, Status *status )
{
    void *vClone = NULL ;
    if ( vSelf != NULL )
    {
        Iterator1D *self  = ( Iterator1D * ) vSelf ;
        Iterator1D *clone = Iterator1D_Allocate ( status ) ;
        Iterator1D_Initialize ( clone, self->offset, self->extent, self->stride, status ) ;
        vClone = ( void * ) clone ;
    }
    return vClone ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
Integer Iterator1D_CurrentIndex ( void *vSelf )
{
    Iterator1D *self = ( Iterator1D * ) vSelf ;
    return self->next ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
Integer Iterator1D_DataOffSet ( void *vSelf )
{
    Integer n = 0 ;
    if ( vSelf != NULL )
    {
        Iterator1D *self  = ( Iterator1D * ) vSelf ;
        n = self->offset ;
    }
    return n ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
void Iterator1D_Deallocate ( void **vSelf )
{
    if ( ( vSelf != NULL ) && ( (*vSelf) != NULL ) ) Memory_Deallocate ( (*vSelf) ) ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
Integer *Iterator1D_Dump (       void    *vSelf  ,
                           const Integer  n0     ,
                                 Integer *n      ,
                                 Status  *status )
{
    Integer *state = NULL ;
    if ( ( vSelf != NULL ) && ( n0 >= 0 ) && ( n != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer         m ; 
        auto Iterator1D *self = ( Iterator1D * ) vSelf ;
        m     = ( n0 + 3 ) ;
        state = Integer_Allocate ( m, status ) ;
        if ( state != NULL )
        {
            state[n0  ] = self->extent ;
            state[n0+1] = self->offset ;
            state[n0+2] = self->stride ;
        }
        (*n) = m ;
    }
    return state ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
void Iterator1D_Initialize ( Iterator1D *self, const Integer offset, const Integer extent, const Integer stride, Status *status )
{
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        if ( extent < 0 ) Status_Set ( status, Status_InvalidArgument ) ;
        else
        {
            self->extent = extent ;
            self->offset = offset ;
            self->stride = stride ;
        }
    }
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
void *Iterator1D_Load ( const Integer  n0     ,
                        const Integer  n      ,
                        const Integer *state  ,
                              Status  *status )
{
    void *vSelf = NULL ;
    if ( ( n0 >= 0 ) && ( state != NULL ) && Status_IsOK ( status ) )
    {
        if ( n != ( n0 + 3 ) ) Status_Set ( status, Status_InvalidArgument ) ;
        else
        {
            auto Iterator1D *self = Iterator1D_Allocate ( status ) ;
            if ( self != NULL )
            {
                self->extent = state[n0  ] ;
                self->offset = state[n0+1] ;
                self->stride = state[n0+2] ;
                vSelf = ( void * ) self ;
            }
        }
    }
    return vSelf ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
void Iterator1D_MakeIterator ( const Iterator1D *self, Iterator *iterator )
{
    if ( ( self != NULL ) && ( iterator != NULL ) )
    {
/*        iterator->readOnly      = False                   ; */
        iterator->extent        = self->extent            ;
        iterator->isRegular     = ( self->extent > 1 )    ;
        iterator->numberOfLoops = 1                       ;
        iterator->size          = self->extent            ;
        iterator->type          = &IteratorType_Regular1D ;
        iterator->vSelf         = ( void * ) self         ;
    }
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
Integer Iterator1D_NextIndex ( void *vSelf )
{
    Iterator1D *self = ( Iterator1D * ) vSelf ;
    Integer     next = self->next ;
    if ( next >= 0 )
    {
        self->counter += 1 ;
        if ( self->counter >= self->extent ) self->next  = -1 ;
        else                                 self->next += self->stride ;
    }
    return next ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
Boolean Iterator1D_NextInnerLoop ( void *vSelf, Integer *first, Integer *extent, Integer *stride )
{
    Iterator1D *self = ( Iterator1D * ) vSelf ;
    Boolean isActive = False ;
    Integer next     = self->next ;
    if ( next >= 0 )
    {
        isActive   = True         ; 
        (*extent)  = self->extent ;
        (*first )  = next         ; 
        (*stride)  = self->stride ;
        self->next = -1           ;
    }
    return isActive ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
void Iterator1D_Reset ( void *vSelf )
{
    Iterator1D *self = ( Iterator1D * ) vSelf ;
    self->counter = 0 ;
    if ( self->extent <= 0 ) self->next = -1 ;
    else                     self->next =  0 ; /*self->offset ;*/
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
