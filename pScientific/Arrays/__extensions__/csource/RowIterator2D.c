/*==================================================================================================================================
! . Row 2-D iterators.
!=================================================================================================================================*/

# include "IntegerUtilities.h"
# include "Memory.h"
# include "NumericalMacros.h"
# include "RowIterator2D.h"

/*==================================================================================================================================
! . Type definitions.
!=================================================================================================================================*/
const IteratorType IteratorType_Row2D = { &RowIterator2D_Clone         ,
                                          &RowIterator2D_CurrentIndex  ,
                                          &RowIterator2D_DataOffSet    ,
                                          &RowIterator2D_Deallocate    ,
                                          &RowIterator2D_Dump          ,
                                          &RowIterator2D_Load          ,
                                          &RowIterator2D_NextIndex     ,
                                          &RowIterator2D_NextInnerLoop ,
                                          &RowIterator2D_Reset         } ;

/*==================================================================================================================================
! . Regular N-D iterator.
!=================================================================================================================================*/
/*--------------------------------------------------------------------------------------------------------------------------------*/
RowIterator2D *RowIterator2D_Allocate ( const Integer extent0, Status *status )
{
    RowIterator2D *self = NULL ;
    if ( ( extent0 >= 0 ) && Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( RowIterator2D ) ;
        if ( self != NULL )
        {
            self->counter0 = 0 ;
            self->counter1 = 0 ;
            self->extent0  = extent0 ;
            self->extent1  = 0 ;
            self->next     = 0 ;
            self->offset   = 0 ;
            self->size     = 0 ;
            self->stride0  = 0 ;
            self->stride1  = 0 ;
            self->rows     = NULL ;
            if ( extent0 > 0 )
            {
                self->rows = Integer_Allocate ( extent0, status ) ;
                if ( self->rows == NULL ) RowIterator2D_Deallocate ( ( void ** ) (&self) ) ;
            }
        }
        if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    }
    return self ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
void *RowIterator2D_Clone ( void *vSelf, Status *status )
{
    void *vClone = NULL ;
    if ( ( vSelf != NULL ) && Status_IsOK ( status ) )
    {
        RowIterator2D *self  = ( RowIterator2D * ) vSelf ;
        RowIterator2D *clone = RowIterator2D_Allocate ( self->extent0, status ) ;
        RowIterator2D_Initialize ( clone, self->extent1, self->offset, self->stride0, self->stride1, self->rows, status ) ;
        vClone = ( void * ) clone ;
    }
    return vClone ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
Integer RowIterator2D_CurrentIndex ( void *vSelf )
{
    RowIterator2D *self = ( RowIterator2D * ) vSelf ;
    return self->next ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
Integer RowIterator2D_DataOffSet ( void *vSelf )
{
    Integer n = 0 ;
    if ( vSelf != NULL )
    {
        auto RowIterator2D *self = ( RowIterator2D * ) vSelf ;
        n = self->offset ;
    }
    return n ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
void RowIterator2D_Deallocate ( void **vSelf )
{
    if ( ( vSelf != NULL ) && ( (*vSelf) != NULL ) )
    {
        RowIterator2D **self = ( RowIterator2D ** ) vSelf ;
        RowIterator2D_Finalize ( (*self) ) ;
        Memory_Deallocate ( (*vSelf) ) ;
    }
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
Integer *RowIterator2D_Dump (       void    *vSelf  ,
                              const Integer  n0     ,
                                    Integer *n      ,
                                    Status  *status )
{
    Integer *state = NULL ;
    if ( ( vSelf != NULL ) && ( n0 >= 0 ) && ( n != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer        m ; 
        auto RowIterator2D *self = ( RowIterator2D * ) vSelf ;
        m     = ( n0 + 6 + self->extent0 ) ;
        state = Integer_Allocate ( m, status ) ;
        if ( state != NULL )
        {
            auto Integer i, o = n0 ; 
            state[o] = self->extent0 ; o++ ;
            state[o] = self->extent1 ; o++ ;
            state[o] = self->offset  ; o++ ;
            state[o] = self->size    ; o++ ;
            state[o] = self->stride0 ; o++ ;
            state[o] = self->stride1 ; o++ ;
            for ( i = 0 ; i < self->extent0 ; i++, o++ ) { state[o] = self->rows[i] ; }
        }
        (*n) = m ;
    }
    return state ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
void RowIterator2D_Finalize ( RowIterator2D *self )
{
    if ( self != NULL ) { Integer_Deallocate ( &(self->rows) ) ; }
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
void RowIterator2D_Initialize (       RowIterator2D *self    ,
                                const Integer        extent1 ,
                                const Integer        offset  ,
                                const Integer        stride0 ,
                                const Integer        stride1 ,
                                const Integer       *rows    ,
                                      Status        *status  )
{
    if (   ( self != NULL ) && Status_IsOK ( status ) )
    {
        if ( ( extent1 < 0 ) || ( ( rows == NULL ) && ( self->extent0 != 0 ) ) ) Status_Set ( status, Status_InvalidArgument ) ;
        else
        {
            self->extent1 = extent1 ;
            self->offset  = offset  ;
            self->size    = self->extent0 * self->extent1 ;
            self->stride0 = stride0 ;
            self->stride1 = stride1 ;
            Integer_CopyTo ( rows, self->extent0, self->rows, status ) ;
        }
    }
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
void *RowIterator2D_Load ( const Integer  n0     ,
                           const Integer  n      ,
                           const Integer *state  ,
                                 Status  *status )
{
    void *vSelf = NULL ;
    if ( ( n0 >= 0 ) && ( state != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer extent0 = state[n0] ;
        if ( n != ( n0 + 6 + extent0 ) ) Status_Set ( status, Status_InvalidArgument ) ;
        else
        {
            auto RowIterator2D *self = RowIterator2D_Allocate ( extent0, status ) ;
            if ( self != NULL )
            {
                auto Integer i, o = n0+1 ;
                self->extent1 = state[o] ; o++ ;
                self->offset  = state[o] ; o++ ;
                self->size    = state[o] ; o++ ;
                self->stride0 = state[o] ; o++ ;
                self->stride1 = state[o] ; o++ ;
                for ( i = 0 ; i < extent0 ; i++, o++ ) { self->rows[i] = state[o] ; }
                vSelf = ( void * ) self ;
            }
        }
    }
    return vSelf ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
void RowIterator2D_MakeIterator ( const RowIterator2D *self, Iterator *iterator )
{
    if ( ( self != NULL ) && ( iterator != NULL ) )
    {
/*        iterator->readOnly      = False                    ; */
        iterator->extent        = self->extent1            ;
        iterator->isRegular     = ( iterator->extent > 1 ) ;
        iterator->numberOfLoops = self->extent0            ;
        iterator->size          = self->size               ;
        iterator->type          = &IteratorType_Row2D      ;
        iterator->vSelf         = ( void * ) self          ;
    }
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
Integer RowIterator2D_NextIndex ( void *vSelf )
{
    RowIterator2D *self = ( RowIterator2D * ) vSelf ;
    Integer next = self->next ;
    if ( next >= 0 )
    {
        auto Integer i ;
        /* . Columns then rows. */
        i = self->counter1 + 1 ;
        if ( i >= self->extent1 )
        {
            self->counter1 = 0 ;
            i = self->counter0 + 1 ;
            if ( i >= self->extent0 ) { self->next = -1 ; }
            else
            {
                self->counter0 = i ;
                self->next     = self->rows[i] * self->stride0 ; /*+ self->offset ;*/
            }
        }
        else
        {
            self->counter1 = i ;
            self->next    += self->stride1 ;
        }
    }
    return next ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
Boolean RowIterator2D_NextInnerLoop ( void *vSelf, Integer *first, Integer *extent, Integer *stride )
{
    RowIterator2D *self = ( RowIterator2D * ) vSelf ;
    Boolean isActive = False ;
    Integer next     = self->next ;
    if ( next >= 0 )
    {
        auto Integer i ;
        /* . Set the inner loop. */
        isActive  = True          ; 
        (*extent) = self->extent1 ;
        (*first ) = next          ; 
        (*stride) = self->stride1 ;
        /* . Find the next loop. */
        i = self->counter0 + 1 ;
        if ( i >= self->extent0 ) { self->next = -1 ; }
        else
        {
            self->counter0 = i ;
            self->next     = self->rows[i] * self->stride0 ; /*+ self->offset ;*/
        }
    }
    return isActive ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
void RowIterator2D_Reset ( void *vSelf )
{
    RowIterator2D *self = ( RowIterator2D * ) vSelf ;
    self->counter0 = 0 ;
    self->counter1 = 0 ;
    if ( self->size <= 0 ) self->next = -1 ;
    else                   self->next = self->rows[0] * self->stride0 ; /*+ self->offset ;*/
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
