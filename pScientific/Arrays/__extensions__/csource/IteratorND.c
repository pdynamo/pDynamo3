/*==================================================================================================================================
! . Regular N-D iterators.
!=================================================================================================================================*/

/*# include <stdio.h>*/

# include "IntegerUtilities.h"
# include "IteratorND.h"
# include "Memory.h"
# include "NumericalMacros.h"

/*==================================================================================================================================
! . Type definitions.
!=================================================================================================================================*/
const IteratorType IteratorType_RegularND = { &IteratorND_Clone          ,
                                              &IteratorND_CurrentIndex   ,
                                              &IteratorND_DataOffSet     ,
                                              &IteratorND_Deallocate     ,
                                              &IteratorND_Dump           ,
                                              &IteratorND_Load           ,
                                              &IteratorND_NextIndex      ,
                                              &IteratorND_NextInnerLoop  ,
                                              &IteratorND_Reset          } ;

/*==================================================================================================================================
! . Regular N-D iterator.
!=================================================================================================================================*/
/*--------------------------------------------------------------------------------------------------------------------------------*/
IteratorND *IteratorND_Allocate ( const Integer rank, Status *status )
{
    IteratorND *self = NULL ;
    if ( ( rank > 0 ) && Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( IteratorND ) ;
        if ( self != NULL )
        {
            self->counters = NULL ;
            self->extents  = NULL ;
            self->strides  = NULL ;
            self->rank     = rank ;
            self->offset   = 0    ;
            self->size     = 0    ;
            if ( rank > 0 )
            {
                self->counters = Integer_Allocate ( rank, status ) ;
                self->extents  = Integer_Allocate ( rank, status ) ;
                self->strides  = Integer_Allocate ( rank, status ) ;
                if ( ( self->counters == NULL ) || ( self->extents == NULL ) || ( self->strides == NULL ) ) IteratorND_Deallocate ( ( void ** ) (&self) ) ;
            }
        }
        if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    }
    return self ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
void *IteratorND_Clone ( void *vSelf, Status *status )
{
    void *vClone = NULL ;
    if ( ( vSelf != NULL ) && Status_IsOK ( status ) )
    {
        IteratorND *self  = ( IteratorND * ) vSelf ;
        IteratorND *clone = IteratorND_Allocate ( self->rank, status ) ;
        IteratorND_Initialize ( clone, self->offset, self->extents, self->strides, status ) ;
        vClone = ( void * ) clone ;
    }
    return vClone ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
Integer IteratorND_CurrentIndex ( void *vSelf )
{
    IteratorND *self = ( IteratorND * ) vSelf ;
    return self->next ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
Integer IteratorND_DataOffSet ( void *vSelf )
{
    Integer n = 0 ;
    if ( vSelf != NULL )
    {
        auto IteratorND *self = ( IteratorND * ) vSelf ;
        n = self->offset ;
    }
    return n ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
void IteratorND_Deallocate ( void **vSelf )
{
    if ( ( vSelf != NULL ) && ( (*vSelf) != NULL ) )
    {
        IteratorND **self = ( IteratorND ** ) vSelf ;
        IteratorND_Finalize ( (*self) ) ;
        Memory_Deallocate ( (*vSelf) ) ;
    }
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
Integer *IteratorND_Dump (       void    *vSelf  ,
                           const Integer  n0     ,
                                 Integer *n      ,
                                 Status  *status )
{
    Integer *state = NULL ;
    if ( ( vSelf != NULL ) && ( n0 >= 0 ) && ( n != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer         m ; 
        auto IteratorND *self = ( IteratorND * ) vSelf ;
        m     = ( n0 + 3 + 2 * self->rank ) ;
        state = Integer_Allocate ( m, status ) ;
        if ( state != NULL )
        {
            auto Integer i, o = n0 ; 
            state[o] = self->offset ; o++ ;
            state[o] = self->rank   ; o++ ;
            state[o] = self->size   ; o++ ;
            for ( i = 0 ; i < self->rank ; i++, o++ ) { state[o] = self->extents[i] ; }
            for ( i = 0 ; i < self->rank ; i++, o++ ) { state[o] = self->strides[i] ; }
        }
        (*n) = m ;
    }
    return state ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
void IteratorND_Finalize ( IteratorND *self )
{
    if ( self != NULL )
    {
        Integer_Deallocate ( &(self->counters) ) ;
        Integer_Deallocate ( &(self->extents ) ) ;
        Integer_Deallocate ( &(self->strides ) ) ;
    }
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
void IteratorND_Initialize ( IteratorND *self, const Integer offset, const Integer *extents, const Integer *strides, Status *status )
{
    if ( self != NULL )
    {
        auto Integer d, extent, n ;
        for ( d = 0, n = 1 ; d < self->rank ; d++ )
        {
            extent = extents[d] ;
            if ( extent < 0 )
            {
                extent = 0 ;
                Status_Set ( status, Status_InvalidArgument ) ;
            }
            self->extents[d] = extent     ;
            self->strides[d] = strides[d] ;
            n *= extent ;
        }
        self->offset = offset ;
        self->size   = n      ;
    }
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
void *IteratorND_Load ( const Integer  n0     ,
                        const Integer  n      ,
                        const Integer *state  ,
                              Status  *status )
{
    void *vSelf = NULL ;
    if ( ( n0 >= 0 ) && ( state != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer rank = state[n0+1] ;
        if ( n != ( n0 + 3 + 2 * rank ) ) Status_Set ( status, Status_InvalidArgument ) ;
        else
        {
            auto IteratorND *self = IteratorND_Allocate ( rank, status ) ;
            if ( self != NULL )
            {
            auto Integer i, o = n0 ; 
                self->offset = state[o] ; o += 2 ;
                self->size   = state[o] ; o++    ;
                for ( i = 0 ; i < self->rank ; i++, o++ ) { self->extents[i] = state[o] ; }
                for ( i = 0 ; i < self->rank ; i++, o++ ) { self->strides[i] = state[o] ; }
                vSelf = ( void * ) self ;
            }
        }
    }
    return vSelf ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
void IteratorND_MakeIterator ( const IteratorND *self, Iterator *iterator )
{
    if ( ( self != NULL ) && ( iterator != NULL ) )
    {
        auto Integer d, innerLoops = 1, rank = self->rank ;
        for ( d = 0 ; d < rank - 1 ; d++ ) innerLoops *= self->extents[d] ;
/*        iterator->readOnly      = False                   ; */
        iterator->extent        = self->extents[rank-1]    ;
        iterator->isRegular     = ( iterator->extent > 1 ) ;
        iterator->numberOfLoops = innerLoops               ;
        iterator->size          = self->size               ;
        iterator->type          = &IteratorType_RegularND  ;
        iterator->vSelf         = ( void * ) self          ;
    }
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
Integer IteratorND_NextIndex ( void *vSelf )
{
    IteratorND *self = ( IteratorND * ) vSelf ;
    Integer next = self->next ;
    if ( next >= 0 )
    {
        auto Integer d, i ;
        for ( d = self->rank-1 ; d >= 0 ; d-- )
        {
            i = self->counters[d] + 1 ;
            if ( i >= self->extents[d] )
            {
                if ( d == 0 ) self->next = -1 ;
                else
                {
                    self->counters[d] = 0 ;
                    self->next       -= ( self->extents[d] - 1 ) * self->strides[d] ;
                }
            }
            else
            {
                self->counters[d] = i ;
                self->next       += self->strides[d] ;
                break ;
            }
        }
    }
    return next ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
Boolean IteratorND_NextInnerLoop ( void *vSelf, Integer *first, Integer *extent, Integer *stride )
{
    IteratorND *self = ( IteratorND * ) vSelf ;
    Boolean isActive = False ;
    Integer next     = self->next ;
    if ( next >= 0 )
    {
        auto Integer d, i ;
        /* . Set the inner loop. */
        d = self->rank-1 ;        
        isActive  = True             ; 
        (*extent) = self->extents[d] ;
        (*first ) = next             ; 
        (*stride) = self->strides[d] ;
        /* . Find the next loop. */
        if ( self->rank > 1 )
        {
            for ( d = self->rank-2 ; d >= 0 ; d-- )
            {
                i = self->counters[d] + 1 ;
                if ( i >= self->extents[d] )
                {
                    if ( d == 0 ) self->next = -1 ;
                    else
                    {
                        self->counters[d] = 0 ;
                        self->next       -= ( self->extents[d] - 1 ) * self->strides[d] ;
                    }
                }
                else
                {
                    self->counters[d] = i ;
                    self->next       += self->strides[d] ;
                    break ;
                }
            }
        }
        else self->next = -1 ;
    }
    return isActive ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
void IteratorND_Reset ( void *vSelf )
{
    Integer d ;
    IteratorND *self = ( IteratorND * ) vSelf ;
    for ( d = 0 ; d < self->rank ; d++ ) self->counters[d] = 0 ;
    if ( self->size <= 0 ) self->next = -1 ;
    else                   self->next =  0 ; /*self->offset ;*/
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
