/*==================================================================================================================================
! . Slice handling.
!=================================================================================================================================*/

# include "Memory.h"
# include "NumericalMacros.h"
# include "Slice.h"

/*==================================================================================================================================
! . Local functions.
!=================================================================================================================================*/
static void Slice_Initialize      (       Slice   *self     ) ;
static void Slice_SetSliceIndices (       Slice   *self     ,
                                    const Boolean  isScalar ,
                                    const Integer  start    ,
                                    const Integer  stop     ,
                                    const Integer  stride   ,
                                    const Integer  n        ) ;

/*==================================================================================================================================
! . A multislice.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
MultiSlice *MultiSlice_Allocate ( const Integer capacity, Status *status )
{
    MultiSlice *self = NULL ;
    if ( ( capacity > 0 ) && Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( MultiSlice ) ;
        if ( self != NULL )
        {
            self->capacity = 0    ;
            self->items    = NULL ;
            self->rank     = 0    ;
            self->items    = Memory_AllocateArrayOfTypes ( capacity, Slice ) ;
            if ( self->items == NULL ) MultiSlice_Deallocate ( &self ) ;
            else
            {
                self->capacity = capacity ;
                auto Integer i ;
                for ( i = 0 ; i < capacity ; i++ ) Slice_Initialize ( &(self->items[i]) ) ;
            }
        }
        if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MultiSlice_Deallocate ( MultiSlice **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        Memory_Deallocate ( (*self)->items ) ;
        Memory_Deallocate ( (*self)        ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the capacity of the multislice.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer MultiSlice_GetCapacity ( const MultiSlice *self ) { return ( (self) == NULL ? 0 : (self)->capacity ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the extent of a dimension of the multislice.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer MultiSlice_GetExtent ( const MultiSlice *self, const Integer dimension, Status *status )
{
    Integer extent = -1 ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        if ( dimension < self->capacity ) extent = self->items[dimension].extent ;
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
    return extent ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the rank of the multislice.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer MultiSlice_GetRank ( const MultiSlice *self ) { return ( (self) == NULL ? 0 : (self)->rank ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the size of the multislice (scalars are rank 0 but have size 1).
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer MultiSlice_GetSize ( const MultiSlice *self )
{
    Integer size = 0 ;
    if ( self != NULL )
    {
        auto Integer i ;
        for ( i = 0, size = 1 ; i < self->capacity ; i++ ) size *= self->items[i].extent ;
    }
    return size ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set the rank from the constituent slices.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MultiSlice_SetRank ( MultiSlice *self )
{
    if ( self != NULL )
    {
        auto Integer i, rank = 0 ;
        for ( i = 0 ; i < self->capacity ; i++ ) { if ( ! self->items[i].isScalar ) rank += 1 ; }
        self->rank = rank ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set the data for a slice.
! . All slice arguments are assumed to have been checked beforehand ...
!---------------------------------------------------------------------------------------------------------------------------------*/
void MultiSlice_SetSliceIndices (       MultiSlice *self      ,
                                  const Integer     dimension ,
                                  const Boolean     isScalar  ,
                                  const Integer     start     ,
                                  const Integer     stop      ,
                                  const Integer     stride    ,
                                  const Integer     n         )
{
    if ( ( self != NULL ) && ( Integer_IsInRange ( dimension, 0, self->capacity ) ) )
    {
        Slice_SetSliceIndices ( &(self->items[dimension]), isScalar, start, stop, stride, n ) ;
    }
}

/*==================================================================================================================================
! . A single slice.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void Slice_Initialize ( Slice *self )
{
    if ( self != NULL )
    {
        self->isScalar = False ;
        self->extent   =  0 ;
        self->start    =  0 ;
        self->stop     = -1 ;
        self->stride   =  1 ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set a slice from input data - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void Slice_SetSliceIndices (       Slice   *self     ,
                                    const Boolean  isScalar ,
                                    const Integer  start    ,
                                    const Integer  stop     ,
                                    const Integer  stride   ,
                                    const Integer  n        )
{
    if ( self != NULL )
    {
        if ( isScalar )
        {
            self->isScalar = True      ;
            self->extent   = 1         ;
            self->start    = start     ;
            self->stop     = start + 1 ;
            self->stride   = 1         ;
        }
        else
        {
            self->isScalar = False     ;
            self->extent   = n         ;
            self->start    = start     ;
            self->stop     = stop      ;
            self->stride   = stride    ;
        }
    }
}

/*==================================================================================================================================
! . Slice index checking functions.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Check the scalar index.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer SliceIndices_CheckScalar ( const Integer  index  ,
                                   const Integer  extent ,
                                         Status  *status )
{
    Integer start = index ;
    if ( Status_IsOK ( status ) )
    {
        if ( start < 0 ) start += extent ;
        if ( ! Integer_IsInRange ( start, 0, extent ) ) Status_Set ( status, Status_IndexOutOfRange ) ;
    }
    return start ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check slice indices given an extent.
! . It seems that Python only allows stop to be -1 when it has the value None (pStop == NULL).
! . However, when stop is explicitly specified it's lower limit is 0.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer SliceIndices_CheckSlice  ( const Integer  *pStart  ,
                                   const Integer  *pStop   ,
                                   const Integer  *pStride ,
                                   const Integer   extent  ,
                                         Integer  *qStart  ,
                                         Integer  *qStop   ,
                                         Integer  *qStride ,
                                         Status   *status  )
{
    Integer n = 0 ;
    if ( Status_IsOK ( status ) )
    {
        auto Integer start = 0, stop = 0, stride = 1 ;
        if ( extent > 0 )
        {
            if ( pStride == NULL ) { stride = 1 ; } else { stride = (*pStride) ; }
            if ( stride == 0 ) Status_Set ( status, Status_IndexOutOfRange ) ;
            else
            {
                if ( pStart == NULL ) { if ( stride > 0 ) { start = 0 ; } else { start = extent-1 ; } }
                else
                {
                    start = (*pStart) ;
                    while ( start < 0 ) start += extent ;
                    if ( stride > 0 ) start = Minimum ( start, extent   ) ;
                    else              start = Minimum ( start, extent-1 ) ;
                }
                if ( pStop == NULL ) { if ( stride > 0 ) { stop = extent ; } else { stop = -1 ; } }
                else
                {
                    stop = (*pStop) ;
                    while ( stop < 0 ) stop += extent ;
                    if ( stride > 0 ) stop = Minimum ( stop, extent   ) ;
                    else              stop = Minimum ( stop, extent-1 ) ;
                }
                if      ( ( stride > 0 ) && ( stop > start ) ) n = ( stop - start - 1 ) / stride + 1 ;
                else if ( ( stride < 0 ) && ( stop < start ) ) n = ( stop - start + 1 ) / stride + 1 ;
            }
        }
        if ( qStart  != NULL ) (*qStart ) = start  ;
        if ( qStop   != NULL ) (*qStop  ) = stop   ;
        if ( qStride != NULL ) (*qStride) = stride ;
    }
    return n ;
}
