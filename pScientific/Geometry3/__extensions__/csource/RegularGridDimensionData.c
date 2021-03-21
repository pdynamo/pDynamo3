/*==================================================================================================================================
! . 1-D indexed value arrays.
!=================================================================================================================================*/

# include <stdio.h>
# include <stdlib.h>

# include "Memory.h"
# include "NumericalMacros.h"
# include "RegularGridDimensionData.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Integer Value_CompareAscending  ( const void *vterm1, const void *vterm2 ) ;
static Integer Value_CompareDescending ( const void *vterm1, const void *vterm2 ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
RegularGridDimensionData *RegularGridDimensionData_Allocate ( const Integer extent, Status *status )
{
    RegularGridDimensionData *self = Memory_AllocateType ( RegularGridDimensionData ) ;
    if ( self != NULL )
    {
        self->items   = NULL   ;
        self->isOwner = True   ;
        self->isSlice = False  ;
        self->extent  = Maximum ( extent, 0 ) ;
        self->offset  = 0      ;
        self->size    = Maximum ( extent, 0 ) ;
        self->stride  = 1      ;
        if ( extent > 0 )
        {
            self->items = Memory_AllocateArrayOfTypes ( extent, RegularGridDimensionDatum ) ;
            if ( self->items == NULL ) RegularGridDimensionData_Deallocate ( &self ) ;
        }
        else if ( extent < 0 ) Status_Set ( status, Status_InvalidArgument ) ;
    }
    if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check the data indices to ensure that they fall within a specific range.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer RegularGridDimensionData_CheckIndices ( RegularGridDimensionData *self, const Integer range, const Boolean isPeriodic )
{
    Integer numberOutsideRange = 0 ;
    if ( self != NULL )
    {
        auto Integer i, index ;
        if ( isPeriodic )
        {
            for ( i = 0 ; i < self->extent ; i++ )
            {
                index = self->items[i].index ;
                if ( index < 0 )
                {
                    while ( index <  0     ) index += range ;
                    self->items[i].index = index ;
                }
                else if ( index >= range )
                {
                    while ( index >= range ) index -= range ;
                    self->items[i].index = index ;
                }
            }
        }
        else
        {
            for ( i = 0 ; i < self->extent ; i++ )
            {
                index = self->items[i].index ;
                if ( index < 0 )
                {
                    numberOutsideRange += 1 ;
                    self->items[i].index = 0 ;
                }
                else if ( index >= range )
                {
                    numberOutsideRange += 1 ;
                    self->items[i].index = range - 1 ;
                }
            }
        }
    }
    return numberOutsideRange ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RegularGridDimensionData_Deallocate ( RegularGridDimensionData **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        if ( (*self)->isOwner ) Memory_Deallocate ( (*self)->items ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Length.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer RegularGridDimensionData_Length ( const RegularGridDimensionData *self )
{
    if ( self == NULL ) return 0 ;
    else                return self->extent ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Sorting by value - ascending in-place.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RegularGridDimensionData_SortByValue ( RegularGridDimensionData *self, const Boolean doAscending )
{
    if ( ( self != NULL ) && ( self->extent > 1 ) )
    {
        if ( doAscending ) qsort ( ( void * ) RegularGridDimensionData_Items ( self ), ( size_t ) self->extent, self->stride * sizeof ( RegularGridDimensionDatum ), ( void * ) Value_CompareAscending  ) ;
        else               qsort ( ( void * ) RegularGridDimensionData_Items ( self ), ( size_t ) self->extent, self->stride * sizeof ( RegularGridDimensionDatum ), ( void * ) Value_CompareDescending ) ;
    }
}

/*==================================================================================================================================
! . Private procedures.
!=================================================================================================================================*/
static Integer Value_CompareAscending ( const void *vterm1, const void *vterm2 )
{
    Integer i ;
    RegularGridDimensionDatum *term1, *term2 ;
    term1 = ( RegularGridDimensionDatum * ) vterm1 ;
    term2 = ( RegularGridDimensionDatum * ) vterm2 ;
         if ( term1->value < term2->value ) i = -1 ;
    else if ( term1->value > term2->value ) i =  1 ;
    else i = 0 ;
    return i ;
}

static Integer Value_CompareDescending ( const void *vterm1, const void *vterm2 )
{
    Integer i ;
    RegularGridDimensionDatum *term1, *term2 ;
    term1 = ( RegularGridDimensionDatum * ) vterm1 ;
    term2 = ( RegularGridDimensionDatum * ) vterm2 ;
         if ( term1->value < term2->value ) i =  1 ;
    else if ( term1->value > term2->value ) i = -1 ;
    else i = 0 ;
    return i ;
}
