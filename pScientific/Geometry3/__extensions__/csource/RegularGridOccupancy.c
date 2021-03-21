/*==================================================================================================================================
! . Regular grid occupancy procedures.
!=================================================================================================================================*/

# include <stdio.h>

# include "Memory.h"
# include "RegularGridOccupancy.h"

/*# define DEBUGPRINTING*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
RegularGridOccupancy *RegularGridOccupancy_Allocate ( const Integer numberOfCells, const Integer numberOfPoints, Status *status )
{
    RegularGridOccupancy *self = NULL ;
# ifdef DEBUGPRINTING
printf ( "\nRegular Grid Occupancy Allocate:\n" ) ;
printf ( "Number Of Cells      = %6d\n", numberOfCells  ) ;
printf ( "Number Of Points     = %6d\n", numberOfPoints ) ;
# endif
    if ( ( numberOfCells > 0 ) && ( numberOfPoints > 0 ) )
    {
        self = Memory_AllocateType ( RegularGridOccupancy ) ;
        if ( self != NULL )
        {
            /* . Initialization. */
            self->numberOfCells   = numberOfCells  ;
            self->numberOfPoints  = numberOfPoints ;
            self->cellFirstPoints = NULL ;
            self->cellPoints      = NULL ;
            self->cellTotalPoints = NULL ;
            self->pointCells      = NULL ;
            /* . Array allocation. */
            self->cellFirstPoints = IntegerArray1D_AllocateWithExtent ( numberOfCells , NULL ) ;
            self->cellPoints      = IntegerArray1D_AllocateWithExtent ( numberOfPoints, NULL ) ;
            self->cellTotalPoints = IntegerArray1D_AllocateWithExtent ( numberOfCells , NULL ) ;
            self->pointCells      = IntegerArray1D_AllocateWithExtent ( numberOfPoints, NULL ) ;
            /* . Check. */
            if ( ( self->cellFirstPoints == NULL ) || ( self->cellPoints == NULL ) ||
                 ( self->cellTotalPoints == NULL ) || ( self->pointCells == NULL ) ) RegularGridOccupancy_Deallocate ( &self ) ;
            /* . Array initialization. */
            RegularGridOccupancy_Initialize ( self ) ;
        }
        if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    }
    else Status_Set ( status, Status_InvalidArgument ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RegularGridOccupancy_Deallocate ( RegularGridOccupancy **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        IntegerArray1D_Deallocate ( &((*self)->cellFirstPoints) ) ;
        IntegerArray1D_Deallocate ( &((*self)->cellPoints     ) ) ;
        IntegerArray1D_Deallocate ( &((*self)->cellTotalPoints) ) ;
        IntegerArray1D_Deallocate ( &((*self)->pointCells     ) ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

void RegularGridOccupancy_DeallocateVoid ( void **vSelf ) { RegularGridOccupancy **self = ( RegularGridOccupancy ** ) vSelf ; RegularGridOccupancy_Deallocate ( self ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Fill the structure given a grid and a set of points.
! . Note that points that are not on the grid are ignored.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer RegularGridOccupancy_Fill ( RegularGridOccupancy *self, const RegularGrid *grid, const RealArray2D *points, Status *status )
{
    Integer totalOffGrid = -1 ;
    if ( ( self != NULL ) && ( grid != NULL ) && ( points != NULL ) )
    {
# ifdef DEBUGPRINTING
printf ( "\nRegular Grid Occupancy Fill Start:\n" ) ;
printf ( "Number Of Cells      = %6d (%6d)\n", self->numberOfCells  , RegularGrid_NumberOfGridPoints ( grid   ) ) ;
printf ( "Number Of Dimensions = %6d (%6d)\n", grid->ndimensions    , View2D_Columns            ( points ) ) ;
printf ( "Number Of Points     = %6d (%6d)\n", self->numberOfPoints , View2D_Rows               ( points ) ) ;
# endif
        if ( ( self->numberOfCells  == RegularGrid_NumberOfGridPoints ( grid   ) ) &&
             ( self->numberOfPoints == View2D_Rows                    ( points ) ) &&
             ( grid->ndimensions    == View2D_Columns                 ( points ) ) )
        {
            auto Integer c, n, p ;

            /* . Determine the cell position of each point and the total number of points in each cell. */
            for ( p = totalOffGrid = 0 ; p < self->numberOfPoints ; p++ )
            {
                c = RegularGrid_FindCellIDOfPoint ( grid, Array2D_RowPointer ( points, p ) ) ;
                if ( c >= 0 )
                {
                    Array1D_Item ( self->cellTotalPoints , c ) += 1 ;
                    Array1D_Item ( self->pointCells      , p )  = c ;
                }
                else totalOffGrid++ ;
            }

            /* . Determine the index of the first point for each cell. */
            for ( c = n = 0 ; c < self->numberOfCells ; c++ )
            {
                Array1D_Item ( self->cellFirstPoints, c ) = n ;
                n += Array1D_Item ( self->cellTotalPoints, c ) ;
            }

            /* . Fill the cell points index array. */
            for ( p = 0 ; p < self->numberOfPoints ; p++ )
            {
                c = Array1D_Item ( self->pointCells , p ) ;
                if ( c >= 0 )
                {
                    n = Array1D_Item ( self->cellFirstPoints , c ) ;
                    Array1D_Item ( self->cellPoints          , n )  = p ;
                    Array1D_Item ( self->cellFirstPoints     , c ) += 1 ;
                }
            }

            /* . Reset cellFirstPoints. */
            IntegerArray1D_Add ( self->cellFirstPoints, -1, self->cellTotalPoints, NULL ) ;
# ifdef DEBUGPRINTING
printf ( "\nRegular Grid Occupancy Fill Stop:\n" ) ;
printf ( "Number Of Cells      = %6d\n", self->numberOfCells  ) ;
printf ( "Number Of Dimensions = %6d\n", grid->ndimensions    ) ;
printf ( "Number Of Points     = %6d\n", self->numberOfPoints ) ;
printf ( "Total Off Grid       = %6d\n", totalOffGrid         ) ;
printf ( "Total On  Grid       = %6d\n", IntegerArray1D_Sum ( self->cellTotalPoints ) ) ;
printf ( "cellFirstPoints\n" ) ;
IntegerArray1D_Print ( self->cellFirstPoints  ) ;
printf ( "\ncellPoints\n" ) ;
IntegerArray1D_Print ( self->cellPoints       ) ;  
printf ( "\ncellTotalPoints\n" ) ;
IntegerArray1D_Print ( self->cellTotalPoints  ) ;
printf ( "\npointCells\n" ) ;
IntegerArray1D_Print ( self->pointCells       ) ;
# endif
        }
        else Status_Set ( status, Status_InvalidArgument ) ;
    }
    return totalOffGrid ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Constructor given a grid and associated points.
!---------------------------------------------------------------------------------------------------------------------------------*/
RegularGridOccupancy *RegularGridOccupancy_FromGridAndPoints ( const RegularGrid *grid, const RealArray2D *points, Status *status )
{
    RegularGridOccupancy *self = NULL ;
    if ( ( grid != NULL ) && ( points != NULL ) )
    {
        self = RegularGridOccupancy_Allocate ( RegularGrid_NumberOfGridPoints ( grid ), View2D_Rows ( points ), status ) ;
        RegularGridOccupancy_Fill ( self, grid, points, status ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RegularGridOccupancy_Initialize ( RegularGridOccupancy *self )
{
    if ( self != NULL )
    {
        IntegerArray1D_Set ( self->cellFirstPoints ,  0 ) ;
        IntegerArray1D_Set ( self->cellPoints      , -1 ) ;
        IntegerArray1D_Set ( self->cellTotalPoints ,  0 ) ;
        IntegerArray1D_Set ( self->pointCells      , -1 ) ;
    }
}
