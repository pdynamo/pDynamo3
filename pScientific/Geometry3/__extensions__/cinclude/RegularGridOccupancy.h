# ifndef _REGULARGRIDOCCUPANCY
# define _REGULARGRIDOCCUPANCY

# include "Integer.h"
# include "IntegerArray1D.h"
# include "RealArray2D.h"
# include "RegularGrid.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The occupancy type. */
typedef struct {
    Integer         numberOfCells   ;
    Integer         numberOfPoints  ;
    IntegerArray1D *cellFirstPoints ; /* . The indices of the first points in each cell in the cellPoints array. */
    IntegerArray1D *cellPoints      ; /* . The indices of the points within each cell. */
    IntegerArray1D *cellTotalPoints ; /* . The total number of points per cell. */
    IntegerArray1D *pointCells      ; /* . The grid cell within which a point is found. */
} RegularGridOccupancy ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern RegularGridOccupancy *RegularGridOccupancy_Allocate          ( const Integer                numberOfCells  ,
                                                                      const Integer                numberOfPoints ,
                                                                            Status                *status         ) ;
extern void                  RegularGridOccupancy_Deallocate        (       RegularGridOccupancy **self           ) ;
extern void                  RegularGridOccupancy_DeallocateVoid    (       void                 **vSelf          ) ;
extern Integer               RegularGridOccupancy_Fill              (       RegularGridOccupancy  *self           ,
                                                                      const RegularGrid           *grid           ,
                                                                      const RealArray2D           *points         ,
                                                                            Status                *status         ) ;
extern RegularGridOccupancy *RegularGridOccupancy_FromGridAndPoints ( const RegularGrid           *grid           ,
                                                                      const RealArray2D           *points         ,
                                                                            Status                *status         ) ;
extern void                  RegularGridOccupancy_Initialize        (       RegularGridOccupancy  *self           ) ;

# endif
