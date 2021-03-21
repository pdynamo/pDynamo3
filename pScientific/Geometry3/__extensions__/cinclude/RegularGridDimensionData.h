# ifndef _REGULARGRIDDIMENSIONDATA
# define _REGULARGRIDDIMENSIONDATA

# include "Boolean.h"
# include "Integer.h"
# include "Real.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The element type. */
typedef struct {
    Integer index ;
    Real    value ;
} RegularGridDimensionDatum ;

/* . The array type. */
typedef struct {
    Boolean  isOwner ;
    Boolean  isSlice ;
    Integer  extent  ;
    Integer  offset  ;
    Integer  size    ;
    Integer  stride  ;
    RegularGridDimensionDatum *items ;
} RegularGridDimensionData ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . A pointer to the start of the data. */
# define RegularGridDimensionData_Items( self ) ( &((self)->items[(self)->offset]) )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern RegularGridDimensionData *RegularGridDimensionData_Allocate     ( const Integer                    extent      ,
                                                                               Status                    *status      ) ;
extern void                      RegularGridDimensionData_Deallocate   (       RegularGridDimensionData **self        ) ;
extern Integer                   RegularGridDimensionData_CheckIndices (       RegularGridDimensionData  *self        ,
                                                                         const Integer                    range       ,
                                                                         const Boolean                    isPeriodic  ) ;
extern Integer                   RegularGridDimensionData_Length       ( const RegularGridDimensionData  *self        ) ;
extern void                      RegularGridDimensionData_SortByValue  (       RegularGridDimensionData  *self        ,
                                                                         const Boolean                    doAscending ) ;

# endif
