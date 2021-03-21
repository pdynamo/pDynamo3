# ifndef _DFTGRID
# define _DFTGRID

# include "Coordinates3.h"
# include "DFTGridWeights.h"
# include "GridFunctionDataBlock.h"
# include "IntegerArray1D.h"
# include "List.h"
# include "RealArray1D.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The number of accuracy classes. */
# define NDFTGRID_ACCURACY 5

/* . DFT grid accuracy definitions. */
typedef enum {
    DFTGridAccuracy_VeryLow  = 0 ,
    DFTGridAccuracy_Low      = 1 ,
    DFTGridAccuracy_Medium   = 2 ,
    DFTGridAccuracy_High     = 3 ,
    DFTGridAccuracy_VeryHigh = 4
} DFTGridAccuracy ;

/* . The type for storing grid points. */
typedef struct {
    Integer                atom           ;
    Integer                numberOfPoints ;
    Coordinates3          *coordinates3   ;
    RealArray1D           *weights        ;
    GridFunctionDataBlock *functionData   ;
} DFTGridPointBlock ;

/* . The grid type. */
typedef struct {
    DFTGridAccuracy     accuracy        ;
    Integer             blockSize       ;
    Integer             numberOfPoints  ;
    Integer             numberOfRecords ;
    Real                bfTolerance     ;
    Real                rhoTolerance    ;
    DFTGridPointBlock **records         ;
    DFTGridWeights     *weights         ;
    List               *points          ;
} DFTGrid ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern DFTGrid           *DFTGrid_Allocate               ( const DFTGridAccuracy  accuracy       ,
                                                                 Status          *status         ) ;
extern DFTGrid           *DFTGrid_Construct              ( const DFTGridAccuracy  accuracy       ,
                                                           const IntegerArray1D  *atomicNumbers  ,
                                                           const Coordinates3    *qcCoordinates3 ,
                                                                 Status          *status         ) ;
extern void               DFTGrid_Deallocate             (       DFTGrid        **self           ,
                                                                 Status          *status         ) ;
extern void               DFTGrid_DeallocateFunctionData (       DFTGrid         *self           ,
                                                                 Status          *status         ) ;
extern Integer            DFTGrid_EstimatedPoints        ( const DFTGridAccuracy  accuracy       ,
                                                           const IntegerArray1D  *atomicNumbers  ,
                                                                 Status          *status         ) ;
extern Real               DFTGrid_FunctionByteSize       (       DFTGrid         *self           ,
                                                                 Status          *status         ) ;
extern Boolean            DFTGrid_HasFunctionData        (       DFTGrid         *self           ,
                                                                 Status          *status         ) ;
extern DFTGridPointBlock *DFTGrid_Iterate                (       DFTGrid         *self           ) ;
extern void               DFTGrid_MakeRecords            (       DFTGrid         *self           ,
                                                                 Status          *status         ) ;
extern Integer            DFTGrid_NumberOfFunctionValues (       DFTGrid         *self           ,
                                                                 Status          *status         ) ;
extern Integer            DFTGrid_NumberOfPoints         (       DFTGrid         *self           ) ;
extern Integer            DFTGrid_NumberOfRecords        (       DFTGrid         *self           ) ;

# endif
