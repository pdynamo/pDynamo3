# ifndef _REGULARGRID
# define _REGULARGRID

# include "Boolean.h"
# include "BooleanArray1D.h"
# include "Integer.h"
# include "IntegerArray1D.h"
# include "IntegerArray2D.h"
# include "NumericalMacros.h"
# include "Real.h"
# include "RealArrayND.h"
# include "Status.h"
# include "Vector3.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The regular grid dimension type. */
typedef struct {
    Boolean isPeriodic    ;
    Integer bins          ;
    Integer stride        ;
    Real    binSize       ;
    Real    lower         ;
    Real    midPointLower ;
    Real    period        ;
    Real    upper         ;
} RegularGridDimension ;

/* . The regular grid search range type. */
typedef struct {
    Integer         numberOfCells       ;
    Integer         numberOfCells0      ;
    Real            cutOff              ;
    Real            cutOffSquared       ;
    BooleanArray1D *isFullyWithinRange  ;
    BooleanArray1D *isFullyWithinRange0 ;
    Integer        *workI               ;
    IntegerArray1D *cellIDs             ;
    IntegerArray2D *cellIndices0        ;
} RegularGridSearchRange ;

/* . The regular grid type. */
typedef struct {
    Integer               ndimensions ;
    Integer              *workI       ;
    Real                 *workR       ;
    RegularGridDimension *dimensions  ;
} RegularGrid ;


/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Adjust a distance to take account of periodicity. */
# define RegularGridDimension_AdjustPeriodicDistance( self, x ) { x -= ( Real_RoundToInteger ( x / self.period ) ) * self.period ; }

/* . Find the number of bins within a cutOff. */
/* . This number cannot be larger than the number of bins for a periodic dimension. */
# define RegularGridDimension_BinsWithinCutOff( self, cutOff, l, u ) \
    { \
        u = ( Integer ) ceil ( cutOff / self.binSize ) ; \
        if ( ( self.isPeriodic ) && ( ( 2 * u + 1 ) > self.bins ) ) { l = 0 ; u = self.bins - 1 ; } \
        else l = - u ; \
    }

/* . Find the integer bin index for a real datum for the dimension. */
# define RegularGridDimension_FindBinIndex( self, x, f, i ) \
    { \
        f = ( x - self.lower ) / self.binSize ; \
        i = ( Integer ) floor ( f ) ; \
    }

/* . Regularize the index for a dimension. */
# define RegularGridDimension_RegularizeIndex( self, i, isOK ) \
    { \
        isOK = True ; \
        if ( self.isPeriodic ) i = Modulo ( i, self.bins ) ; \
        else if ( ( i < 0 ) || ( i >= self.bins ) ) isOK = False ; \
    }

/* . Regularize the index for a dimension, making it negative if it is invalid. */
# define RegularGridDimension_RegularizeIndexToMinusOne( self, i ) \
    { \
        if ( self.isPeriodic ) i = Modulo ( i, self.bins ) ; \
        else if ( ( i < 0 ) || ( i >= self.bins ) ) i = -1 ; \
    }

/* . Regularize the index for a dimension. . */
# define RegularGridDimension_RegularizeIndexToBoundary( self, i, isOK ) \
    { \
        isOK = True ; \
        if ( self.isPeriodic ) i = Modulo ( i, self.bins ) ; \
        else if ( i <  0         ) { i = 0             ; isOK = False ; } \
        else if ( i >= self.bins ) { i = self.bins - 1 ; isOK = False ; } \
    }

/* . Determine a cellID in 3D. */
# define RegularGrid_CellIndicesToID3D( self, i, j, k ) ( (i)*(self)->dimensions[0].stride+(j)*(self)->dimensions[1].stride+(k)*(self)->dimensions[2].stride )

/* . Copy data for a dimension between source and destination. */
# define RegularGrid_CopyDimensionData( self, source, destination ) \
    { \
        auto Integer d ; \
        for ( d = 0 ; d < self->ndimensions ; d++ ) destination[d] = source[d] ; \
    }

/* . Decompose a point with results stored internally. */
# define RegularGrid_DecomposePoint( self, coordinates ) \
    { \
        auto Integer d ; \
        for ( d = 0 ; d < self->ndimensions ; d++ ) RegularGridDimension_FindBinIndex ( self->dimensions[d], coordinates[d], self->workR[d], self->workI[d] ) ; \
    }

/* . Determine a cellID from internal indices. */
# define RegularGrid_MakeCellID( self, c ) \
    { \
        auto Integer d ; \
        for ( c = d = 0 ; d < self->ndimensions ; d++ ) c += self->dimensions[d].stride * self->workI[d] ; \
    }

/* . Regularize internally stored indices. */
# define RegularGrid_RegularizeIndices( self, outsideGrid ) \
    { \
        auto Boolean isOK ; \
        auto Integer d    ; \
        for ( d = outsideGrid = 0 ; d < self->ndimensions ; d++ ) \
        { \
            RegularGridDimension_RegularizeIndex ( self->dimensions[d], self->workI[d], isOK ) ; \
            if ( ! isOK ) outsideGrid += 1 ; \
        } \
    }

/* . Regularize internally stored indices with indices outside the grid forced to -1. */
# define RegularGrid_RegularizeIndicesToMinusOne( self, outsideGrid ) \
    { \
        auto Integer d ; \
        for ( d = outsideGrid = 0 ; d < self->ndimensions ; d++ ) \
        { \
            RegularGridDimension_RegularizeIndexToMinusOne ( self->dimensions[d], self->workI[d] ) ; \
            if ( self->workI[d] < 0 ) outsideGrid += 1 ; \
        } \
    }

/* . Regularize internally stored indices with indices outside the grid forced to the boundary. */
# define RegularGrid_RegularizeIndicesToBoundary( self, outsideGrid ) \
    { \
        auto Boolean isOK ; \
        auto Integer d    ; \
        for ( d = outsideGrid = 0 ; d < self->ndimensions ; d++ ) \
        { \
            RegularGridDimension_RegularizeIndexToBoundary ( self->dimensions[d], self->workI[d], isOK ) ; \
            if ( ! isOK ) outsideGrid += 1 ; \
        } \
    }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Regular grid dimension.*/
extern void                    RegularGridDimension_CopyTo                        ( const RegularGridDimension    *self                   ,
                                                                                          RegularGridDimension    *other                  ) ;
extern void                    RegularGridDimension_Initialize                    (       RegularGridDimension    *self                   ) ;

/* . Regular grid search range. */
extern RegularGridSearchRange *RegularGridSearchRange_Allocate                    ( const Integer                  cells                  ,
                                                                                    const Integer                  dimensions             ,
                                                                                          Status                  *status                 ) ;
extern void                    RegularGridSearchRange_Deallocate                  (       RegularGridSearchRange **self                   ) ;
extern void                    RegularGridSearchRange_Print                       ( const RegularGridSearchRange  *self                   ) ;

/* . Regular grid. */
extern RegularGrid            *RegularGrid_Allocate                               ( const Integer                  ndimensions            ,
                                                                                          Status                  *status                 ) ;
extern void                    RegularGrid_CellIDToIndices                        ( const RegularGrid             *self                   ,
                                                                                    const Integer                  cellID                 ,
                                                                                          Integer                 *indices                ,
                                                                                          Status                  *status                 ) ;
extern Integer                 RegularGrid_CellIndicesToID                        ( const RegularGrid             *self                   ,
                                                                                    const Integer                 *indices                ) ;
extern RegularGrid            *RegularGrid_Clone                                  ( const RegularGrid             *self                   ,
                                                                                          Status                  *status                 ) ;
extern void                    RegularGrid_Deallocate                             (       RegularGrid            **self                   ) ;
extern void                    RegularGrid_DeallocateVoid                         (       void                   **vSelf                  ) ;
extern void                    RegularGrid_DistancesSquaredBetweenCellAndPoint    ( const RegularGrid             *self                   ,
                                                                                    const Integer                 *indices                ,
                                                                                    const Real                    *coordinates            ,
                                                                                          Real                    *maximumDistanceSquared ,
                                                                                          Real                    *minimumDistanceSquared ,
                                                                                          Status                  *status                 ) ;
extern void                    RegularGrid_DistancesSquaredBetweenCells           ( const RegularGrid             *self                   ,
                                                                                    const Integer                 *indices1               ,
                                                                                    const Integer                 *indices2               ,
                                                                                          Real                    *maximumDistanceSquared ,
                                                                                          Real                    *minimumDistanceSquared ,
                                                                                          Status                  *status                 ) ;
extern Integer                 RegularGrid_FindCellIDOfPoint                      ( const RegularGrid             *self                   ,
                                                                                    const Real                    *coordinates            ) ;
extern Integer                 RegularGrid_FindCellIndicesOfPoint                 ( const RegularGrid             *self                   ,
                                                                                    const Real                    *coordinates            ,
                                                                                    const Boolean                  regularizeToBoundary   ,
                                                                                          Integer                 *indices                ,
                                                                                          Real                    *fractional             ) ;
extern Integer                 RegularGrid_FindCellsWithinRangeOfCell             ( const RegularGrid             *self                   ,
                                                                                    const Integer                  cellID                 ,
                                                                                          RegularGridSearchRange  *range                  ,
                                                                                          Status                  *status                 ) ;
extern Integer                 RegularGrid_FindCellsWithinRangeOfPoint            ( const RegularGrid             *self                   ,
                                                                                    const Integer                  cellID                 ,
                                                                                    const Real                    *coordinates            ,
                                                                                          RegularGridSearchRange  *range                  ,
                                                                                          Status                  *status                 ) ;
extern Integer                 RegularGrid_FindConformingCellsWithinRangeOfCell   ( const RegularGrid             *self                   ,
                                                                                    const Integer                  cellID                 ,
                                                                                    const RegularGrid             *other                  ,
                                                                                    const IntegerArray1D          *offset                 ,
                                                                                          RegularGridSearchRange  *range                  ,
                                                                                          Status                  *status                 ) ;
extern Integer                 RegularGrid_FindConformingCellsWithinRangeOfPoint  ( const RegularGrid             *self                   ,
                                                                                    const Integer                  cellID                 ,
                                                                                    const Real                    *coordinates            ,
                                                                                    const RegularGrid             *other                  ,
                                                                                    const IntegerArray1D          *offset                 ,
                                                                                          RegularGridSearchRange  *range                  ,
                                                                                          Status                  *status                 ) ;
extern void                    RegularGrid_GetGridPointCoordinates                ( const RegularGrid             *self                   ,
                                                                                    const Integer                 *indices                ,
                                                                                          Real                    *coordinates            ) ;
extern Boolean                 RegularGrid_IsConformingRealArrayND                ( const RegularGrid             *self                   ,
                                                                                    const RealArrayND             *data                   ) ;
extern Boolean                 RegularGrid_IsCubic                                ( const RegularGrid             *self                   ) ;
extern Boolean                 RegularGrid_IsPeriodic                             ( const RegularGrid             *self                   ) ;
extern Integer                 RegularGrid_NumberOfGridPoints                     ( const RegularGrid             *self                   ) ;
extern RegularGrid            *RegularGrid_MakePeriodicGrid3                      ( const Vector3                 *boxSize                ,
                                                                                    const Real                     approximateGridSize    ,
                                                                                          Status                  *status                 ) ;
extern RegularGridSearchRange *RegularGrid_MakeSearchRange                        ( const RegularGrid             *self                   ,
                                                                                    const Real                     cutOff                 ,
                                                                                          Status                  *status                 ) ;
extern void                    RegularGrid_Print                                  ( const RegularGrid             *self                   ) ;

# endif

