/*==================================================================================================================================
! . Procedures for pairlist generation from coordinates3.
! . Both direct and grid-based versions are provided.
!=================================================================================================================================*/

/*
! . The cell/cell generation method produces lists that are not sorted with respect to i.
! . Sorting with respect to j is optional although it can contribute substantially to the time.
*/

# include <math.h>
# include <stdio.h>

# include "BooleanBlock.h"
# include "BooleanUtilities.h"
# include "IntegerUtilities.h"
# include "Memory.h"
# include "PairListGenerator.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Self-pairlist generation macros.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Cross interaction - direct with radii. */
# define CheckForDirectCrossInteractionWithRadii \
    if ( isIncluded[j] && ( QORI || QOR2[j] ) ) \
    { \
        rij = ri + Array1D_Item ( radii2, j ) ; \
        Coordinates3_GetRow ( coordinates32, j, xj, yj, zj ) ; \
        dx = xi - xj ; \
        dy = yi - yj ; \
        dz = zi - zj ; \
        if ( ( dx * dx + dy * dy + dz * dz ) <= ( rij * rij ) ) { indices[n] = j ; n++ ; } \
    }

/* . Cross interaction - direct with no radii. */
# define CheckForDirectCrossInteractionWithNoRadii \
    if ( isIncluded[j] && ( QORI || QOR2[j] ) ) \
    { \
        Coordinates3_GetRow ( coordinates32, j, xj, yj, zj ) ; \
        dx = xi - xj ; \
        dy = yi - yj ; \
        dz = zi - zj ; \
        if ( ( dx * dx + dy * dy + dz * dz ) <= cutOffSquared ) { indices[n] = j ; n++ ; } \
    }

/* . Cross interaction - grid with radii. */
# define CheckForGridCrossInteractionWithRadii( p, pStart ) \
    { \
        j = Array1D_Item ( occupancy2->cellPoints, p+pStart ) ; \
        if ( isIncluded[j] && ( QORI || QOR2[j] ) ) \
        { \
	    rij = ri + Array1D_Item ( radii2, j ) ; \
            Coordinates3_GetRow ( coordinates32, j, xj, yj, zj ) ; \
            dx = xi - xj ; \
	    dy = yi - yj ; \
	    dz = zi - zj ; \
	    if ( ( dx * dx + dy * dy + dz * dz ) <= ( rij * rij ) ) { indices[n] = j ; n++ ; } \
        } \
    }

/* . Cross interaction - grid with no radii. */
# define CheckForGridCrossInteractionWithNoRadii( p, pStart ) \
    { \
        j = Array1D_Item ( occupancy2->cellPoints, p+pStart ) ; \
        if ( isIncluded[j] && ( QORI || QOR2[j] ) ) \
        { \
            Coordinates3_GetRow ( coordinates32, j, xj, yj, zj ) ; \
            dx = xi - xj ; \
	    dy = yi - yj ; \
	    dz = zi - zj ; \
	    if ( ( dx * dx + dy * dy + dz * dz ) <= cutOffSquared ) { indices[n] = j ; n++ ; } \
	} \
    }

/* . Cross interaction - grid with all interactions. */
# define IncludeAllCellCrossInteractions( p, pStart ) \
    { \
        j = Array1D_Item ( occupancy2->cellPoints, p+pStart ) ; \
        if ( isIncluded[j] && ( QORI || QOR2[j] ) ) { indices[n] = j ; n++ ; } \
    }

/* . Get information about a cell. */
# define GetCellInformation( m, occupancy, c, includeAll, pStart, pTotal ) \
    { \
        c          = Array1D_Item ( gridSearchRange->cellIDs            , m ) ; \
        includeAll = Array1D_Item ( gridSearchRange->isFullyWithinRange , m ) ; \
        pStart     = Array1D_Item ( occupancy->cellFirstPoints          , c ) ; \
        pTotal     = Array1D_Item ( occupancy->cellTotalPoints          , c ) ; \
    }

/* . Initialize isIncluded. */
# define InitializeIsIncluded( numberOfPoints, QAND ) \
    { \
        for ( i = 0 ; i < numberOfPoints ; i++ ) isIncluded[i] = QAND[i] ; \
    }

/* . Non-self exclusions. */
# define FlagNonSelfExclusions \
    { \
        xFirst = exclusions->itemsI[i]   ; \
        xLast  = exclusions->itemsI[i+1] ; \
        for ( x = xFirst ; x < xLast ; x++ ) { j = exclusions->itemsJ[x] ; isIncluded[j] = False ; } \
    }

# define UnflagNonSelfExclusions( QAND ) \
    { \
        for ( x = xFirst ; x < xLast ; x++ ) { j = exclusions->itemsJ[x] ; if ( QAND[j] ) isIncluded[j] = True ; } \
    }

/* . Non-self interaction - grid with radii. */
# define CheckForGridNonSelfInteractionWithRadii( p, pStart ) \
    { \
        j = Array1D_Item ( occupancy->cellPoints, p+pStart ) ; \
        if ( isIncluded[j] && ( QORI || QOR[j] ) ) \
        { \
	    rij = ri + Array1D_Item ( radii, j ) ; \
            Coordinates3_DifferenceRow ( coordinates3, i, j, xij, yij, zij ) ; \
            if ( ( ( xij * xij ) + ( yij * yij ) + ( zij * zij ) ) <= ( rij * rij ) ) { indices[n] = j ; n++ ; } \
        }  \
    }

/* . Non-self interaction - grid with no radii. */
# define CheckForGridNonSelfInteractionWithNoRadii( p, pStart ) \
    { \
        j = Array1D_Item ( occupancy->cellPoints, p+pStart ) ; \
        if ( isIncluded[j] && ( QORI || QOR[j] ) ) \
        { \
            Coordinates3_DifferenceRow ( coordinates3, i, j, xij, yij, zij ) ; \
            if ( ( ( xij * xij ) + ( yij * yij ) + ( zij * zij ) ) <= cutOffSquared ) { indices[n] = j ; n++ ; } \
        } \
    }

/* . Non-self interaction - grid with all interactions. */
# define IncludeAllCellNonSelfInteractions( p, pStart ) \
    { \
        j = Array1D_Item ( occupancy->cellPoints, p+pStart ) ; \
        if ( isIncluded[j] && ( QORI || QOR[j] ) ) { indices[n] = j ; n++ ; } \
    }

/* . Self exclusions. */
# define FlagSelfExclusions \
    { \
        xFirst = exclusions->itemsI[i]   ; \
        xLast  = exclusions->itemsI[i+1] ; \
        for ( x = xFirst ; x < xLast ; x++ ) { j = exclusions->itemsJ[x] ; if ( j >= i ) break ; isIncluded[j] = False ; } \
    }

# define UnflagSelfExclusions \
    { \
        for ( x = xFirst ; x < xLast ; x++ ) { j = exclusions->itemsJ[x] ; if ( j >= i ) break ; if ( QAND[j] ) isIncluded[j] = True ; } \
    }

/* . Self interaction - direct with radii. */
# define CheckForDirectSelfInteractionWithRadii \
    if ( isIncluded[j] && ( QORI || QOR[j] ) ) \
    { \
        rij = ri + Array1D_Item ( radii, j ) ; \
        Coordinates3_DifferenceRow ( coordinates3, i, j, xij, yij, zij ) ; \
        if ( ( ( xij * xij ) + ( yij * yij ) + ( zij * zij ) ) <= ( rij * rij ) ) { indices[n] = j ; n++ ; } \
    }

/* . Self interaction - direct with no radii. */
# define CheckForDirectSelfInteractionWithNoRadii \
    if ( isIncluded[j] && ( QORI || QOR[j] ) ) \
    { \
        Coordinates3_DifferenceRow ( coordinates3, i, j, xij, yij, zij ) ; \
        if ( ( ( xij * xij ) + ( yij * yij ) + ( zij * zij ) ) <= cutOffSquared ) { indices[n] = j ; n++ ; } \
    }

/* . Self interaction - grid with radii. */
# define CheckForGridSelfInteractionWithRadii( p, pStart ) \
    { \
        j = Array1D_Item ( occupancy->cellPoints, p+pStart ) ; \
        if ( j >= i ) break ; \
        else if ( isIncluded[j] && ( QORI || QOR[j] ) ) \
        { \
	    rij = ri + Array1D_Item ( radii, j ) ; \
            Coordinates3_DifferenceRow ( coordinates3, i, j, xij, yij, zij ) ; \
            if ( ( ( xij * xij ) + ( yij * yij ) + ( zij * zij ) ) <= ( rij * rij ) ) { indices[n] = j ; n++ ; } \
        }  \
    }

/* . Self interaction - grid with no radii. */
# define CheckForGridSelfInteractionWithNoRadii( p, pStart ) \
    { \
        j = Array1D_Item ( occupancy->cellPoints, p+pStart ) ; \
        if ( j >= i ) break ; \
        else if ( isIncluded[j] && ( QORI || QOR[j] ) ) \
        { \
            Coordinates3_DifferenceRow ( coordinates3, i, j, xij, yij, zij ) ; \
            if ( ( ( xij * xij ) + ( yij * yij ) + ( zij * zij ) ) <= cutOffSquared ) { indices[n] = j ; n++ ; } \
        } \
    }

/* . Self interaction - grid with all interactions. */
# define IncludeAllCellSelfInteractions( p, pStart ) \
    { \
        j = Array1D_Item ( occupancy->cellPoints, p+pStart ) ; \
        if ( j >= i ) break ; \
        else if ( isIncluded[j] && ( QORI || QOR[j] ) ) { indices[n] = j ; n++ ; } \
    }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void    AllocateTemporarySpace ( const Integer    numberOfPoints ,
                                        Integer        **indices        ,
                                        Boolean        **isIncluded     ,
                                        Status          *status         ) ;

static Boolean CheckRadii             ( const Integer      numberOfPoints ,
                                        const RealArray1D *radii          ,
                                        Real              *maximumRadius  ,
                                        Real              *minimumRadius  ,
                                        Status            *status         ) ;

static void    CheckSelections        ( const Integer     numberOfPoints ,
                                        Selection        *andSelection   ,
                                        Selection        *orSelection    ,
                                        Boolean         **QAND           ,
                                        Boolean         **QOR            ,
                                        Status           *status         ) ;

static PairList *MakeCrossPairListFromCoordinates3 ( const PairListGenerator *self          ,
                                                     const Coordinates3      *coordinates31 ,
                                                     const Coordinates3      *coordinates32 ,
                                                     const RealArray1D       *radii1        ,
                                                     const RealArray1D       *radii2        ,
                                                     Selection               *andSelection1 ,
                                                     Selection               *andSelection2 ,
                                                     Selection               *orSelection1  ,
                                                     Selection               *orSelection2  ,
                                                     const PairConnections   *exclusions    ,
                                                     const Boolean            excludeSelf   ,
                                                     RegularGrid             *grid1         ,
                                                     RegularGrid             *grid2         ,
                                                     RegularGridOccupancy    *occupancy1    ,
                                                     RegularGridOccupancy    *occupancy2    ,
                                                     IntegerArray1D          *offSet        ,
                                                     Status                  *status        ) ;

static Boolean MakeCrossPairListFromCoordinates3CellCell ( PairList *pairList, const Boolean                 sortIndices     ,
                                                                               const Real                    cutOff          ,
                                                                               const Coordinates3           *coordinates31   ,      
                                                                               const Coordinates3           *coordinates32   ,
                                                                               const RealArray1D            *radii1          ,
                                                                               const RealArray1D            *radii2          ,
                                                                               const PairConnections        *exclusions      ,
                                                                               const Boolean                 excludeSelf     ,
                                                                               const Boolean                *QAND1           ,
                                                                               const Boolean                *QAND2           ,
                                                                               const Boolean                *QOR1            ,
                                                                               const Boolean                *QOR2            ,
                                                                               const RegularGrid            *grid1           ,
                                                                               const RegularGrid            *grid2           ,
                                                                               const RegularGridOccupancy   *occupancy1      ,
                                                                               const RegularGridOccupancy   *occupancy2      ,
                                                                               const IntegerArray1D         *offSet          ,
                                                                                     RegularGridSearchRange *gridSearchRange ,
                                                                                     Boolean                *isIncluded      ,
                                                                                     Integer                *indices         ) ;

static Boolean MakeCrossPairListFromCoordinates3Direct ( PairList *pairList, const Real             cutOff          ,
                                                                             const Coordinates3    *coordinates31   ,   
                                                                             const Coordinates3    *coordinates32   ,
                                                                             const RealArray1D     *radii1          ,
                                                                             const RealArray1D     *radii2          ,
                                                                             const PairConnections *exclusions      ,
                                                                             const Boolean          excludeSelf     ,
                                                                             const Boolean         *QAND1           ,
                                                                             const Boolean         *QAND2           ,
                                                                             const Boolean         *QOR1            ,
                                                                             const Boolean         *QOR2            ,
                                                                                   Boolean         *isIncluded      ,
                                                                                   Integer         *indices         ) ;

static Boolean MakeCrossPairListFromCoordinates3PointCell ( PairList *pairList, const Boolean                 sortIndices     ,
                                                                                const Real                    cutOff          ,
                                                                                const Coordinates3           *coordinates31   ,      
                                                                                const Coordinates3           *coordinates32   ,
                                                                                const RealArray1D            *radii1          ,
                                                                                const RealArray1D            *radii2          ,
                                                                                const PairConnections        *exclusions      ,
                                                                                const Boolean                 excludeSelf     ,
                                                                                const Boolean                *QAND1           ,
                                                                                const Boolean                *QAND2           ,
                                                                                const Boolean                *QOR1            ,
                                                                                const Boolean                *QOR2            ,
                                                                                const RegularGrid            *grid1           ,
                                                                                const RegularGrid            *grid2           ,
                                                                                const RegularGridOccupancy   *occupancy1      ,
                                                                                const RegularGridOccupancy   *occupancy2      ,
                                                                                const IntegerArray1D         *offSet          ,
                                                                                      RegularGridSearchRange *gridSearchRange ,
                                                                                      Boolean                *isIncluded      ,
                                                                                      Integer                *indices         ) ;

static PairList *MakeSelfPairListFromCoordinates3 ( const PairListGenerator    *self         ,
                                                    const Coordinates3         *coordinates3 ,
                                                    const RealArray1D          *radii        ,
                                                          Selection            *andSelection ,
                                                          Selection            *orSelection  ,
                                                    const PairConnections      *exclusions   ,
                                                          RegularGrid          *grid         ,
                                                          RegularGridOccupancy *occupancy    ,
                                                          Status               *status       ) ;

static Boolean MakeSelfPairListFromCoordinates3CellCell ( PairList *pairList, const Boolean                 sortIndices     ,
                                                                              const Real                    cutOff          ,
                                                                              const Coordinates3           *coordinates3    ,
                                                                              const RealArray1D            *radii           ,
                                                                              const PairConnections        *exclusions      ,
                                                                              const Boolean                *QAND            ,
                                                                              const Boolean                *QOR             ,
                                                                              const RegularGrid            *grid            ,
                                                                              const RegularGridOccupancy   *occupancy       ,
                                                                                    RegularGridSearchRange *gridSearchRange ,
                                                                                    Boolean                *isIncluded      ,
                                                                                    Integer                *indices         ) ;

static Boolean MakeSelfPairListFromCoordinates3Direct ( PairList *pairList, const Real             cutOff       ,
                                                                            const Coordinates3    *coordinates3 ,
                                                                            const RealArray1D     *radii        ,
                                                                            const PairConnections *exclusions   ,
                                                                            const Boolean         *QAND         ,
                                                                            const Boolean         *QOR          ,
                                                                                  Boolean         *isIncluded   ,
                                                                                  Integer         *indices      ) ;

static Boolean MakeSelfPairListFromCoordinates3PointCell ( PairList *pairList, const Boolean                 sortIndices     ,
                                                                               const Real                    cutOff          ,
                                                                               const Coordinates3           *coordinates3    ,
                                                                               const RealArray1D            *radii           ,
                                                                               const PairConnections        *exclusions      ,
                                                                               const Boolean                *QAND            ,
                                                                               const Boolean                *QOR             ,
                                                                               const RegularGrid            *grid            ,
                                                                               const RegularGridOccupancy   *occupancy       ,
                                                                                     RegularGridSearchRange *gridSearchRange ,
                                                                                     Boolean                *isIncluded      ,
                                                                                     Integer                *indices         ) ;

static Boolean SaveInteractions (       PairList *pairList    , 
                                  const Integer   i           , 
                                  const Integer   n           , 
                                  const Integer  *indices     , 
                                  const Boolean   sortIndices ) ;

/*==================================================================================================================================
! . Public procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairListGenerator *PairListGenerator_Allocate ( void )
{
    PairListGenerator *self = Memory_AllocateType ( PairListGenerator ) ;
    if ( self != NULL )
    {
        /* . Set defaults. */
        self->sortIndices          = True    ;
        self->useGridByCell        = False   ;
        self->minimumCellExtent    =   2     ;
        self->minimumPoints        = 500     ;
        self->cellSize             = 0.0e+00 ;
        self->cutOff               = 0.0e+00 ;
        self->cutOffCellSizeFactor = 0.5e+00 ;
        self->minimumCellSize      = 1.0e+00 ;
        self->minimumExtentFactor  = 1.5e+00 ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairListGenerator *PairListGenerator_Clone ( const PairListGenerator *self )
{
    PairListGenerator *new = NULL ;
    if ( self != NULL )
    {
        new = PairListGenerator_Allocate ( ) ;
        new->sortIndices          = self->sortIndices          ;
        new->useGridByCell        = self->useGridByCell        ;
        new->minimumCellExtent    = self->minimumCellExtent    ;
        new->minimumPoints        = self->minimumPoints        ;
        new->cellSize             = self->cellSize             ;
        new->cutOff               = self->cutOff               ;
        new->cutOffCellSizeFactor = self->cutOffCellSizeFactor ;
        new->minimumCellSize      = self->minimumCellSize      ;
        new->minimumExtentFactor  = self->minimumExtentFactor  ;
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairListGenerator_Deallocate ( PairListGenerator **self )
{
    if ( (*self) != NULL ) Memory_Deallocate ( (*self) ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Create a cross pairlist from two sets of coordinates.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Double case. */
PairList *PairListGenerator_CrossPairListFromDoubleCoordinates3 ( const PairListGenerator    *self          ,
                                                                  const Coordinates3         *coordinates31 ,
                                                                  const Coordinates3         *coordinates32 ,
                                                                  const RealArray1D          *radii1        ,
                                                                  const RealArray1D          *radii2        ,
                                                                        Selection            *andSelection1 ,
                                                                        Selection            *andSelection2 ,
                                                                        Selection            *orSelection1  ,
                                                                        Selection            *orSelection2  ,
                                                                        PairList             *exclusions    ,
                                                                        RegularGrid          *grid1         ,
                                                                        RegularGridOccupancy *occupancy1    ,
                                                                        Status               *status        )
{
    PairList *pairList = NULL ;
    if ( ( self != NULL ) &&
         ( ! Selection_IsEmpty ( andSelection1, True ) ) &&
         ( ! Selection_IsEmpty ( andSelection2, True ) ) &&
         ( ! Selection_IsEmpty (  orSelection1, True ) ) &&
         ( ! Selection_IsEmpty (  orSelection2, True ) ) )
    {
        /* . The exclusion pairlist must be a cross-pairlist. */
        if ( ( exclusions == NULL ) || ( ( exclusions != NULL ) && ( ! exclusions->isSelf ) ) )
        {
            auto Integer               numberOfPoints = 0 ;
            auto IntegerArray1D       *offSet         = NULL ;
            auto PairConnections      *connections    = NULL ;
            auto RegularGrid          *grid2          = NULL ;
            auto RegularGridOccupancy *occupancy2     = NULL ;
            if ( coordinates31 != NULL ) numberOfPoints = coordinates31->extent0 ;
            connections = CrossPairList_MakeConnections ( exclusions, numberOfPoints, status ) ;
            if ( grid1 != NULL ) Coordinates3_MakeConformingGridAndOccupancy ( coordinates32 ,
                                                                               andSelection2 ,
                                                                               grid1         ,
                                                                               &grid2        ,
                                                                               &occupancy2   ,
                                                                               &offSet       ,
                                                                               status        ) ;
            pairList = MakeCrossPairListFromCoordinates3 ( self          ,
                                                           coordinates31 ,
                                                           coordinates32 ,
                                                           radii1        ,
                                                           radii2        ,
                                                           andSelection1 ,
                                                           andSelection2 ,
                                                           orSelection1  ,
                                                           orSelection2  ,
                                                           connections   ,
                                                           False         ,
                                                           grid1         ,
                                                           grid2         ,
                                                           occupancy1    ,
                                                           occupancy2    ,
                                                           offSet        ,
                                                           status        ) ;
            if ( grid1 != NULL )
            {
                IntegerArray1D_Deallocate       ( &offSet     ) ;
                RegularGrid_Deallocate          ( &grid2      ) ;
                RegularGridOccupancy_Deallocate ( &occupancy2 ) ;
            }
        }
        else Status_Set ( status, Status_InvalidArgument ) ;
    }
    return pairList ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Create a cross pairlist from one set of coordinates.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The second grid is not constructed if the first grid and its occupancy already contain all the points. */
PairList *PairListGenerator_CrossPairListFromSingleCoordinates3 ( const PairListGenerator    *self          ,
                                                                  const Coordinates3         *coordinates3  ,
                                                                  const RealArray1D          *radii         ,
                                                                        Selection            *andSelection1 ,
                                                                        Selection            *andSelection2 ,
                                                                        Selection            *orSelection   ,
                                                                        PairList             *exclusions    ,
                                                                  const Boolean               excludeSelf   ,
                                                                        RegularGrid          *grid1         ,
                                                                        RegularGridOccupancy *occupancy1    ,
                                                                        Status               *status        )
{
    PairList *pairList = NULL ;
    if ( ( self != NULL ) &&
         ( ! Selection_IsEmpty ( andSelection1, True ) ) &&
         ( ! Selection_IsEmpty ( andSelection2, True ) ) &&
         ( ! Selection_IsEmpty (  orSelection , True ) ) )
    {
        /* . The exclusion pairlist must be a self-pairlist. */
        if ( ( exclusions == NULL ) || ( ( exclusions != NULL ) && ( exclusions->isSelf ) ) )
        {
            auto Boolean               madeGrid2      = False ;
            auto Integer               numberOfPoints = 0 ;
            auto IntegerArray1D       *offSet         = NULL ;
            auto PairConnections      *connections    = NULL ;
            auto RegularGrid          *grid2          = NULL ;
            auto RegularGridOccupancy *occupancy2     = NULL ;
            if ( coordinates3 != NULL ) numberOfPoints = coordinates3->extent0 ;
            connections = SelfPairList_MakeConnections ( exclusions, numberOfPoints, status ) ;
            if ( grid1 != NULL )
            {
                madeGrid2 = ( Coordinates3_Rows ( coordinates3 ) != occupancy1->numberOfPoints ) ;
                if ( madeGrid2 ) Coordinates3_MakeConformingGridAndOccupancy ( coordinates3  ,
                                                                               andSelection2 ,
                                                                               grid1         ,
                                                                               &grid2        ,
                                                                               &occupancy2   ,
                                                                               &offSet       ,
                                                                               status        ) ;
                else
                {
                    grid2      = grid1      ;
                    occupancy2 = occupancy1 ;
                    offSet     = IntegerArray1D_AllocateWithExtent ( 3, status ) ;
                    IntegerArray1D_Set ( offSet, 0 ) ;
                }
            }
            pairList = MakeCrossPairListFromCoordinates3 ( self          ,
                                                           coordinates3  ,
                                                           coordinates3  ,
                                                           radii         ,
                                                           radii         ,
                                                           andSelection1 ,
                                                           andSelection2 ,
                                                           orSelection   ,
                                                           orSelection   ,
                                                           connections   ,
                                                           excludeSelf   ,
                                                           grid1         ,
                                                           grid2         ,
                                                           occupancy1    ,
                                                           occupancy2    ,
                                                           offSet        ,
                                                           status        ) ;
            if ( grid1 != NULL ) IntegerArray1D_Deallocate ( &offSet ) ;
            if ( madeGrid2 )
            {
                RegularGrid_Deallocate          ( &grid2      ) ;
                RegularGridOccupancy_Deallocate ( &occupancy2 ) ;
            }
        }
        else Status_Set ( status, Status_InvalidArgument ) ;
    }
    return pairList ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Decide on whether to use a direct or a grid-based search.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean PairListGenerator_DetermineMethod ( const PairListGenerator *self         ,
                                            const Coordinates3      *coordinates3 ,
                                            const Selection         *andSelection )
{    
    Boolean useGridSearch = False ;
    if ( ( self         != NULL ) &&
         ( coordinates3 != NULL ) &&
         ( ! Selection_IsEmpty ( andSelection, True ) ) )
    {
        auto Integer  maximumCells, numberOfPoints ;
        auto Real     maximumExtent ;
        auto Vector3 *extents = NULL, *origin = NULL ;

        /* . Allocation. */
        extents = Vector3_Allocate ( ) ;
        origin  = Vector3_Allocate ( ) ;

        /* . Find the extents of the coordinates. */
        Coordinates3_EnclosingOrthorhombicBox ( coordinates3, andSelection, NULL, origin, extents ) ;

        /* . Determine whether to use a direct or a grid-based search. */
        maximumExtent  = Vector3_AbsoluteMaximum ( extents ) ;
        maximumCells   = ( Integer  ) ceil ( maximumExtent / self->cellSize ) ;
        numberOfPoints = Coordinates3_Rows ( coordinates3 ) ; 
        useGridSearch  = ( ( maximumExtent  > self->cutOffCellSizeFactor * self->cutOff ) &&
                           ( maximumCells   > self->minimumCellExtent                   ) &&
                           ( numberOfPoints > self->minimumPoints                       ) ) ;

        /* . Finish up. */
        Vector3_Deallocate ( &extents ) ;
        Vector3_Deallocate ( &origin  ) ;
    }
    return useGridSearch ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Create a self pairlist from a set of coordinates.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairList *PairListGenerator_SelfPairListFromCoordinates3 ( const PairListGenerator    *self         ,
                                                           const Coordinates3         *coordinates3 ,
                                                           const RealArray1D          *radii        ,
                                                                 Selection            *andSelection ,
                                                                 Selection            *orSelection  ,
                                                                 PairList             *exclusions   ,
                                                                 RegularGrid          *grid         ,
                                                                 RegularGridOccupancy *occupancy    ,
                                                                 Status               *status       )
{
    PairList *pairList = NULL ;
    if ( ( self != NULL ) &&
         ( ! Selection_IsEmpty ( andSelection, True ) ) &&
         ( ! Selection_IsEmpty (  orSelection, True ) ) )
    {
        /* . The exclusion pairlist must be a self pairlist. */
        if ( ( exclusions == NULL ) || ( ( exclusions != NULL ) && ( exclusions->isSelf ) ) )
        {
            auto Integer          numberOfPoints = 0 ;
            auto PairConnections *connections    = NULL ;
            if ( coordinates3 != NULL ) numberOfPoints = coordinates3->extent0 ;
            connections = SelfPairList_MakeConnections ( exclusions, numberOfPoints, status ) ;
            pairList    = MakeSelfPairListFromCoordinates3 ( self         ,
                                                             coordinates3 ,
                                                             radii        ,
                                                             andSelection ,
                                                             orSelection  ,
                                                             connections  ,
                                                             grid         ,
                                                             occupancy    ,
                                                             status       ) ;
        }
        else Status_Set ( status, Status_InvalidArgument ) ;
    }
    return pairList ;
}

/*==================================================================================================================================
! . Local procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocate some temporary space.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void AllocateTemporarySpace ( const Integer    numberOfPoints ,
                                           Integer  **indices        ,
                                           Boolean  **isIncluded     ,
                                           Status    *status         )
{
    (*indices   ) = Integer_Allocate ( numberOfPoints, NULL ) ;
    (*isIncluded) = Boolean_Allocate  ( numberOfPoints, NULL ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check a set of input radii.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Boolean CheckRadii ( const Integer      numberOfPoints ,
                            const RealArray1D *radii          ,
                                  Real        *maximumRadius  ,
                                  Real        *minimumRadius  ,
                                  Status      *status         )
{
    auto Boolean hasRadii ;
    hasRadii = ( radii != NULL ) ;
    if ( hasRadii )
    {
        if ( radii->extent != numberOfPoints ) Status_Set ( status, Status_InvalidArgument ) ;
        (*maximumRadius) = RealArray1D_Maximum ( radii ) ;
        (*minimumRadius) = RealArray1D_Minimum ( radii ) ;
    }
    else { (*maximumRadius) = 0.0e+00 ; (*minimumRadius) = 0.0e+00 ; }
    return hasRadii ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check a set of input selections.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CheckSelections ( const Integer     numberOfPoints ,
                                    Selection  *andSelection   ,
                                    Selection  *orSelection    ,
                                    Boolean   **QAND           ,
                                    Boolean   **QOR            ,
                                    Status     *status         )
{
    /* . Initialization. */
    (*QAND) = NULL ;
    (*QOR ) = NULL ;
    /* . Check the AND selection. */
    if ( andSelection == NULL )
    {
        (*QAND) = Boolean_Allocate ( numberOfPoints, NULL ) ;
        Boolean_Set ( (*QAND), numberOfPoints, True ) ;
    }
    else
    {
        auto BooleanBlock *flags = Selection_MakeFlags ( andSelection, numberOfPoints, NULL ) ;
        if ( flags != NULL ) (*QAND) = Block_Items ( flags ) ;
    }
    /* . Check the OR selection. */
    if ( orSelection == NULL )
    {
        (*QOR) = Boolean_Allocate ( numberOfPoints, NULL ) ;
        Boolean_Set ( (*QOR), numberOfPoints, True ) ;
    }
    else
    {
        auto BooleanBlock *flags = Selection_MakeFlags ( orSelection, numberOfPoints, NULL ) ;
        if ( flags != NULL ) (*QOR) = Block_Items ( flags ) ;
    }
    /* . Finish up. */
    if ( ( (*QAND) == NULL ) || ( (*QOR) == NULL ) ) Status_Set ( status, Status_OutOfMemory ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Create a cross-pairlist from two sets of coordinates3.
!---------------------------------------------------------------------------------------------------------------------------------*/
static PairList *MakeCrossPairListFromCoordinates3 ( const PairListGenerator    *self          ,
                                                     const Coordinates3         *coordinates31 ,          
                                                     const Coordinates3         *coordinates32 ,
                                                     const RealArray1D          *radii1        ,
                                                     const RealArray1D          *radii2        ,
                                                           Selection            *andSelection1 ,
                                                           Selection            *andSelection2 ,
                                                           Selection            *orSelection1  ,
                                                           Selection            *orSelection2  ,
                                                     const PairConnections      *exclusions    ,
                                                     const Boolean               excludeSelf   ,
                                                           RegularGrid          *grid1         ,
                                                           RegularGrid          *grid2         ,
                                                           RegularGridOccupancy *occupancy1    ,
                                                           RegularGridOccupancy *occupancy2    ,
                                                           IntegerArray1D       *offSet        ,
                                                           Status               *status        )
{
    PairList *pairList = NULL ;
    if ( ( self != NULL ) && ( coordinates31 != NULL ) && ( coordinates32 != NULL ) )
    {
        auto Boolean   emptyList, isOK, useGridSearch ;
        auto Integer   numberOfPoints1, numberOfPoints2 ;
        auto Real      maximumCutOff, maximumRadius1, maximumRadius2, minimumRadius1, minimumRadius2 ;
        auto Status    localStatus ;
        auto Boolean  *QAND1 = NULL, *QAND2 = NULL, *isIncluded = NULL, *QOR1 = NULL, *QOR2 = NULL ;
        auto Integer  *indices = NULL ;
        auto RegularGridSearchRange *gridSearchRange = NULL ;

        /* . Initialization. */
        localStatus = Status_OK ;

        /* . Get the dimensions of the problem. */
        numberOfPoints1 = coordinates31->extent0 ;
        numberOfPoints2 = coordinates32->extent0 ;

        /* . Check the exclusions, radii and selections. */
        CheckRadii ( numberOfPoints1, radii1, &maximumRadius1, &minimumRadius1, &localStatus ) ;
        CheckRadii ( numberOfPoints2, radii2, &maximumRadius2, &minimumRadius2, &localStatus ) ;
        CheckSelections ( numberOfPoints1, andSelection1, orSelection1, &QAND1, &QOR1, &localStatus ) ;
        CheckSelections ( numberOfPoints2, andSelection2, orSelection2, &QAND2, &QOR2, &localStatus ) ;

        /* . Check for an empty list. */
        maximumCutOff = self->cutOff + maximumRadius1 + maximumRadius2 ;
        emptyList     = ( maximumCutOff < 0.0e+00 ) ;

        /* . Set up for a grid search. */
        /* . Possible to use something else here appropriate for conforming grids? */
        useGridSearch = ( grid1 != NULL ) && ( grid2 != NULL ) && ( occupancy1 != NULL ) && ( occupancy2 != NULL ) && ( offSet != NULL ) ;
        if ( useGridSearch ) gridSearchRange = RegularGrid_MakeSearchRange ( grid1, maximumCutOff, &localStatus ) ;

        /* . Create the pairlist and allocate some temporary space. */
        pairList = PairList_Allocate ( numberOfPoints1, &localStatus ) ;
        AllocateTemporarySpace ( numberOfPoints2, &indices, &isIncluded, &localStatus ) ;

        /* . Check for a memory error. */
        isOK = ( pairList != NULL ) && ( localStatus == Status_OK ) ;
        if ( isOK && ( ! emptyList ) )
        {
            /* . Generate the list. */
            if ( useGridSearch )
            {
                if ( self->useGridByCell ) isOK = MakeCrossPairListFromCoordinates3CellCell  ( pairList          ,
                                                                                               self->sortIndices ,
                                                                                               self->cutOff      ,
                                                                                               coordinates31     ,
                                                                                               coordinates32     ,
                                                                                               radii1            ,
                                                                                               radii2            ,
                                                                                               exclusions        ,
                                                                                               excludeSelf       ,
                                                                                               QAND1             ,
                                                                                               QAND2             ,
                                                                                               QOR1              ,
                                                                                               QOR2              ,
                                                                                               grid1             ,
                                                                                               grid2             ,
                                                                                               occupancy1        ,
                                                                                               occupancy2        ,
                                                                                               offSet            ,
                                                                                               gridSearchRange   ,
                                                                                               isIncluded        ,
                                                                                               indices           ) ;
                else                       isOK = MakeCrossPairListFromCoordinates3PointCell ( pairList          ,
                                                                                               self->sortIndices ,
                                                                                               self->cutOff      ,
                                                                                               coordinates31     ,
                                                                                               coordinates32     ,
                                                                                               radii1            ,
                                                                                               radii2            ,
                                                                                               exclusions        ,
                                                                                               excludeSelf       ,
                                                                                               QAND1             ,
                                                                                               QAND2             ,
                                                                                               QOR1              ,
                                                                                               QOR2              ,
                                                                                               grid1             ,
                                                                                               grid2             ,
                                                                                               occupancy1        ,
                                                                                               occupancy2        ,
                                                                                               offSet            ,
                                                                                               gridSearchRange   ,
                                                                                               isIncluded        ,
                                                                                               indices           ) ;
            }
            else                           isOK = MakeCrossPairListFromCoordinates3Direct    ( pairList          ,
                                                                                               self->cutOff      ,
                                                                                               coordinates31     ,
                                                                                               coordinates32     ,
                                                                                               radii1            ,
                                                                                               radii2            ,
                                                                                               exclusions        ,
                                                                                               excludeSelf       ,
                                                                                               QAND1             ,
                                                                                               QAND2             ,
                                                                                               QOR1              ,
                                                                                               QOR2              ,
                                                                                               isIncluded        ,
                                                                                               indices           ) ;
        }

        /* . Finish up. */
        Memory_Deallocate ( isIncluded ) ;
        Memory_Deallocate ( indices    ) ;
        if ( andSelection1 == NULL ) Memory_Deallocate ( QAND1 ) ;
        if ( andSelection2 == NULL ) Memory_Deallocate ( QAND2 ) ;
        if ( orSelection1  == NULL ) Memory_Deallocate ( QOR1  ) ;
        if ( orSelection2  == NULL ) Memory_Deallocate ( QOR2  ) ;
        RegularGridSearchRange_Deallocate ( &gridSearchRange ) ;
        if ( ! isOK ) PairList_Deallocate ( &pairList ) ;
        Status_Set ( status, localStatus ) ;
    }
    return pairList ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cross-pairlist generation using a cell-cell method.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Boolean MakeCrossPairListFromCoordinates3CellCell ( PairList *pairList, const Boolean                 sortIndices     ,
                                                                               const Real                    cutOff          ,
                                                                               const Coordinates3           *coordinates31   ,      
                                                                               const Coordinates3           *coordinates32   ,
                                                                               const RealArray1D            *radii1          ,
                                                                               const RealArray1D            *radii2          ,
                                                                               const PairConnections        *exclusions      ,
                                                                               const Boolean                 excludeSelf     ,
                                                                               const Boolean                *QAND1           ,
                                                                               const Boolean                *QAND2           ,
                                                                               const Boolean                *QOR1            ,
                                                                               const Boolean                *QOR2            ,
                                                                               const RegularGrid            *grid1           ,
                                                                               const RegularGrid            *grid2           ,
                                                                               const RegularGridOccupancy   *occupancy1      ,
                                                                               const RegularGridOccupancy   *occupancy2      ,
                                                                               const IntegerArray1D         *offSet          ,
                                                                                     RegularGridSearchRange *gridSearchRange ,
                                                                                     Boolean                *isIncluded      ,
                                                                                     Integer                *indices         )
{
    Boolean  hasExclusions, hasRadii1, hasRadii2, includeAll, isOK = True, QORI ;
    Integer  c1, c2, i, j, m2, n, numberOfCells, numberOfPoints2, p1, p1Start, p1Total, p2, p2Start, p2Total, x, xFirst = 0, xLast = 0 ;
    Real     cutOffSquared, dx, dy, dz, ri, rij, xi, xj, yi, yj, zi, zj ;

    /* . Initialization. */
    cutOffSquared   = cutOff * cutOff ;
    hasExclusions   = ( exclusions != NULL ) ;
    hasRadii1       = ( radii1     != NULL ) ;
    hasRadii2       = ( radii2     != NULL ) ;
    numberOfPoints2 = coordinates32->extent0 ;
    InitializeIsIncluded ( numberOfPoints2, QAND2 ) ;

    /* . Outer loop over occupied cells in grid. */
    for ( c1 = 0 ; c1 < RegularGrid_NumberOfGridPoints ( grid1 ); c1++ )
    {
        /* . Set some information for the cell. */
        p1Start = Array1D_Item ( occupancy1->cellFirstPoints, c1 ) ;
        p1Total = Array1D_Item ( occupancy1->cellTotalPoints, c1 ) ;

        /* . Find all conforming cells that are within the cutOff of the current cell. */
        numberOfCells = RegularGrid_FindConformingCellsWithinRangeOfCell ( grid1, c1, grid2, offSet, gridSearchRange, NULL ) ;

        /* . Outer loop over indices. */
        for ( p1 = 0 ; p1 < p1Total ; p1++ )
        {
            /* . Get the point index. */
            i = Array1D_Item ( occupancy1->cellPoints, p1+p1Start ) ;

            /* . AND test. */
            if ( QAND1[i] )
            {
                /* . Get some information for the point. */
                Coordinates3_GetRow ( coordinates31, i, xi, yi, zi ) ;
                QORI = QOR1[i] ;

                /* . Flag excluded points for i. */
                if ( hasExclusions ) FlagNonSelfExclusions ;
                if ( excludeSelf && QAND2[i] ) isIncluded[i] = False ;

                /* . Set the cutOff. */
                ri = cutOff ;
                if ( hasRadii1 ) ri += Array1D_Item ( radii1, i ) ;

                /* . Loop over cells. */
                for ( m2 = n = 0 ; m2 < numberOfCells ; m2++ )
                {
                    /* . Set some information for the cell. */
                    GetCellInformation ( m2, occupancy2, c2, includeAll, p2Start, p2Total ) ;

                    /* . Generate interactions. */
                         if ( hasRadii2  ) { for ( p2 = 0 ; p2 < p2Total ; p2++ ) { CheckForGridCrossInteractionWithRadii   ( p2, p2Start ) ; } }
                    else if ( includeAll ) { for ( p2 = 0 ; p2 < p2Total ; p2++ ) { IncludeAllCellCrossInteractions         ( p2, p2Start ) ; } }
                    else                   { cutOffSquared = ri * ri ;
                                             for ( p2 = 0 ; p2 < p2Total ; p2++ ) { CheckForGridCrossInteractionWithNoRadii ( p2, p2Start ) ; } }
                }

                /* . Save the interactions. */
                isOK = SaveInteractions ( pairList, i, n, indices, sortIndices ) ;
                if ( ! isOK ) break ;

                /* . Unflag excluded points for i. */
                if ( hasExclusions ) UnflagNonSelfExclusions ( QAND2 ) ;
                if ( excludeSelf && QAND2[i] ) isIncluded[i] = True ;
            }
        }
    }
    return isOK ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cross-pairlist generation using a direct method.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Boolean MakeCrossPairListFromCoordinates3Direct ( PairList *pairList, const Real             cutOff        ,
                                                                             const Coordinates3    *coordinates31 ,   
                                                                             const Coordinates3    *coordinates32 ,
                                                                             const RealArray1D     *radii1        ,
                                                                             const RealArray1D     *radii2        ,
                                                                             const PairConnections *exclusions    ,
                                                                             const Boolean          excludeSelf   , 
                                                                             const Boolean         *QAND1         , 
                                                                             const Boolean         *QAND2         , 
                                                                             const Boolean         *QOR1          , 
                                                                             const Boolean         *QOR2          , 
                                                                                   Boolean         *isIncluded    , 
                                                                                   Integer         *indices       ) 
{
    Boolean  hasExclusions, hasRadii1, hasRadii2, isOK = True, QORI ;
    Integer  i, j, n, numberOfPoints1, numberOfPoints2, x, xFirst = 0, xLast = 0 ;
    Real     cutOffSquared, dx, dy, dz, ri, rij, xi, xj, yi, yj, zi, zj ;

    /* . Initialization. */
    cutOffSquared   = cutOff * cutOff ;
    hasExclusions   = ( exclusions != NULL ) ;
    hasRadii1       = ( radii1     != NULL ) ;
    hasRadii2       = ( radii2     != NULL ) ;
    numberOfPoints1 = coordinates31->extent0 ;
    numberOfPoints2 = coordinates32->extent0 ;
    InitializeIsIncluded ( numberOfPoints2, QAND2 ) ;

    /* . Outer loop over indices. */
    for ( i = 0 ; i < numberOfPoints1 ; i++ )
    {
        /* . AND test. */
        if ( QAND1[i] )
        {
            /* . Get some information for the point. */
            Coordinates3_GetRow ( coordinates31, i, xi, yi, zi ) ;
            QORI = QOR1[i] ;

            /* . Flag excluded points for i. */
            if ( hasExclusions ) FlagNonSelfExclusions ;
            if ( excludeSelf && QAND2[i] ) isIncluded[i] = False ;

            /* . Set the cutOff. */
            ri = cutOff ;
            if ( hasRadii1 ) ri += Array1D_Item ( radii1, i ) ;

            /* . Generate interactions. */
            if ( hasRadii2 ) { for ( j = n = 0 ; j < numberOfPoints2 ; j++ ) { CheckForDirectCrossInteractionWithRadii   ; } }
            else             { cutOffSquared = ri * ri ;
                               for ( j = n = 0 ; j < numberOfPoints2 ; j++ ) { CheckForDirectCrossInteractionWithNoRadii ; } }

            /* . Save the interactions. */
            isOK = SaveInteractions ( pairList, i, n, indices, False ) ;
            if ( ! isOK ) break ;

            /* . Unflag excluded points for i. */
            if ( hasExclusions ) UnflagNonSelfExclusions ( QAND2 ) ;
            if ( excludeSelf && QAND2[i] ) isIncluded[i] = True ;
        }
    }
    return isOK ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cross-pairlist generation using a point-cell method.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Boolean MakeCrossPairListFromCoordinates3PointCell ( PairList *pairList, const Boolean                 sortIndices     ,
                                                                                const Real                    cutOff          ,
                                                                                const Coordinates3           *coordinates31   ,      
                                                                                const Coordinates3           *coordinates32   ,
                                                                                const RealArray1D            *radii1          ,
                                                                                const RealArray1D            *radii2          ,
                                                                                const PairConnections        *exclusions      ,
                                                                                const Boolean                 excludeSelf     ,
                                                                                const Boolean                *QAND1           ,
                                                                                const Boolean                *QAND2           ,
                                                                                const Boolean                *QOR1            ,
                                                                                const Boolean                *QOR2            ,
                                                                                const RegularGrid            *grid1           ,
                                                                                const RegularGrid            *grid2           ,
                                                                                const RegularGridOccupancy   *occupancy1      ,
                                                                                const RegularGridOccupancy   *occupancy2      ,
                                                                                const IntegerArray1D         *offSet          ,
                                                                                      RegularGridSearchRange *gridSearchRange ,
                                                                                      Boolean                *isIncluded      ,
                                                                                      Integer                *indices         )
{
    Boolean  hasExclusions, hasRadii1, hasRadii2, includeAll, isOK = True, QORI ;
    Integer  c, i, j, m, n, numberOfCells, numberOfPoints1, numberOfPoints2, p, pStart, pTotal, x, xFirst = 0, xLast = 0 ;
    Real     cutOffSquared, dx, dy, dz, ri, rij, xi, xj, yi, yj, zi, zj ;

    /* . Initialization. */
    cutOffSquared   = cutOff * cutOff ;
    hasExclusions   = ( exclusions != NULL ) ;
    hasRadii1       = ( radii1     != NULL ) ;
    hasRadii2       = ( radii2     != NULL ) ;
    numberOfPoints1 = coordinates31->extent0 ;
    numberOfPoints2 = coordinates32->extent0 ;
    InitializeIsIncluded ( numberOfPoints2, QAND2 ) ;

    /* . Outer loop over indices. */
    for ( i = 0 ; i < numberOfPoints1 ; i++ )
    {
        /* . AND test. */
        if ( QAND1[i] )
        {
            /* . Get some information for the point. */
            Coordinates3_GetRow ( coordinates31, i, xi, yi, zi ) ;
            QORI = QOR1[i] ;

            /* . Flag excluded points for i. */
            if ( hasExclusions ) FlagNonSelfExclusions ;
            if ( excludeSelf && QAND2[i] ) isIncluded[i] = False ;

            /* . Set the cutOff. */
            ri = cutOff ;
            if ( hasRadii1 ) ri += Array1D_Item ( radii1, i ) ;

            /* . Find all cells that are within the cutOff of the current point. */
            numberOfCells = RegularGrid_FindConformingCellsWithinRangeOfPoint ( grid1, Array1D_Item ( occupancy1->pointCells, i ) ,
                                                                                            Coordinates3_RowPointer ( coordinates31, i ) ,
                                                                                                  grid2, offSet, gridSearchRange, NULL ) ;

            /* . Loop over cells. */
            for ( m = n = 0 ; m < numberOfCells ; m++ )
            {
                /* . Set some information for the cell. */
                GetCellInformation ( m, occupancy2, c, includeAll, pStart, pTotal ) ;

                /* . Generate interactions. */
                     if ( hasRadii2  ) { for ( p = 0 ; p < pTotal ; p++ ) { CheckForGridCrossInteractionWithRadii   ( p, pStart ) ; } }
                else if ( includeAll ) { for ( p = 0 ; p < pTotal ; p++ ) { IncludeAllCellCrossInteractions         ( p, pStart ) ; } }
                else                   { cutOffSquared = ri * ri ;
                                         for ( p = 0 ; p < pTotal ; p++ ) { CheckForGridCrossInteractionWithNoRadii ( p, pStart ) ; } }
            }

            /* . Save the interactions. */
            isOK = SaveInteractions ( pairList, i, n, indices, sortIndices ) ;
            if ( ! isOK ) break ;

            /* . Unflag excluded points for i. */
            if ( hasExclusions ) UnflagNonSelfExclusions ( QAND2 ) ;
            if ( excludeSelf && QAND2[i] ) isIncluded[i] = True ;
        }
    }
    return isOK ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Create a self-pairlist from a set of coordinates3.
!---------------------------------------------------------------------------------------------------------------------------------*/
static PairList *MakeSelfPairListFromCoordinates3 ( const PairListGenerator    *self         ,
                                                    const Coordinates3         *coordinates3 ,
                                                    const RealArray1D          *radii        ,
                                                          Selection            *andSelection ,
                                                          Selection            *orSelection  ,
                                                    const PairConnections      *exclusions   ,
                                                          RegularGrid          *grid         ,
                                                          RegularGridOccupancy *occupancy    ,
                                                          Status               *status       )
{
    PairList *pairList = NULL ;
    if ( ( self != NULL ) && ( coordinates3 != NULL ) && Status_IsOK ( status ) )
    {
        auto Boolean   emptyList, isOK, useGridSearch ;
        auto Integer   numberOfPoints ;
        auto Real      maximumCutOff, maximumRadius, minimumRadius ;
        auto Status    localStatus ;
        auto Boolean  *QAND = NULL, *isIncluded = NULL, *QOR = NULL ;
        auto Integer  *indices = NULL ;
        auto RegularGridSearchRange *gridSearchRange = NULL ;

        /* . Initialization. */
        localStatus = Status_OK ;

        /* . Get the dimension of the problem. */
        numberOfPoints = coordinates3->extent0 ;

        /* . Check the exclusions, radii and selections. */
        CheckRadii ( numberOfPoints, radii, &maximumRadius, &minimumRadius, &localStatus ) ;
        CheckSelections ( numberOfPoints, andSelection, orSelection, &QAND, &QOR, &localStatus ) ;

        /* . Check for an empty list. */
        maximumCutOff = self->cutOff + 2.0e+00 * maximumRadius ;
        emptyList     = ( maximumCutOff < 0.0e+00 ) ;

        /* . Set up for a grid search. */
        useGridSearch = ( grid != NULL ) && ( occupancy != NULL ) ;
        if ( useGridSearch ) gridSearchRange = RegularGrid_MakeSearchRange ( grid, maximumCutOff, &localStatus ) ;

        /* . Create the pairlist and allocate some temporary space. */
        pairList = PairList_Allocate ( numberOfPoints, &localStatus ) ;
        AllocateTemporarySpace ( numberOfPoints, &indices, &isIncluded, &localStatus ) ;

        /* . Check for a memory error. */
        isOK = ( pairList != NULL ) && ( localStatus == Status_OK ) ;
        if ( isOK && ( ! emptyList ) )
        {
            /* . Set the type of pair-list. */
            pairList->isSelf = True ;

            /* . Generate the list. */
            if ( useGridSearch )
            {
                if ( self->useGridByCell ) isOK = MakeSelfPairListFromCoordinates3CellCell  ( pairList          ,
                                                                                              self->sortIndices ,
                                                                                              self->cutOff      ,
                                                                                              coordinates3      ,
                                                                                              radii             ,
                                                                                              exclusions        ,
                                                                                              QAND              ,
                                                                                              QOR               ,
                                                                                              grid              ,
                                                                                              occupancy         ,
                                                                                              gridSearchRange   ,
                                                                                              isIncluded        ,
                                                                                              indices           ) ;
                else                       isOK = MakeSelfPairListFromCoordinates3PointCell ( pairList          ,
                                                                                              self->sortIndices ,
                                                                                              self->cutOff      ,
                                                                                              coordinates3      ,
                                                                                              radii             ,
                                                                                              exclusions        ,
                                                                                              QAND              ,
                                                                                              QOR               ,
                                                                                              grid              ,
                                                                                              occupancy         ,
                                                                                              gridSearchRange   ,
                                                                                              isIncluded        ,
                                                                                              indices           ) ;
            }
            else                           isOK = MakeSelfPairListFromCoordinates3Direct    ( pairList          ,
                                                                                              self->cutOff      ,
                                                                                              coordinates3      ,
                                                                                              radii             ,
                                                                                              exclusions        ,
                                                                                              QAND              ,
                                                                                              QOR               ,
                                                                                              isIncluded        ,
                                                                                              indices           ) ;
        }

        /* . Finish up. */
        Memory_Deallocate ( isIncluded ) ;
        Memory_Deallocate ( indices    ) ;
        if ( andSelection == NULL ) Memory_Deallocate ( QAND ) ;
        if ( orSelection  == NULL ) Memory_Deallocate ( QOR  ) ;
        RegularGridSearchRange_Deallocate ( &gridSearchRange ) ;
        if ( ! isOK ) PairList_Deallocate ( &pairList ) ;
        Status_Set ( status, localStatus ) ;
    }
    return pairList ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Self-pairlist generation using a grid cell-cell method.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Boolean MakeSelfPairListFromCoordinates3CellCell ( PairList *pairList, const Boolean                 sortIndices     ,
                                                                              const Real                    cutOff          ,
                                                                              const Coordinates3           *coordinates3    ,
                                                                              const RealArray1D            *radii           ,
                                                                              const PairConnections        *exclusions      ,
                                                                              const Boolean                *QAND            ,
                                                                              const Boolean                *QOR             ,
                                                                              const RegularGrid            *grid            ,
                                                                              const RegularGridOccupancy   *occupancy       ,
                                                                                    RegularGridSearchRange *gridSearchRange ,
                                                                                    Boolean                *isIncluded      ,
                                                                                    Integer                *indices         )
{
    Boolean  hasExclusions, hasRadii, includeAll, isOK = True, QORI ;
    Integer  c1, c2, i, j, m2, n, numberOfCells, numberOfPoints, p1, p1Start, p1Total, p2, p2Start, p2Total, p2Upper, x, xFirst = 0, xLast = 0 ;
    Real     cutOffSquared, ri = 0.0e+00, rij, xij, yij, zij ;

    /* . Initialization. */
    cutOffSquared  = cutOff * cutOff ;
    hasExclusions  = ( exclusions != NULL ) ;
    hasRadii       = ( radii      != NULL ) ;
    numberOfPoints = coordinates3->extent0  ;
    InitializeIsIncluded ( numberOfPoints, QAND ) ;

    /* . Outer loop over occupied cells in grid. */
    for ( c1 = 0 ; c1 < RegularGrid_NumberOfGridPoints ( grid ); c1++ )
    {
        /* . Set some information for the cell. */
        p1Start = Array1D_Item ( occupancy->cellFirstPoints, c1 ) ;
        p1Total = Array1D_Item ( occupancy->cellTotalPoints, c1 ) ;

        /* . Find all cells that are within the cutOff of the current cell. */
        numberOfCells = RegularGrid_FindCellsWithinRangeOfCell ( grid, c1, gridSearchRange, NULL ) ;

        /* . Outer loop over indices. */
        for ( p1 = 0 ; p1 < p1Total ; p1++ )
        {
            /* . Get the point index. */
            i = Array1D_Item ( occupancy->cellPoints, p1+p1Start ) ;

            /* . AND test. */
            if ( QAND[i] )
            {
                /* . Get some information for the point. */
                QORI = QOR[i] ;
                if ( hasRadii ) ri = cutOff + Array1D_Item ( radii, i ) ;

                /* . Flag excluded points for i. */
                if ( hasExclusions ) FlagNonSelfExclusions ;

                /* . Loop over cells. */
                for ( m2 = n = 0 ; m2 < numberOfCells ; m2++ )
                {
                    /* . Set some information for the cell. */
                    GetCellInformation ( m2, occupancy, c2, includeAll, p2Start, p2Total ) ;

                    /* . Skip interactions with a cell of higher index. */
                    if ( c2 > c1 ) continue ;

                    /* . Set the upper limit for interactions within the same cell. */
                    if ( c1 == c2 ) p2Upper = p1      ;
                    else            p2Upper = p2Total ;

                    /* . Generate interactions. */
                         if ( hasRadii   ) { for ( p2 = 0 ; p2 < p2Upper ; p2++ ) { CheckForGridNonSelfInteractionWithRadii   ( p2, p2Start ) ; } }
                    else if ( includeAll ) { for ( p2 = 0 ; p2 < p2Upper ; p2++ ) { IncludeAllCellNonSelfInteractions         ( p2, p2Start ) ; } }
                    else                   { for ( p2 = 0 ; p2 < p2Upper ; p2++ ) { CheckForGridNonSelfInteractionWithNoRadii ( p2, p2Start ) ; } }
                }

                /* . Save the interactions. */
                isOK = SaveInteractions ( pairList, i, n, indices, sortIndices ) ;
                if ( ! isOK ) break ;

                /* . Unflag excluded points for i. */
                if ( hasExclusions ) UnflagNonSelfExclusions ( QAND ) ;
            }
        }
    }
    return isOK ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Self-pairlist generation using a direct method.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Boolean MakeSelfPairListFromCoordinates3Direct ( PairList *pairList, const Real             cutOff       ,
                                                                            const Coordinates3    *coordinates3 ,
                                                                            const RealArray1D     *radii        ,
                                                                            const PairConnections *exclusions   ,
                                                                            const Boolean         *QAND         ,
                                                                            const Boolean         *QOR          ,
                                                                                  Boolean         *isIncluded   ,
                                                                                  Integer         *indices      )
{
    Boolean  hasExclusions, hasRadii, isOK = True, QORI ;
    Integer  i, j, n, numberOfPoints, x, xFirst = 0, xLast = 0 ;
    Real     cutOffSquared, ri = 0.0e+00, rij, xij, yij, zij ;

    /* . Initialization. */
    cutOffSquared  = cutOff * cutOff ;
    hasExclusions  = ( exclusions != NULL ) ;
    hasRadii       = ( radii      != NULL ) ;
    numberOfPoints = coordinates3->extent0  ;
    InitializeIsIncluded ( numberOfPoints, QAND ) ;

    /* . Outer loop over indices. */
    for ( i = 1 ; i < numberOfPoints ; i++ )
    {
        /* . AND test. */
        if ( QAND[i] )
        {
            /* . Get some information for the point. */
            QORI = QOR[i] ;
            if ( hasRadii ) ri = cutOff + Array1D_Item ( radii, i ) ;

            /* . Flag excluded points for i. */
            if ( hasExclusions ) FlagSelfExclusions ;

            /* . Generate interactions. */
            if ( hasRadii ) { for ( j = n = 0 ; j < i ; j++ ) { CheckForDirectSelfInteractionWithRadii   ; } }
            else            { for ( j = n = 0 ; j < i ; j++ ) { CheckForDirectSelfInteractionWithNoRadii ; } }

            /* . Save the interactions. */
            isOK = SaveInteractions ( pairList, i, n, indices, False ) ;
            if ( ! isOK ) break ;

            /* . Unflag excluded points for i. */
            if ( hasExclusions ) UnflagSelfExclusions ;
        }
    }
    return isOK ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Self-pairlist generation using a grid point-cell method.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Boolean MakeSelfPairListFromCoordinates3PointCell ( PairList *pairList, const Boolean                 sortIndices     ,
                                                                               const Real                    cutOff          ,
                                                                               const Coordinates3           *coordinates3    ,
                                                                               const RealArray1D            *radii           ,
                                                                               const PairConnections        *exclusions      ,
                                                                               const Boolean                *QAND            ,
                                                                               const Boolean                *QOR             ,
                                                                               const RegularGrid            *grid            ,
                                                                               const RegularGridOccupancy   *occupancy       ,
                                                                                     RegularGridSearchRange *gridSearchRange ,
                                                                                     Boolean                *isIncluded      ,
                                                                                     Integer                *indices         )
{
    Boolean  hasExclusions, hasRadii, includeAll, isOK = True, QORI ;
    Integer  c, i, j, m, n, numberOfCells, numberOfPoints, p, pStart, pTotal, x, xFirst = 0, xLast = 0 ;
    Real     cutOffSquared, ri, rij, xij, yij, zij ;

    /* . Initialization. */
    cutOffSquared  = cutOff * cutOff ;
    hasExclusions  = ( exclusions != NULL ) ;
    hasRadii       = ( radii      != NULL ) ;
    numberOfPoints = coordinates3->extent0  ;
    InitializeIsIncluded ( numberOfPoints, QAND ) ;

    /* . Outer loop over indices. */
    for ( i = 1 ; i < numberOfPoints ; i++ )
    {
        /* . AND test. */
        if ( QAND[i] )
        {
            /* . Get some information for the point. */
            QORI = QOR[i] ;
            if ( hasRadii ) ri = cutOff + Array1D_Item ( radii, i ) ;

            /* . Flag excluded points for i. */
            if ( hasExclusions ) FlagSelfExclusions ;

            /* . Find all cells that are within the cutOff of the cell of the current point. */
            numberOfCells = RegularGrid_FindCellsWithinRangeOfPoint ( grid, Array1D_Item ( occupancy->pointCells, i ) ,
                                                                                 Coordinates3_RowPointer ( coordinates3, i ) ,
                                                                                                     gridSearchRange, NULL ) ;

            /* . Loop over cells. */
            for ( m = n = 0 ; m < numberOfCells ; m++ )
            {
                /* . Set some information for the cell. */
                GetCellInformation ( m, occupancy, c, includeAll, pStart, pTotal ) ;

                /* . Generate interactions. */
                     if ( hasRadii   ) { for ( p = 0 ; p < pTotal ; p++ ) { CheckForGridSelfInteractionWithRadii   ( p, pStart ) ; } }
                else if ( includeAll ) { for ( p = 0 ; p < pTotal ; p++ ) { IncludeAllCellSelfInteractions         ( p, pStart ) ; } }
                else                   { for ( p = 0 ; p < pTotal ; p++ ) { CheckForGridSelfInteractionWithNoRadii ( p, pStart ) ; } }
            }

            /* . Save the interactions. */
            isOK = SaveInteractions ( pairList, i, n, indices, sortIndices ) ;
            if ( ! isOK ) break ;

            /* . Unflag excluded points for i. */
            if ( hasExclusions ) UnflagSelfExclusions ;
        }
    }
    return isOK ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Save interactions for a point on the pairlist.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Boolean SaveInteractions ( PairList *pairList, const Integer  i, const Integer  n, const Integer  *indices, const Boolean sortIndices )
{
    Boolean isOK = True ;
    if ( n > 0 )
    {
        auto PairRecord *record ;
        auto Status      localStatus = Status_OK ;
        record = PairRecord_FromIndices ( i, n, indices, &localStatus ) ;
        if ( sortIndices ) PairRecord_Sort ( record ) ;
        PairList_Append ( pairList, record, &localStatus ) ;
        isOK = ( localStatus == Status_OK ) ;
    }
    return isOK ;
}
