/*==================================================================================================================================
! . This module handles pair lists.
!=================================================================================================================================*/

# include "BooleanBlock.h"
# include "BooleanUtilities.h"
# include "IntegerUtilities.h"
# include "Integer.h"
# include "Memory.h"
# include "NumericalMacros.h"
# include "PairList.h"

/* . Add checks in Make/To PairList functions for actual versus analytic number of pairs? */

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _GrowthFactor    1.1e+00
# define _MinimumCapacity 32

/*==================================================================================================================================
! . Utility functions.
!=================================================================================================================================*/
/* . Finalize and array. */
static void _AndArrayFinalize ( const Selection *selection, Boolean **and )
{
    if ( selection == NULL ) Boolean_Deallocate ( and ) ;
}

/* . Initialize and array. */
static Boolean *_AndArrayInitialize ( const Integer capacity, Selection *selection, Status *status )
{
    auto Boolean *and = NULL ;
    if ( selection == NULL )
    {
        and = Boolean_Allocate ( capacity, status ) ;
        Boolean_Set ( and, capacity, True ) ;
    }
    else
    {
        auto BooleanBlock *flags = Selection_MakeFlags ( selection, capacity, status ) ;
        if ( flags != NULL ) and = Block_Items ( flags ) ;
    }
    return and ;
}

/* . Finalize index array. */
static void _IndexArrayFinalize ( const Selection *selection, Integer **indices )
{
    if ( selection == NULL ) Integer_Deallocate ( indices ) ;
}

/* . Initialize index array. */
static Integer *_IndexArrayInitialize ( const Integer capacity, const Selection *selection, Status *status )
{
    auto Integer *indices = NULL, s ;
    if ( selection == NULL )
    {
        indices = Integer_Allocate ( capacity, status ) ;
        if ( indices != NULL ) { for ( s = 0 ; s < capacity ; s++ ) indices[s] = s ; }
    }
    else indices = selection->indices ;
    return indices ;
}

/* . Initialize or array. */
static Boolean *_OrArrayInitialize ( const Integer capacity, Selection *selection, Status *status )
{
    auto Boolean *or = NULL ;
    if ( selection != NULL )
    {
        auto BooleanBlock *flags = Selection_MakeFlags ( selection, capacity, status ) ;
        if ( flags != NULL ) or = Block_Items ( flags ) ;
    }
    return or ;
}

/*==================================================================================================================================
! . Pair connection functions.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairConnections *PairConnections_Allocate ( const Integer capacityI, const Integer capacityJ, Status *status )
{
    PairConnections *self = NULL ;
    if ( ( capacityI > 0 ) && ( capacityJ > 0 ) && Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( PairConnections ) ;
        if ( self != NULL )
        {
            self->itemsI = Integer_Allocate ( capacityI + 1, NULL ) ;
            self->itemsJ = Integer_Allocate ( capacityJ    , NULL ) ;
            if ( ( self->itemsI == NULL ) || ( self->itemsJ == NULL ) ) PairConnections_Deallocate ( &self ) ;
            else
            {
                self->capacityI = capacityI ;
                self->capacityJ = capacityJ ;
            }
        }
        if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairConnections_Deallocate ( PairConnections **self )
{
    if ( (*self) != NULL )
    {
        Integer_Deallocate ( &((*self)->itemsI) ) ;
        Integer_Deallocate ( &((*self)->itemsJ) ) ;
        Memory_Deallocate  (   (*self)          ) ;
    }
}

/*==================================================================================================================================
! . Pair excluded functions.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairExcluded *PairExcluded_Allocate ( const Integer capacity, Status *status )
{
    PairExcluded *self = NULL ;
    if ( ( capacity > 0 ) && Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( PairExcluded ) ;
        if ( self != NULL )
        {
            self->indices = Integer_Allocate ( capacity, NULL ) ;
            self->work    = Integer_Allocate ( capacity, NULL ) ;
            if ( ( self->indices == NULL ) ||
                 ( self->work    == NULL ) ) PairExcluded_Deallocate ( &self ) ;
            else self->capacity = capacity ;
        }
        if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairExcluded_Deallocate ( PairExcluded **self )
{
    if ( (*self) != NULL )
    {
        Integer_Deallocate ( &((*self)->indices) ) ;
        Integer_Deallocate ( &((*self)->work   ) ) ;
        Memory_Deallocate  (   (*self)           ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Constructor from indices.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairExcluded *PairExcluded_FromIndices ( const Integer capacity, const Integer *indices, Status *status )
{
    PairExcluded *self = PairExcluded_Allocate ( capacity, status ) ;
    if ( self != NULL )
    {
        Integer_CopyTo ( indices, capacity, self->indices, NULL ) ;
        Integer_Sort   ( self->indices, capacity ) ;
    }
    return self ;
}

/*==================================================================================================================================
! . General pair-list functions.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairList *PairList_Allocate ( const Integer capacity, Status *status )
{
    PairList *self = Memory_AllocateType ( PairList ) ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer n ;
        n = Maximum ( capacity, _MinimumCapacity ) ;
        PairList_Initialize ( self ) ;
        self->records = Memory_AllocateArrayOfReferences ( n, PairRecord ) ;
        if ( self->records == NULL ) PairList_Deallocate ( &self ) ;
        else
        {
            auto Integer i ;
            for ( i = 0 ; i < n ; i++ ) self->records[i] = NULL ;
            self->capacity = n ;
        }
    }
    if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Append a record.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairList_Append ( PairList *self, PairRecord *record, Status *status )
{
    if ( ( self != NULL ) && ( record != NULL ) && Status_IsOK ( status ) )
    {
        auto Boolean isOK = True ;
        if ( self->count >= self->capacity )
        {
            isOK = PairList_Reallocate ( self, ( Integer ) ( self->capacity * _GrowthFactor ), status ) ;
        }
        if ( isOK )
        {
            self->records[self->count] = record ;
            self->count               += 1 ;
            self->numberOfPairs       += record->capacity ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Clear representations.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairList_ClearRepresentations ( PairList *self )
{
    if ( self != NULL ) PairConnections_Deallocate ( &(self->connections) ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairList_Deallocate ( PairList **self )
{
    if ( (*self) != NULL )
    {
        auto Integer i ;
        PairExcluded_Deallocate       ( &((*self)->excluded) ) ;
        PairList_ClearRepresentations (   (*self) ) ;
        for ( i = 0 ; i < (*self)->count ; i++ ) PairRecord_Deallocate ( &((*self)->records[i]) ) ;
        Memory_Deallocate ( (*self)->records ) ;
        Memory_Deallocate ( (*self)          ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get a record.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . No checking and excluded pair-lists must be sorted. */
PairRecord *PairList_GetRecord ( PairList *self, const Integer index )
{
    PairRecord *record = self->records[index] ;
    if ( self->excluded != NULL )
    {
        auto PairExcluded *excluded   = self->excluded ;
        auto PairRecord   *exclusions = record ;
        record        = &(excluded->record) ;
        record->index = exclusions->index ;
        if ( exclusions->capacity == 0 )
        {
            if ( self->isSelf )
            {
                auto Integer i = record->index, j, n ;
                for ( n = 0 ; n < excluded->capacity ; n++ ) { j = excluded->indices[n] ; if ( j >= i ) break ; }
                record->capacity = n ;
            }
            else record->capacity = excluded->capacity ;
            record->indices = excluded->indices ;
        }
        else
        {
            auto Integer c, e, eMaximum, eNext, j, jMaximum, n ;
            eMaximum = excluded->indices[excluded->capacity-1] + 1 ;
            if ( self->isSelf ) jMaximum = record->index ;
            else                jMaximum = eMaximum ;
            e = 1 ; eNext = exclusions->indices[0] ;
            for ( c = n = 0 ; c < excluded->capacity ; c++ )
            {
                j = excluded->indices[c] ;
                if ( j >= jMaximum ) break ;
                else if ( j <  eNext ) { excluded->work[n] = j ; n++ ; }
                else if ( j == eNext )
                {
                    if ( e < exclusions->capacity ) { eNext = exclusions->indices[e] ; e++ ; }
                    else { eNext = eMaximum ; }
                }
            }
            record->capacity = n ;
            record->indices  = excluded->work ;
        }
    }
    return record ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairList_Initialize ( PairList *self )
{
    if ( self != NULL )
    {
        self->isSelf        = False ;
        self->isSorted      = False ;
        self->capacity      = 0 ;
        self->count         = 0 ;
        self->numberOfPairs = 0 ;
        self->connections   = NULL ;
        self->excluded      = NULL ;
        self->records       = NULL ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . The maximum record size.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer PairList_MaximumRecordSize ( const PairList *self )
{
    Integer n = 0 ;
    if ( self != NULL )
    {
        auto Integer r ;
        if ( self->excluded )
        {
            auto PairExcluded *excluded   = self->excluded ;
            auto PairRecord   *exclusions = NULL ;
            for ( r = 0 ; r < self->count ; r++ )
            {
                exclusions = self->records[r] ;
                if ( exclusions->capacity == 0 )
                {
                    if ( self->isSelf )
                    {
                        auto Integer i = exclusions->index, j, m ;
                        for ( m = 0 ; m < excluded->capacity ; m++ ) { j = excluded->indices[m] ; if ( j >= i ) break ; }
                        n = Maximum ( n, m ) ;
                    }
                    else
                    {
                        n = Maximum ( n, excluded->capacity ) ;
                    }
                }
                else
                {
                    auto Integer     c, e, eMaximum, eNext, j, jMaximum, m ;
                    auto PairRecord *record = &(excluded->record) ;
                    eMaximum = excluded->indices[excluded->capacity-1] + 1 ;
                    if ( self->isSelf ) jMaximum = record->index ;
                    else                jMaximum = eMaximum ;
                    e = 1 ; eNext = exclusions->indices[0] ;
                    for ( c = m = 0 ; c < excluded->capacity ; c++ )
                    {
                        j = excluded->indices[c] ;
                        if ( j >= jMaximum ) break ;
                        else if ( j <  eNext ) { excluded->work[m] = j ; m++ ; }
                        else if ( j == eNext )
                        {
                            if ( e < exclusions->capacity ) { eNext = exclusions->indices[e] ; e++ ; }
                            else { eNext = eMaximum ; }
                        }
                    }
                    n = Maximum ( n, m ) ;
                }
            }
        }
        else
        {
            for ( r = 0 ; r < self->count ; r++ ) n = Maximum ( n, self->records[r]->capacity ) ;
        }
    }
    return n ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . The number of pairs.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer PairList_NumberOfPairs ( const PairList *self ) { return ( ( self == NULL ) ? 0 : self->numberOfPairs ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . The number of (active) records.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer PairList_NumberOfRecords ( const PairList *self ) { return ( ( self == NULL ) ? 0 : self->count ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Reallocate records.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The list can never be shortened below the current value of count. */
Boolean PairList_Reallocate ( PairList *self, const Integer capacity, Status *status )
{
    Boolean isOK = True ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer n ;
             if ( capacity > self->capacity ) n = capacity    ;
        else if ( capacity < self->count    ) n = self->count ;
        n = Maximum ( capacity, _MinimumCapacity ) ;
        if ( n != self->capacity )
        {
            auto PairRecord **records  ;
            records = Memory_ReallocateArrayOfReferences ( self->records, n, PairRecord ) ;
            if ( records != NULL )
            {
                if ( n > self->count ) { auto Integer i ; for ( i = self->count ; i < n ; i++ ) records[i] = NULL ; }
                self->capacity = n ;
                self->records  = records ;
            }
            else isOK = False ;
        }
    }
    return isOK ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Sorting.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Integer _RecordCompare ( const void *vSelf, const void *vOther )
{
    Integer result ;
    PairRecord **self, **other ;
    self  = ( PairRecord ** ) vSelf  ;
    other = ( PairRecord ** ) vOther ;
         if ( (*self)->index < (*other)->index ) result = -1 ;
    else if ( (*self)->index > (*other)->index ) result =  1 ;
    else result = 0 ;
    return result ;
}

void PairList_Sort ( PairList *self )
{
    if ( ( self != NULL ) && ( self->count > 1 ) && ( ! self->isSorted ) )
    {
        auto Integer i ;
        qsort ( ( void * ) self, ( Size ) self->count, SizeOf ( PairRecord * ), ( void * ) _RecordCompare ) ;
        for ( i = 0 ; i < self->count ; i++ ) Integer_Sort ( self->records[i]->indices, self->capacity ) ;
        self->isSorted = True ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . The upper bound for the i interactions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . This function is really for internal use to ensure that a list is in the correct format for specific operations. */
Integer PairList_UpperBound ( const PairList *self, const Boolean isSelf )
{
    Integer upperBound = 0 ;
    if ( ( self != NULL ) && ( self->count > 0 ) && ( self->isSelf == isSelf ) && ( self->isSorted ) ) upperBound = self->records[self->count-1]->index + 1 ;
    return upperBound ;
}

/*==================================================================================================================================
! . Cross pair-list functions.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the connection representation of the pair-list.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairConnections *CrossPairList_MakeConnections ( PairList *self, const Integer upperBound, Status *status )
{
    PairConnections *connections = NULL ;
    Integer          upper       = PairList_UpperBound ( self, False );
    if ( ( upper > 0 ) && Status_IsOK ( status ) )
    {
        connections = self->connections ;
        if ( ( connections == NULL ) || ( connections->capacityI < upperBound ) )
        {
            PairConnections_Deallocate ( &(self->connections) ) ;
            connections = PairConnections_Allocate ( Maximum ( upper, upperBound ), self->numberOfPairs, status ) ;
            if ( connections != NULL )
            {
                auto Integer    i, m, n, r ;
                auto PairRecord *pairRecord ;
                Integer_Set ( connections->itemsI, connections->capacityI + 1, 0 ) ;
                for ( n = r = 0 ; r < self->count ; r++ )
                {
                    pairRecord = self->records[r] ;
                    i = pairRecord->index    ;
                    m = pairRecord->capacity ;
                    connections->itemsI[i] = m ;
                    Integer_CopyTo ( pairRecord->indices, m, &(connections->itemsJ[n]), NULL ) ;
                    n += m ;
                }
                for ( i = n = 0 ; i < connections->capacityI ; i++ )
                {
                    m = connections->itemsI[i] ;
                    connections->itemsI[i] = n ;
                    n += m ;
                }
                connections->itemsI[connections->capacityI] = n ;
                self->connections = connections ;
            }
        }
    }
    return connections ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a full pair-list given index information.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern PairList *CrossPairList_MakeFull ( const Integer    capacity1     ,
                                                Selection *andSelection1 ,
                                          const Integer    capacity2     ,
                                                Selection *andSelection2 ,
                                                Status    *status        )
{
    PairList *self = PairList_Allocate ( capacity1, status ) ;
    if ( ( self != NULL ) && ( capacity1 > 0 ) && ( capacity2 > 0 ) )
    {
        auto Boolean  isOK ;
        auto Integer *indices1, *indices2, s ;
        indices1 = _IndexArrayInitialize ( capacity1, andSelection1, status ) ;
        indices2 = _IndexArrayInitialize ( capacity2, andSelection2, status ) ;
        isOK     = ( indices1 != NULL ) && ( indices2 != NULL ) ;
        if ( isOK )
        {
            auto PairRecord *record ;
            auto Status      localStatus = Status_OK ;
            for ( s = 0 ; s < capacity1 ; s++ )
            {
                record = PairRecord_FromIndices ( indices1[s], capacity2, indices2, &localStatus ) ;
                PairList_Append ( self, record, &localStatus ) ;
                isOK = Status_IsValueOK ( localStatus ) ;
                if ( ! isOK ) break ;
            }
            self->isSorted = True ;
            if ( ! isOK ) Status_Set ( status, localStatus ) ;
        }
        _IndexArrayFinalize ( andSelection1, &indices1 ) ;
        _IndexArrayFinalize ( andSelection2, &indices2 ) ;
        if ( ! isOK ) PairList_Deallocate ( &self ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a full pair-list given index information.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern PairList *CrossPairList_MakeFullExcluded ( const Integer    capacity1     ,
                                                        Selection *andSelection1 ,
                                                  const Integer    capacity2     ,
                                                        Selection *andSelection2 ,
                                                        Status    *status        )
{
    PairList *self = PairList_Allocate ( capacity1, status ) ;
    if ( ( self != NULL ) && ( capacity1 > 0 ) && ( capacity2 > 0 ) )
    {
        auto Boolean  isOK ;
        auto Integer *indices1, *indices2, s ;
        indices1 = _IndexArrayInitialize ( capacity1, andSelection1, status ) ;
        indices2 = _IndexArrayInitialize ( capacity2, andSelection2, status ) ;
        self->excluded = PairExcluded_FromIndices ( capacity2, indices2, status ) ;
        isOK = ( indices1 != NULL ) && ( indices2 != NULL ) && ( self->excluded != NULL ) ;
        if ( isOK )
        {
            auto PairRecord *record ;
            auto Status      localStatus = Status_OK ;
            for ( s = 0 ; s < capacity1 ; s++ )
            {
                record = PairRecord_FromIndices ( indices1[s], 0, NULL, &localStatus ) ;
                PairList_Append ( self, record, &localStatus ) ;
                isOK = Status_IsValueOK ( localStatus ) ;
                if ( ! isOK ) break ;
            }
            self->isSorted      = True ;
            self->numberOfPairs = ( capacity1 * capacity2 ) ;
            if ( ! isOK ) Status_Set ( status, localStatus ) ;
        }
        _IndexArrayFinalize ( andSelection1, &indices1 ) ;
        _IndexArrayFinalize ( andSelection2, &indices2 ) ;
        if ( ! isOK ) PairList_Deallocate ( &self ) ;
    }
    return self ;
}

/*==================================================================================================================================
! . Self pair-list functions.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the connected components of the pair-list.
! . Upperbound is needed here for the case where there are connected components consisting of a single index
! . (i.e. those that are absent from the pair-list).
!---------------------------------------------------------------------------------------------------------------------------------*/
SelectionContainer *SelfPairList_GetConnectedComponents ( PairList *self, const Integer upperBound, Status *status )
{
    SelectionContainer *new         = NULL ;
    PairConnections    *connections = SelfPairList_MakeConnections ( self, upperBound, status ) ;
    if ( connections != NULL )
    {
        auto Boolean   *isAssigned, isOK ;
        auto Integer   *indicesI, *indicesJ, n ;
        auto Selection *selection ;
        /* . Allocate space. */
        n           = connections->capacityI ;
        isAssigned  = Boolean_Allocate ( n    , status ) ;
        indicesI    = Integer_Allocate ( n + 1, status ) ;
        indicesJ    = Integer_Allocate ( n    , status ) ;
        /* . Check for memory. */
        isOK = ( indicesI != NULL ) && ( indicesJ != NULL) && ( isAssigned != NULL ) && Status_IsOK ( status ) ;
        if ( isOK )
        {
            auto Integer c, i, j, numberOfComponents, s, start ;
            /* . Initialization. */
            Boolean_Set ( isAssigned, n, False ) ;
            /* . Loop over all indices. */
            for ( n = numberOfComponents = s = 0 ; s < connections->capacityI ; s++ )
            {
                if ( ! isAssigned[s] )
                {
                    /* . Start the new isolate. */
                    indicesJ[n]   = s    ;
                    isAssigned[s] = True ;
                    start         = n    ;
                    n++ ;
                    /* . Assign all indices in the new isolate. */
                    for ( i = start ; i < n ; i++ )
                    {
                        for ( c = connections->itemsI[indicesJ[i]] ; c < connections->itemsI[indicesJ[i]+1] ; c++ )
                        {
                            j = connections->itemsJ[c] ;
                            if ( ! isAssigned[j]  )
                            {
                                indicesJ[n]   = j    ;
                                isAssigned[j] = True ;
                                n++ ;
                            }
                        }
                    }
                    /* . Create the new isolate. */
                    if ( n > start )
                    {
                        indicesI[numberOfComponents] = start ;
                        numberOfComponents++ ;
                    }
                }
            }
            indicesI[numberOfComponents] = n ;
            /* . Create the isolates. */
            new = SelectionContainer_Allocate ( numberOfComponents, status ) ;
            if ( new != NULL )
            {
                for ( s = 0 ; s < numberOfComponents ; s++ )
                {
                    n = indicesI[s+1] - indicesI[s] ;
                    selection = Selection_FromIntegers ( n, &(indicesJ[indicesI[s]]), status ) ;
                    if ( selection == NULL ) { isOK = False ; break ; }
                    new->items[s] = selection ;
                }
            }
            else isOK = False ;
        }
        /* . Finish up. */
        Boolean_Deallocate  ( &isAssigned ) ;
        Integer_Deallocate ( &indicesI   ) ;
        Integer_Deallocate ( &indicesJ   ) ;
        if ( ! isOK ) SelectionContainer_Deallocate ( &new ) ;
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the connection representation of the pair-list.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairConnections *SelfPairList_MakeConnections ( PairList *self, const Integer upperBound, Status *status )
{
    PairConnections *connections = NULL ;
    Integer          upper       = PairList_UpperBound ( self, True ) ;
    if ( ( upper > 0 ) && Status_IsOK ( status ) )
    {
        connections = self->connections ;
        if ( ( connections == NULL ) || ( connections->capacityI < upperBound ) )
        {
            auto Integer *itemsN ;
            PairConnections_Deallocate ( &(self->connections) ) ;
            upper       = Maximum ( upper, upperBound ) ;
            connections = PairConnections_Allocate ( upper, 2 * self->numberOfPairs, status ) ;
            itemsN      = Integer_Allocate        ( upper,                          status ) ;
            if ( ( connections != NULL ) && ( itemsN != NULL ) )
            {
                auto Integer     i, j, m, n, r ;
                auto PairRecord *pairRecord ;
                Integer_Set ( connections->itemsI, connections->capacityI + 1, 0 ) ;
                Integer_Set (              itemsN, connections->capacityI    , 0 ) ;
                for ( r = 0 ; r < self->count ; r++ )
                {
                    pairRecord = self->records[r]  ;
                    i          = pairRecord->index ;
       	            for ( m = 0 ; m < pairRecord->capacity ; m++ )
	            {
	                j          = pairRecord->indices[m] ;
                        itemsN[i] += 1 ;
                        itemsN[j] += 1 ;
                    }
                }
                for ( i = n = 0 ; i < connections->capacityI ; i++ )
                {
                    n += itemsN[i] ;
                    connections->itemsI[i+1] = n ;
                    itemsN[i]                = 0 ;
                }
                for ( r = 0 ; r < self->count ; r++ )
                {
                    pairRecord = self->records[r]  ;
                    i          = pairRecord->index ;
       	            for ( m = 0 ; m < pairRecord->capacity ; m++ )
	            {
	                j = pairRecord->indices[m] ;
                        connections->itemsJ[connections->itemsI[i]+itemsN[i]] = j ;
                        connections->itemsJ[connections->itemsI[j]+itemsN[j]] = i ;
                        itemsN[i] += 1 ;
                        itemsN[j] += 1 ;
                    }
                }
                self->connections = connections ;
            }
            Integer_Deallocate ( &itemsN ) ;
        }
    }
    return connections ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Renumbering based on an input selection.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SelfPairList_Renumber ( PairList  *self    ,
                             Selection *mapping ,
                             Status    *status  )
{
    Integer upper = PairList_UpperBound ( self, True ) ;
    if ( ( mapping != NULL ) && ( upper > 0 ) )
    {
        IntegerBlock *positions = Selection_MakePositions ( mapping, upper, status ) ;
        if ( positions != NULL )
        {
            auto Integer     i, j, m, r ;
            auto Integer    *indices ;
            auto PairRecord *record  ;
            indices = Block_Items ( positions ) ;
            for ( r = 0 ; r < self->count ; r++ )
            {
                record = self->records[r]  ;
                i      = record->index ;
                record->index = indices[i] ;
       	        for ( m = 0 ; m < record->capacity ; m++ )
	        {
	            j = record->indices[m] ;
                    record->indices[m] = indices[j] ;
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Convert a self pair-list into a cross pair-list.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairList *SelfPairList_ToCrossPairList ( PairList  *self          ,
                                         Selection *andSelection1 ,
                                         Selection *andSelection2 ,
                                         Selection *orSelection   ,
                                         Status    *status        )
{
    Integer          upper       = PairList_UpperBound ( self, True ) ;
    PairList        *new         = NULL ;
    PairConnections *connections = SelfPairList_MakeConnections ( self, upper, status ) ;
    if ( ( connections != NULL ) && Status_IsOK ( status ) )
    {
        auto Boolean *and1 = NULL, *and2 = NULL, isOK, *or = NULL, orTest ;
        auto Integer *indices ;
        /* . Check the AND selections. */
        and1 = _AndArrayInitialize ( upper, andSelection1, status ) ;
        and2 = _AndArrayInitialize ( upper, andSelection2, status ) ;
        /* . Check the OR selection. */
        orTest = ( orSelection != NULL ) ;
        if ( orTest ) or = _OrArrayInitialize ( upper, orSelection, status ) ;
        /* . Create the pair-list and other temporary space. */
        indices = Integer_Allocate ( upper, status ) ;
        new     = PairList_Allocate ( upper, status ) ;
        /* . Check for a memory error. */
        isOK = ( and1 != NULL ) && ( and2 != NULL ) && ( indices != NULL ) && ( new != NULL ) && ( ( ! orTest ) || ( orTest && ( or != NULL ) ) ) ;
        if ( isOK )
        {
            auto Integer     count, i, j, m, n ;
            auto PairRecord *record ;
            auto Status      localStatus = Status_OK ;
            /* . Iterate over the list. */
            for ( i = 0 ; i < Minimum ( connections->capacityI, upper ) ; i++ )
            {
                if ( and1[i] )
                {
                    /* . Loop over the indices. */
       	            for ( m = connections->itemsI[i], n = 0 ; m < connections->itemsI[i+1] ; m++ )
	            {
                        j = connections->itemsJ[m] ;
                        if ( and2[j] ) { indices[n] = j ; n++ ; }
                    }
                    count = n ;
                    /* . Apply OR test. */
                    if ( orTest && ( ! or[i] ) )
                    {
                        /* . Loop over the indices. */
       	                for ( m = n = 0 ; m < count ; m++ )
	                {
                            j = indices[m] ;
                            if ( or[j] ) { indices[n] = j ; n++ ; }
                        }
                        count = n ;
                    }
                    /* . Save the data. */
                    if ( count > 0 )
                    {
                        record = PairRecord_FromIndices ( i, count, indices, &localStatus ) ;
                        PairRecord_Sort ( record ) ;
                        PairList_Append ( new, record, &localStatus ) ;
                        isOK = Status_IsValueOK ( localStatus ) ;
                        if ( ! isOK ) break ;
                    }
                }
            }
            new->isSorted = True ;
            if ( ! isOK ) Status_Set ( status, localStatus ) ;
        }
        /* . Finish up. */
        Integer_Deallocate ( &indices ) ;
        _AndArrayFinalize ( andSelection1, &and1 ) ;
        _AndArrayFinalize ( andSelection2, &and2 ) ;
        if ( ! isOK ) PairList_Deallocate ( &new ) ;
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Convert a self pair-list into an excluded cross pair-list.
! . Both the AND selections have to be present.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairList *SelfPairList_ToCrossPairListExcluded ( PairList  *self          ,
                                                 Selection *andSelection1 ,
                                                 Selection *andSelection2 ,
                                                 Selection *orSelection   ,
                                                 Status    *status        )
{
    Boolean   hasSelf = ( self != NULL ) ;
    PairList *new     = NULL ;
    if ( ( ( ! hasSelf ) || ( hasSelf && ( self->isSelf ) ) ) && ( andSelection1 != NULL ) && ( andSelection2 != NULL ) && Status_IsOK ( status ) )
    {
        auto Boolean         *and2, isOK, *or = NULL, orTest ;
        auto Integer          capacity1, capacity2, *indices, *indices1, *indices2, upper ;
        auto PairConnections *connections ;
        upper            = Maximum ( PairList_UpperBound  ( self,    True ) ,                         
                           Maximum ( Selection_UpperBound ( andSelection1 ) ,                         
                                     Selection_UpperBound ( andSelection2 ) ) ) ;                     
        connections      = SelfPairList_MakeConnections ( self, upper, status ) ;
        capacity1        = Selection_Capacity ( andSelection1 ) ; indices1 = andSelection1->indices ;
        capacity2        = Selection_Capacity ( andSelection2 ) ; indices2 = andSelection2->indices ;
        and2             = _AndArrayInitialize ( upper, andSelection2, status ) ;
        orTest           = ( orSelection != NULL ) ;
        if ( orTest ) or = _OrArrayInitialize ( upper, orSelection, status ) ;
        indices          = Integer_Allocate ( 2 * capacity2, status ) ;
        new              = PairList_Allocate (     capacity1, status ) ;
        new->excluded    = PairExcluded_FromIndices ( capacity2, indices2, status ) ;
        isOK = ( and2 != NULL ) && ( ( ! hasSelf ) || ( hasSelf && ( connections != NULL ) ) ) &&
               ( indices != NULL ) && ( new != NULL ) && ( new->excluded != NULL ) &&
               ( ( ! orTest ) || ( orTest && ( or != NULL ) ) ) ;
        if ( isOK )
        {
            auto Integer     count, i, j, m, n, r ;
            auto PairRecord *record ;
            auto Status      localStatus = Status_OK ;
            for ( r = 0 ; r < capacity1 ; r++ )
            {
                i = indices1[r] ;
                n = 0 ;
                if ( hasSelf )
                {
                    for ( m = connections->itemsI[i] ; m < connections->itemsI[i+1] ; m++ )
	            {
                        j = connections->itemsJ[m] ;
                        if ( and2[j] ) { indices[n] = j ; n++ ; }
                    }
                }
                if ( orTest && ( ! or[i] ) )
                {
                    for ( m = 0 ; m < capacity2 ; m++ )
                    {
                        j = indices2[m] ;
                        if ( ! or[j] ) { indices[n] = j ; n++ ; }
                    }
                }
                count = Integer_SortUnique ( indices, n ) ;
                if ( count < capacity2 )
                {
                    record = PairRecord_FromIndices ( i, count, indices, &localStatus ) ;
                    PairList_Append ( new, record, &localStatus ) ;
                    isOK = Status_IsValueOK ( localStatus ) ;
                    if ( ! isOK ) break ;
                }
            }
            new->isSorted      = True ;
            new->numberOfPairs = ( ( new->count * capacity2 ) - new->numberOfPairs ) ;
            if ( ! isOK ) Status_Set ( status, localStatus ) ;
        }
        /* . Finish up. */
        Integer_Deallocate ( &indices ) ;
        if ( ! isOK ) PairList_Deallocate ( &new ) ;
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Convert a self pair-list into another one.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairList *SelfPairList_ToSelfPairList ( PairList  *self         ,
                                        Selection *andSelection ,
                                        Selection *orSelection  ,
                                        Status    *status       )
{
    PairList *new   = NULL ;
    Integer   upper = PairList_UpperBound ( self, True ) ;
    if ( ( upper > 0 ) && Status_IsOK ( status ) )
    {
        auto Boolean *and = NULL, isOK, *or = NULL, orTest ;
        auto Integer *indices ;
        /* . Check the AND selection. */
        and = _AndArrayInitialize ( upper, andSelection, status ) ;
        /* . Check the OR selection. */
        orTest = ( orSelection != NULL ) ;
        if ( orTest ) or = _OrArrayInitialize ( upper, orSelection, status ) ;
        /* . Create the pair-list and other temporary space. */
        new     = PairList_Allocate ( self->count, status ) ;
        indices = Integer_Allocate ( upper - 1  , status ) ;
        /* . Check for a memory error. */
        isOK = ( and != NULL ) && ( indices != NULL ) && ( new != NULL ) && ( ( ! orTest ) || ( orTest && ( or != NULL ) ) ) ;
        if ( isOK )
        {
            auto Integer     count, i, j, m, n, r ;
            auto PairRecord *newRecord, *oldRecord ;
            auto Status      localStatus = Status_OK ;
            /* . Iterate over the list. */
            for ( r = 0 ; r < self->count ; r++ )
            {
                oldRecord = self->records[r] ;
                i = oldRecord->index ;
                if ( and[i] )
                {
                    /* . Loop over the indices. */
       	            for ( m = n = 0 ; m < oldRecord->capacity ; m++ )
	            {
                        j = oldRecord->indices[m] ;
                        if ( and[j] ) { indices[n] = j ; n++ ; }
                    }
                    count = n ;
                    /* . Apply OR test. */
                    if ( orTest && ( ! or[i] ) )
                    {
                        /* . Loop over the indices. */
       	                for ( m = n = 0 ; m < count ; m++ )
	                {
                            j = indices[m] ;
                            if ( or[j] ) { indices[n] = j ; n++ ; }
                        }
                        count = n ;
                    }
                    /* . Save the data. */
                    if ( count > 0 )
                    {
                        newRecord = PairRecord_FromIndices ( i, count, indices, &localStatus ) ;
                        PairRecord_Sort ( newRecord ) ;
                        PairList_Append ( new, newRecord, &localStatus ) ;
                        isOK = Status_IsValueOK ( localStatus ) ;
                        if ( ! isOK ) break ;
                   }
                }
            }
            new->isSelf   = True ;
            new->isSorted = True ;
            if ( ! isOK ) Status_Set ( status, localStatus ) ;
        }
        /* . Finish up. */
        Integer_Deallocate ( &indices ) ;
        _AndArrayFinalize ( andSelection, &and ) ;
        if ( ! isOK ) PairList_Deallocate ( &new ) ;
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Convert a self pair-list into an excluded self pair-list.
! . The AND selection is optional if all indices are to be included.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairList *SelfPairList_ToSelfPairListExcluded ( PairList  *self         ,
                                                Integer    capacity     ,
                                                Selection *andSelection ,
                                                Selection *orSelection  ,
                                                Status    *status       )
{
    Boolean   hasSelf = ( self != NULL ) ;
    PairList *new     = NULL ;
    if ( ( ( ! hasSelf ) || ( hasSelf && ( self->isSelf ) ) ) && ( capacity > 0 ) && Status_IsOK ( status ) )
    {
        auto Boolean         *and, isOK, *or = NULL, orTest ;
        auto Integer         *indices, *indicesE ;
        auto PairConnections *connections ;
        if ( andSelection != NULL ) capacity = Selection_Capacity ( andSelection ) ;
        connections      = SelfPairList_MakeConnections ( self, capacity, status ) ;
        and              = _AndArrayInitialize   ( capacity, andSelection, status ) ;
        indices          = _IndexArrayInitialize ( capacity, andSelection, status ) ;
        orTest           = ( orSelection != NULL ) ;
        if ( orTest ) or = _OrArrayInitialize ( capacity, orSelection, status ) ;
        indicesE         = Integer_Allocate ( 2 * capacity, status ) ;
        new              = PairList_Allocate (     capacity, status ) ;
        new->excluded    = PairExcluded_FromIndices ( capacity, indices, status ) ;
        isOK = ( and != NULL ) && ( ( ! hasSelf ) || ( hasSelf && ( connections != NULL ) ) ) &&
               ( indices != NULL ) && ( new != NULL ) && ( new->excluded != NULL ) &&
               ( ( ! orTest ) || ( orTest && ( or != NULL ) ) ) ;
        if ( isOK )
        {
            auto Integer     count, i, j, m, n, p = 0, r ;
            auto PairRecord *record ;
            auto Status      localStatus = Status_OK ;
            for ( r = 0 ; r < capacity ; r++ )
            {
                i = indices[r] ;
                n = 0 ;
                if ( hasSelf )
                {
       	            for ( m = connections->itemsI[i] ; m < connections->itemsI[i+1] ; m++ )
	            {
                        j = connections->itemsJ[m] ;
                        if ( ( j < i ) && and[j] ) { indicesE[n] = j ; n++ ; }
                    }
                }
                if ( orTest && ( ! or[i] ) )
                {
                    for ( m = 0 ; m < capacity ; m++ )
                    {
                        j = indices[m] ;
                        if ( ( j < i ) &&  ( ! or[j] ) ) { indicesE[n] = j ; n++ ; }
                    }
                }
                count = Integer_SortUnique ( indicesE, n ) ;
                if ( count < r ) /* . The maximum number of interactions for i is r. */
                {
                    p     += r ;
                    record = PairRecord_FromIndices ( i, count, indicesE, &localStatus ) ;
                    PairList_Append ( new, record, &localStatus ) ;
                    isOK = Status_IsValueOK ( localStatus ) ;
                    if ( ! isOK ) break ;
                }
            }
            new->isSelf        = True ;
            new->isSorted      = True ;
            new->numberOfPairs = ( p - new->numberOfPairs ) ;
            if ( ! isOK ) Status_Set ( status, localStatus ) ;
        }
        /* . Finish up. */
        _AndArrayFinalize   ( andSelection, &and     ) ;
        _IndexArrayFinalize ( andSelection, &indices ) ;
        Integer_Deallocate ( &indicesE ) ;
        if ( ! isOK ) PairList_Deallocate ( &new ) ;
    }
    return new ;
}

/*==================================================================================================================================
! . Pair-list iterator.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairListIterator_Initialize ( PairListIterator *self, PairList *target )
{
    if ( self != NULL )
    {
        self->current = 0    ;
        self->target  = NULL ;
        if ( ( target != NULL ) && ( target->count > 0 ) ) self->target = target ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Next iteration.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairRecord *PairListIterator_Next ( PairListIterator *self )
{
    PairRecord *next = NULL ;
    if ( ( self != NULL ) && ( self->target != NULL ) )
    {
        if ( self->current < self->target->count )
        {
            next = PairList_GetRecord ( self->target, self->current ) ; /* self->target->records[self->current] ; */
            self->current += 1 ;
        }
    }
    return next ;
}

/*==================================================================================================================================
! . Pair record functions.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairRecord *PairRecord_Allocate ( const Integer capacity, Status *status )
{
    PairRecord *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( PairRecord ) ;
        if ( self != NULL )
        {
            PairRecord_Initialize ( self ) ;
            if ( capacity > 0 )
            {
                self->indices = Integer_Allocate ( capacity, status ) ;
                if ( self->indices == NULL ) PairRecord_Deallocate ( &self ) ;
                else self->capacity = capacity ;
            }
            if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairRecord_Deallocate ( PairRecord **self )
{
    if ( (*self) != NULL )
    {
        Integer_Deallocate ( &((*self)->indices) ) ;
        Memory_Deallocate   (   (*self)           ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Constructor from indices.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairRecord *PairRecord_FromIndices ( const Integer index, const Integer capacity, const Integer *indices, Status *status )
{
    PairRecord *self = PairRecord_Allocate ( capacity, status ) ;
    if ( self != NULL )
    {
        Integer_CopyTo ( indices, capacity, self->indices, NULL ) ;
        self->index = index ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairRecord_Initialize ( PairRecord *self )
{
    if ( self != NULL )
    {
        self->capacity = 0 ;
        self->index    = 0 ;
        self->indices  = NULL ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Sorting.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairRecord_Sort ( PairRecord *self ) { if ( self != NULL ) Integer_Sort ( self->indices, self->capacity ) ; }
