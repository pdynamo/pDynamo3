"""Handle cross and self pair-lists."""

from  collections   import defaultdict
from .CoreError     import CoreError
from .Serialization import RawObjectConstructor

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class PairList:

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            PairList_Deallocate ( &self.cObject )
            self.isOwner = False

    def __getstate__ ( self ):
        """Return the state."""
        state = None
        if len ( self ) > 0:
            items = []
            for ( i, j ) in self: items.extend ( [ i, j ] )
            state = { "items" : items, "shape" : [ len ( items ) // 2, 2 ] }
            if self.label is not None: state["label"] = self.label
        return state

    def __init__ ( self, capacity = 0, label = None ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( capacity )
        if label is not None: self.label = label

    def __iter__ ( self ):
        """Return an iterator."""
        return PairListIterator ( self )

    def __len__ ( self ):
        """The number of pairs."""
        return PairList_NumberOfPairs ( self.cObject )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        if state is not None:
            self._CObjectFromIndexPairs ( state["items"] )
            if "label" in state: self.label = state["label"]

    def _Allocate ( self, CInteger capacity ):
        """Constructor."""
        self.cObject = PairList_Allocate ( capacity, NULL )
        self.isOwner = True
        if self.cObject == NULL: raise CoreError ( "Unable to allocate pair-list." )
        self._PostAllocate ( )

    def _CObjectFromIndexPairs ( self, indices ):
        """Make a list from a list of index pairs."""
        cdef CInteger     i, m, maximumJ, n, numberI
        cdef CInteger    *cIndices
        cdef CPairRecord *cRecord
        cdef CStatus      cStatus = CStatus_OK
        PairList_Deallocate ( &self.cObject )
        if len ( indices ) > 0:
            if len ( indices ) % 2 != 0: raise CoreError ( "Odd number of items in index pair-list." )
            sets = self._MakeSetsFromIndexPairs ( indices )
            ( numberI, maximumJ, records ) = self._MakeRecordsFromSets ( sets )
            cIndices = Integer_Allocate ( maximumJ, &cStatus )
            self._Allocate ( numberI )
            if cIndices != NULL:
                for i in sorted ( records.keys ( ) ):
                    js = records[i]
                    n  = len ( js )
                    for m from 0 <= m < n: cIndices[m] = js[m]
                    cRecord = PairRecord_FromIndices ( i, n, cIndices, &cStatus )
                    PairList_Append ( self.cObject, cRecord, &cStatus )
                    if cStatus != CStatus_OK: break
                Integer_Deallocate ( &cIndices )
            if cStatus != CStatus_OK:
                raise CoreError ( "Unable to make pair-list from index pairs." )
            else:
                self.cObject.isSorted = CTrue

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False
        self.label   = self.defaultLabel

    def _PostAllocate ( self ):
        """Post-allocation."""
        pass

    def _MakeRecordsFromSets ( self, sets ):
        """Make records from sets."""
        maximumJ = 0
        minimumJ = 0
        records  = {}
        for ( i, js0 ) in sets.items ( ):
            js         = sorted ( js0 )
            maximumJ   = max ( maximumJ, len ( js ) )
            minimumJ   = min ( minimumJ, js[0] )
            records[i] = js
        minimumI = min ( records.keys ( ) )
        numberI  = len ( records )
        if ( minimumI < 0 ) or ( minimumJ < 0 ): raise CoreError ( "Negative indices in pair-list generation." )
        return ( numberI, maximumJ, records )

    def _MakeSetsFromIndexPairs ( self, indices ):
        """Make sets from index pairs."""
        sets = defaultdict ( set )
        for ( i, j ) in zip ( indices[::2], indices[1::2] ):
            sets[i].add ( j )
        return sets

    @classmethod
    def FromIndexPairs ( selfClass, indices, label = None ):
        """Constructor from index pairs."""
        self = selfClass.Raw ( )
        state = { "items" : indices, "shape" : [ len ( indices ) // 2, 2 ] }
        if label is not None: state["label"] = label
        self.__setstate__ ( state )
        return self

    @classmethod
    def Merge ( selfClass, items, information = {} ):
        """Merging - slow version."""
        # . None allowed.
        new        = None
        increments = information.get ( "Index Increments", None )
        if increments is not None:
            indices = []
            label   = None
            for ( increment, item ) in zip ( increments, items ):
                if item is not None:
                    state = item.__getstate__ ( )
                    if state is not None:
                        label = state.get ( "label", None )
                        old   = state["items"]
                        for i in old: indices.append ( i + increment )
            n0  = len ( indices ) // 2
            if n0 > 0:
                new   = selfClass.Raw ( )
                state = { "items" : indices, "shape" : [ n0, 2 ] }
                if label is not None: state["label"] = label
                new.__setstate__ ( state )
        return new

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def SummaryItems ( self ):
        """Summary items."""
        return [ ( self.label, "{:d}".format ( len ( self ) ) ) ]

    @property
    def defaultLabel  ( self ): return "Pair-List"
    @property
    def numberOfPairs ( self ): return len ( self )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class CrossPairList ( PairList ):
    """Cross pairlist."""

    @classmethod
    def FromSelfPairList ( selfClass, PairList  pairList      ,
                                      Selection andSelection1 ,
                                      Selection andSelection2 ,
                                      Selection orSelection   , excluded = False ):
        """Make a cross pairlist from a self pairlist."""
        cdef CrossPairList  self
        cdef CPairList     *cPairList0   = NULL
        cdef CPairList     *cPairList    = NULL
        cdef CSelection    *cSelection1  = NULL
        cdef CSelection    *cSelection2  = NULL
        cdef CSelection    *cOrSelection = NULL
        cdef CStatus        status       = CStatus_OK
        if andSelection1 is not None: cSelection1  = andSelection1.cObject
        if andSelection2 is not None: cSelection2  = andSelection2.cObject
        if orSelection   is not None: cOrSelection = orSelection.cObject
        if pairList      is not None: cPairList0   = pairList.cObject
        if excluded:
            cPairList = SelfPairList_ToCrossPairListExcluded ( cPairList0   ,
                                                               cSelection1  ,
                                                               cSelection2  ,
                                                               cOrSelection ,
                                                               &status      )
        else:
            cPairList = SelfPairList_ToCrossPairList ( cPairList0   ,
                                                       cSelection1  ,
                                                       cSelection2  ,
                                                       cOrSelection ,
                                                       &status      )
        if status != CStatus_OK: raise CoreError ( "Error making cross pairlist from self pairlist." )
        self         = selfClass.Raw ( )
        self.cObject = cPairList
        self.isOwner = True
        return self

    @classmethod
    def Full ( selfClass, capacity1, Selection andSelection1, capacity2, Selection andSelection2, excluded = False ):
        """Make a full cross pairlist."""
        cdef CrossPairList  self
        cdef CPairList     *cPairList   = NULL
        cdef CSelection    *cSelection1 = NULL
        cdef CSelection    *cSelection2 = NULL
        cdef CStatus        status      = CStatus_OK
        if andSelection1 is not None: cSelection1 = andSelection1.cObject
        if andSelection2 is not None: cSelection2 = andSelection2.cObject
        if excluded:
            cPairList = CrossPairList_MakeFullExcluded ( capacity1   ,
                                                         cSelection1 ,
                                                         capacity2   ,
                                                         cSelection2 ,
                                                         &status     )
        else:
            cPairList = CrossPairList_MakeFull ( capacity1   ,
                                                 cSelection1 ,
                                                 capacity2   ,
                                                 cSelection2 ,
                                                 &status     )
        if status != CStatus_OK: raise CoreError ( "Error making full cross pairlist." )
        self         = selfClass.Raw ( )
        self.cObject = cPairList
        self.isOwner = True
        return self

    property defaultLabel:
        def __get__ ( self ): return "Cross Pair-List"

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class SelfPairList ( PairList ):
    """Self pairlist."""

    def _MakeSetsFromIndexPairs ( self, indices ):
        """Make sets from index pairs."""
        sets = defaultdict ( set )
        for ( i, j ) in zip ( indices[::2], indices[1::2] ):
            if   i > j: sets[i].add ( j )
            elif i < j: sets[j].add ( i )
        return sets

    def _PostAllocate ( self ):
        """Post-allocation."""
        if self.cObject != NULL: self.cObject.isSelf = CTrue

    @classmethod
    def FromSelfPairList ( selfClass, PairList  pairList     ,
                                                capacity     ,
                                      Selection andSelection ,
                                      Selection orSelection  , excluded = False ):
        """Make a cross pairlist from a self pairlist."""
        cdef SelfPairList  self
        cdef CPairList    *cPairList0    = NULL
        cdef CPairList    *cPairList     = NULL
        cdef CSelection   *cAndSelection = NULL
        cdef CSelection   *cOrSelection  = NULL
        cdef CStatus       status        = CStatus_OK
        if andSelection is not None: cAndSelection = andSelection.cObject
        if orSelection  is not None: cOrSelection  = orSelection.cObject
        if pairList     is not None: cPairList0    = pairList.cObject
        if excluded:
            cPairList = SelfPairList_ToSelfPairListExcluded ( cPairList0    ,
                                                              capacity      ,
                                                              cAndSelection ,
                                                              cOrSelection  ,
                                                              &status       )
        else:
            cPairList = SelfPairList_ToSelfPairList ( cPairList0    ,
                                                      cAndSelection ,
                                                      cOrSelection  ,
                                                      &status       )
        if status != CStatus_OK: raise CoreError ( "Error making self pairlist from self pairlist." )
        self         = selfClass.Raw ( )
        self.cObject = cPairList
        self.isOwner = True
        return self

    def GetConnectedComponents ( self, CInteger upperBound = 0 ):
        """Generate the connected components of the pair-list."""
        cdef CStatus            cStatus = CStatus_OK
        cdef SelectionContainer new
        new         = SelectionContainer.Raw ( )
        new.cObject = SelfPairList_GetConnectedComponents ( self.cObject, upperBound, &cStatus )
        new.isOwner = True
        if cStatus != CStatus_OK: raise CoreError ( "Unable to generate connected components from self pair-list." )
        return new

    def Prune ( self, Selection selection, information = {} ):
        """Pruning."""
        cdef CStatus      cStatus = CStatus_OK
        cdef SelfPairList new
        new         = SelfPairList.Raw ( )
        new.cObject = SelfPairList_ToSelfPairList ( self.cObject, selection.cObject, NULL, &cStatus )
        SelfPairList_Renumber ( new.cObject, selection.cObject, &cStatus )
        new.isOwner = True
        new.label   = self.label 
        if cStatus != CStatus_OK: raise CoreError ( "Unable to prune self pair-list." )
        return new

    property defaultLabel:
        def __get__ ( self ): return "Self Pair-List"

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class PairListIterator:

    def __init__ ( self, PairList target ):
        """Constructor."""
        cdef CPairList *cTarget
        if target is None: cTarget = NULL
        else:              cTarget = target.cObject
        PairListIterator_Initialize ( &self.cIterator, cTarget )
        self.cRecord = PairListIterator_Next ( &self.cIterator )
        self.current = 0

    def __next__ ( self ):
        """Get the next pair."""
        cdef CInteger i, j
        if ( self.cRecord != NULL ) and ( self.current >= self.cRecord.capacity ):
            self.cRecord = PairListIterator_Next ( &self.cIterator )
            self.current = 0
        if self.cRecord == NULL: raise StopIteration
        else:
            i = self.cRecord.index
            j = self.cRecord.indices[self.current]
            self.current += 1
            return ( i, j )
