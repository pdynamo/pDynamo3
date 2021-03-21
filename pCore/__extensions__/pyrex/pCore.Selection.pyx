"""Selections."""

# . Immutable ordered cardinal sets.

from .CoreError     import CoreError
from .LogFileWriter import logFile, LogFileActive
from .Serialization import RawObjectConstructor

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
_DefaultLabel = "Selection"

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class Selection:

    def __and__ ( Selection self, Selection other ):
        """Intersection."""
        cdef CSelection *old[2]
        cdef CStatus     cStatus = CStatus_OK
        cdef Selection   new
        old[0]      = self.cObject
        old[1]      = other.cObject
        new         = self.__class__.Raw ( )
        new.cObject = Selection_Intersection ( 2, old, &cStatus )
        new.isOwner = True
        if cStatus != CStatus_OK: raise CoreError ( "Unable to create intersection of selections." )
        return new

    def __copy__ ( self ):
        """Copying."""
        cdef CStatus   cStatus = CStatus_OK
        cdef Selection new
        new         = self.__class__.Raw ( )
        new.cObject = Selection_Clone ( self.cObject, &cStatus )
        new.isOwner = True
        new.label   = self.label
        if cStatus != CStatus_OK: raise CoreError ( "Unable to clone Selection." )
        return new

    def __contains__ ( self, CInteger value ):
       """Membership."""
       cdef CBoolean flag
       cdef CStatus  cStatus = CStatus_OK
       flag = Selection_HasItem ( self.cObject, value, &cStatus )
       if cStatus != CStatus_OK: raise CoreError ( "Unable to determine membership of Selection." )
       return ( flag == CTrue )

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            Selection_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getitem__ ( self, CInteger index ):
        """Get an item."""
        if ( index < 0 ) or ( index >= len ( self ) ): raise IndexError
        else: return self.cObject.indices[index]

    def __getstate__ ( self ):
        """Return the state."""
        cdef CInteger i
        items = []
        for i from 0 <= i < self.capacity:
            items.append ( self.cObject.indices[i] )
        return { "items" : items, "label" : self.label }

    def __init__ ( self, iterable ):
        """Constructor given an iterable."""
        cdef CInteger capacity, i, v
        indices  = self._ProcessIterable ( iterable )
        capacity = len ( indices )
        self._Initialize ( )
        self._Allocate ( capacity )
        for ( i, v ) in enumerate ( indices ):
            self.cObject.indices[i] = v

    def __invert__ ( self ):
        """Complement."""
        return self.Complement ( )

    def __len__ ( self ):
        """Length."""
        return self.capacity

    def __or__ ( Selection self, Selection other ):
        """Union."""
        cdef CSelection *old[2]
        cdef CStatus     cStatus = CStatus_OK
        cdef Selection   new
        old[0]      = self.cObject
        old[1]      = other.cObject
        new         = self.__class__.Raw ( )
        new.cObject = Selection_Union ( 2, old, &cStatus )
        new.isOwner = True
        if cStatus != CStatus_OK: raise CoreError ( "Unable to create union of selections." )
        return new

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        cdef CInteger i
        items = self._ProcessIterable ( state["items"] )
        self._Allocate ( len ( items ) )
        for i from 0 <= i < self.capacity:
            self.cObject.indices[i] = items[i]
        self.label = state.get ( "label", _DefaultLabel )

    def __sub__ ( Selection self, Selection other ):
        """Difference."""
        cdef CSelection *cNew
        cdef CStatus     cStatus = CStatus_OK
        cdef Selection   new
        cNew = Selection_Difference ( self.cObject, 1, &other.cObject, &cStatus )
        if cStatus != CStatus_OK: raise CoreError ( "Unable to create difference of selections." )
        if cNew == self.cObject:
            new = self
        else:
            new         = self.__class__.Raw ( )
            new.cObject = cNew
            new.isOwner = True
        return new

    def __xor__ ( Selection self, Selection other ):
        """Symmetric difference."""
        cdef CSelection *old[2]
        cdef CStatus     cStatus = CStatus_OK
        cdef Selection   new
        old[0]      = self.cObject
        old[1]      = other.cObject
        new         = self.__class__.Raw ( )
        new.cObject = Selection_SymmetricDifference ( 2, old, &cStatus )
        new.isOwner = True
        if cStatus != CStatus_OK: raise CoreError ( "Unable to create symmetric difference of selections." )
        return new

    def _Allocate ( self, CInteger capacity ):
        """Allocation."""
        cdef CStatus cStatus = CStatus_OK
        self.cObject = Selection_Allocate ( capacity, &cStatus )
        self.isOwner = True
        if cStatus != CStatus_OK: raise CoreError ( "Unable to allocate Selection." )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False
        self.label   = _DefaultLabel

    @staticmethod
    def _ProcessIterable ( iterable ):
        """Process the indices in an iterable."""
        indices = sorted ( set ( iterable ) )
        if ( len ( indices ) > 0 ) and ( indices[0] < 0 ): raise CoreError ( "Invalid member in Selection." )
        return indices

    def Complement ( self, CInteger upperBound = 0 ):
        """Return a complementary selection."""
        cdef CStatus   cStatus = CStatus_OK
        cdef Selection new
        new         = Selection.Raw ( )
        new.cObject = Selection_Complement ( self.cObject, upperBound, &cStatus )
        new.isOwner = True
        if cStatus != CStatus_OK: raise CoreError ( "Unable to take the complement of Selection." )
        return new

    @classmethod
    def FromIterable ( selfClass, iterable ):
        """Constructor given an iterable."""
        return selfClass ( iterable )

    @classmethod
    def Merge ( selfClass, items, information = {} ):
        """Merging - slow version."""
        # . None allowed.
        new        = None
        increments = information.get ( "Index Increments", None )
        if increments is not None:
            indices = []
            label   = _DefaultLabel
            for ( increment, item ) in zip ( increments, items ):
                if item is not None:
                    state = item.__getstate__ ( )
                    old   = state["items"]
                    label = state.get ( "label", _DefaultLabel )
                    for i in old: indices.append ( i + increment )
            if len ( indices ) > 0:
                new       = selfClass.FromIterable ( indices )
                new.label = label
        return new

    def Position ( self, CInteger value ):
       """Return the position of value in the selection."""
       cdef CInteger position
       cdef CStatus  cStatus = CStatus_OK
       if value in self:
           position = Selection_PositionOfItem ( self.cObject, value, &cStatus )
           if cStatus != CStatus_OK: raise CoreError ( "Unable to determine position in Selection." )
           return position
       else:
           return -1

    def Prune ( self, Selection selection, information = {} ):
        """Pruning."""
        cdef Selection new
        new         = Selection.Raw ( )
        new.cObject = Selection_Prune ( self.cObject, selection.cObject, NULL )
        new.isOwner = True
        new.label   = self.label
        if new.cObject == NULL: raise CoreError ( "Unable to prune Selection." )
        return new

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def Summary ( self, log = logFile ):
        """Summary."""
        n = len ( self )
        if LogFileActive ( log ) and ( n > 0 ):
            log.Paragraph ( "There are {:d} {:s} indices.".format ( n, self.label.lower ( ) ) )

    def SummaryItems ( self ):
        """Summary items."""
        return [ ( self.label, "{:d}".format ( len ( self ) ) ) ]

    property capacity:
        def __get__ ( self ): return Selection_Capacity ( self.cObject )

    property upperBound:
        def __get__ ( self ): return Selection_UpperBound ( self.cObject )
