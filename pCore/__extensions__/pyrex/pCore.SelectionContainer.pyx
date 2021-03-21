"""Handle selection containers."""

from .CoreError     import CoreError
from .Serialization import RawObjectConstructor

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
_DefaultItemName = "Selection"

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class SelectionContainer:

    def __copy__ ( self ):
        """Copying."""
        cdef CStatus            cStatus = CStatus_OK
        cdef SelectionContainer new
        new         = self.__class__.Raw ( )
        new.cObject = SelectionContainer_Clone ( self.cObject, &cStatus )
        new.isOwner = True
        if cStatus != CStatus_OK: raise CoreError ( "Unable to clone SelectionContainer." )
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            SelectionContainer_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getitem__ ( self, CInteger index ):
        """Get an item."""
        cdef Selection new
        if ( index < 0 ) or ( index >= len ( self ) ):
            raise IndexError
        else:
            new         = Selection.Raw ( )
            new.cObject = self.cObject.items[index]
            new.isOwner = False
            return new

    def __getstate__ ( self ):
        """Return the state."""
        cdef CInteger    i, s
        cdef CSelection *item
        items = []
        n     = len ( self )
        if n > 0:
            for i from 0 <= i < n:
                item    = self.cObject.items[i]
                indices = []
                for s from 0 <= s < item.capacity:
                    indices.append ( item.indices[s] )
                items.append ( indices )
        return { "items" : items }

    def __init__ ( self, capacity, itemName = _DefaultItemName ):
        """Constructor with capacity."""
        self._Initialize ( )
        self._Allocate ( capacity )
        self.itemName = itemName

    def __len__ ( self ):
        """Length."""
        return SelectionContainer_Capacity ( self.cObject )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        cdef CInteger    i, s
        cdef CSelection *item
        cdef CStatus     cStatus = CStatus_OK
        items    = state["items"]
        capacity = len ( items )
        self._Allocate ( capacity )
        if capacity > 0:
            for ( i, indices ) in enumerate ( items ):
                item = Selection_Allocate ( len ( indices ), &cStatus )
                if item == NULL: break
                for ( s, index ) in enumerate ( indices ):
                    item.indices[s] = index
                self.cObject.items[i] = item
            if cStatus != CStatus_OK: raise CoreError ( "Unable to set the state of SelectionContainer." )

    def _Allocate ( self, capacity ):
        """Allocation."""
        self.cObject = SelectionContainer_Allocate ( capacity, NULL )
        self.isOwner = True
        if self.cObject == NULL: raise CoreError ( "Unable to allocate SelectionContainer." )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject  = NULL
        self.isOwner  = False
        self.itemName = _DefaultItemName

    @classmethod
    def FromCapacity ( selfClass, capacity ):
        """Create a selection container consisting of single entry selections up to capacity."""
        cdef SelectionContainer self
        cdef CStatus            cStatus = CStatus_OK
        self         = selfClass.Raw ( )
        self.cObject = SelectionContainer_FromCapacity  ( capacity, &cStatus )
        if cStatus != CStatus_OK: raise CoreError ( "Unable to create SelectionContainer from capacity {:d}.".format ( capacity ) )
        return self

    def FuseItems ( self, BooleanBlock toFuse not None ):
        """Fuse items with indices in toFuse."""
        cdef CStatus cStatus = CStatus_OK
        SelectionContainer_FuseItems ( self.cObject, toFuse.cObject, &cStatus )
        if cStatus != CStatus_OK: raise CoreError ( "Unable to fuse items within SelectionContainer." )

    def MakeMembershipFlags ( self, Selection members not None, andTest = False ):
        """Make membership flags."""
        cdef BooleanBlock   flags
        cdef CBoolean       cANDTest
        cdef CBooleanBlock *cFlags  = NULL
        cdef CStatus        cStatus = CStatus_OK
        if andTest: cANDTest = CTrue
        else:       cANDTest = CFalse
        cFlags = SelectionContainer_MakeMembershipFlags ( self.cObject, members.cObject, cANDTest, &cStatus )
        if cStatus != CStatus_OK: raise CoreError ( "Unable to make SelectionContainer membership flags." )
        flags         = BooleanBlock.Raw ( )
        flags.cObject = cFlags
        return flags

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def SummaryItems ( self ):
        """Summary items."""
        return [ ( "Number of " + self.itemName + "s", "{:d}".format ( len ( self ) ) ) ]

    def UnionOfItems ( self, Selection toUnion not None ):
        """Get the union of items with indices in toUnion."""
        cdef Selection new
        cdef CStatus   cStatus = CStatus_OK
        new         = Selection.Raw ( )
        new.cObject = SelectionContainer_UnionOfItems ( self.cObject, toUnion.cObject, &cStatus )
        new.isOwner = True
        if cStatus != CStatus_OK: raise CoreError ( "Unable to determine the union of items for SelectionContainer." )
        return new

    property capacity:
        def __get__ ( self ): return SelectionContainer_Capacity ( self.cObject )

    property upperBound:
        def __get__ ( self ): return SelectionContainer_UpperBound ( self.cObject )
