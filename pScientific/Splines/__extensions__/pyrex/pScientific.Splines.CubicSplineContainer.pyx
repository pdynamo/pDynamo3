"""A container for cubic splines."""

import os, os.path

from  pCore       import Clone                , \
                         RawObjectConstructor
from .SplineError import SplineError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class CubicSplineContainer:

    def __copy__ ( self ):
        """Copying."""
        cdef CubicSplineContainer new
        cdef CStatus cStatus = CStatus_OK
        new               = self.__class__.Raw ( )
        new.cObject       = CubicSplineContainer_Clone ( self.cObject, &cStatus )
        if cStatus != CStatus_OK: raise SplineError ( "Error cloning spline container." )
        new.isOwner       = True
        new.label         = self.label
        new.keys          = Clone ( self.keys          )
        new.uniqueEntries = Clone ( self.uniqueEntries )
        return new

    def __dealloc__ ( self ):
        """Destructor."""
        self.uniqueEntries = None
        if self.isOwner:
            CubicSplineContainer_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getstate__ ( self ):
        """Return the state."""
        state = { "Keys"           : self.keys          ,
                  "Unique Entries" : self.uniqueEntries }
        if self.label is not None: state["Label"] = self.label
        return state

    def __init__ ( self, capacity ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( capacity )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        self._CreateObject ( state["Unique Entries"], state["Keys"] )
        self.label = state.get ( "Label", None )

    def _Allocate ( self, capacity ):
        """Constructor."""
        cdef CStatus cStatus = CStatus_OK
        self.cObject = CubicSplineContainer_Allocate ( capacity, &cStatus )
        if cStatus != CStatus_OK: raise SplineError ( "Error allocating spline container." )
        self.isOwner = True

    def _CreateObject ( self, uniqueEntries, keys ):
        """Create the object."""
        cdef CubicSpline entry
        cdef CInteger    i
        capacity = len ( keys     )
        self._Allocate ( capacity )
        self.keys          = keys
        self.uniqueEntries = uniqueEntries
        for i from 0 <= i < self.cObject.capacity:
            entry = uniqueEntries[keys[i]]
            self.cObject.entries[i] = entry.cObject

    def _Initialize ( self ):
        """Initialization."""
        self.cObject       = NULL
        self.isOwner       = False
        self.keys          = None
        self.label         = None
        self.uniqueEntries = None

    @classmethod
    def FromUniqueEntries ( selfClass, uniqueEntries, keys, label = None ):
        """Constructor from unique entries."""
        cdef CubicSplineContainer self
        self = selfClass.Raw ( )
        self._CreateObject ( uniqueEntries, keys )
        return self

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self
