"""A container for block storage."""

import os, os.path

from  pCore              import Clone                , \
                                RawObjectConstructor
from .GaussianBasisError import GaussianBasisError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class BlockStorageContainer:

    def __dealloc__ ( self ):
        """Destructor."""
        if self.isOwner:
            BlockStorageContainer_Deallocate ( &self.cObject )
            self.isOwner = False

    def __init__ ( self, capacity ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( capacity )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def _Allocate ( self, capacity ):
        """Constructor."""
        cdef CStatus cStatus = CStatus_OK
        self.cObject = BlockStorageContainer_Allocate ( capacity, &cStatus )
        if cStatus  != CStatus_OK: raise GaussianBasisError ( "Error allocating spline container." )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self
