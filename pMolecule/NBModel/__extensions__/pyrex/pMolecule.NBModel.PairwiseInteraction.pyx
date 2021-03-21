"""Base class for pairwise interactions."""

from pCore import logFile       , \
                  LogFileActive , \
                  RawObjectConstructor

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class PairwiseInteraction:
    """Base class for pairwise interactions."""

    def __copy__ ( self ):
        """Copying."""
        return None

    def __dealloc__ ( self ):
        """Finalization."""
        pass

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getstate__ ( self ):
        """Return the state."""
        return { }
                 
    def __init__ ( self, **options ):
        """Constructor with options."""
        self._Initialize ( )
        self._Allocate   ( )
        self.SetOptions ( **options )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        self._Allocate ( )
        self.SetOptions ( **state )

    def _Allocate ( self ):
        """Allocation."""
        pass

    def _Initialize ( self ):
        """Initialization."""
        pass

    @classmethod
    def WithOptions ( selfClass, **options ):
        """Constructor from options."""
        return selfClass ( **options )

    def Interactions ( self, RealArray1D r not None ):
        """Calculate the interaction components for an array of distances."""
        return { "Electrostatic"   : None ,
                 "Lennard-Jones A" : None ,
                 "Lennard-Jones B" : None }

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def SetOptions ( self, **options ):
        """Set options for the model."""
        pass

    def Summary ( self, log = logFile ):
        """Summary."""
        pass

    def SummaryItems ( self ):
        """Summary items."""
        pass

    @property
    def range ( self ): return 0.0
