"""Random number generators."""

import time

from pCore import RawObjectConstructor

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class RandomNumberGenerator:

    # . Public methods.
    def __copy__ ( self ):
        """Copying."""
        cdef RandomNumberGenerator new
        new             = self.__class__.Raw ( )
        new.cObject     = RandomNumberGenerator_Clone ( self.cObject )
        new.initialSeed = self.initialSeed
        new.isOwner     = True
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            RandomNumberGenerator_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getstate__ ( self ):
        """Return the state."""
        return { "initialSeed" : self.initialSeed }

    def __init__ ( self ):
        """Constructor given an iterable."""
        self._Initialize ( )
        self._Allocate   ( )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        self.SetSeed ( state["initialSeed"] )

    def _Allocate ( self ):
        """Allocation."""
        self.cObject = RandomNumberGenerator_Allocate ( RandomNumberGeneratorType_MersenneTwister )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject     = NULL
        self.initialSeed = -1
        self.isOwner     = False

    def NextReal ( self ):
        """Return the next real."""
        return RandomNumberGenerator_NextReal ( self.cObject )

    def NextReals ( self, RealArray1D values ):
        """Return the next reals."""
        cdef CRandomNumberGenerator *randomNumberGenerator
        cdef CRealArray1D           *cValues
        cdef CInteger                 i
        cValues               = values.cObject
        randomNumberGenerator = self.cObject
        for i from 0 <= i < values.extent:
            RealArray1D_SetItem ( cValues, i, RandomNumberGenerator_NextReal ( randomNumberGenerator ), NULL )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def SetSeed ( self, CCardinal seed ):
        """Reset the seed."""
        RandomNumberGenerator_SetSeed ( self.cObject, seed )
        self.initialSeed = seed

    def Test ( self ):
        """Self test."""
        cdef CCardinal n
        cdef CInteger  i
        self.SetSeed ( 4357 )
        for i from 0 <= i < 1000:
            n = RandomNumberGenerator_NextCardinal ( self.cObject )
        return ( n == 1186927261 )

    @classmethod
    def WithRandomSeed ( selfClass ):
        """Constructor given a random seed."""
        cdef CCardinal seed
        seed = time.time ( )
        return selfClass.WithSeed ( seed )

    @classmethod
    def WithSeed ( selfClass, CCardinal seed ):
        """Constructor given a seed."""
        self = selfClass ( )
        self.SetSeed ( seed )
        return self
