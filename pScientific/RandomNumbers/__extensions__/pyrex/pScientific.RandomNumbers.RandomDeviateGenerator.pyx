"""Random number distributions."""

from pCore import RawObjectConstructor

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class NormalDeviateGenerator:

    # . Public methods.
    def __copy__ ( self ):
        """Copying."""
        return self.__class__.FromRandomNumberGenerator ( self.randomNumberGenerator.__copy__ ( ), mu = self.mu, sigma = self.sigma )

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getstate__ ( self ):
        """Return the state."""
        return { "mu"                    : self.mu                    ,
                 "randomNumberGenerator" : self.randomNumberGenerator ,
                 "sigma"                 : self.sigma                 }

    def __init__ ( self ):
        """Constructor given an iterable."""
        self._Initialize ( )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        self.mu                    = state["mu"                   ]
        self.randomNumberGenerator = state["randomNumberGenerator"]
        self.sigma                 = state["sigma"                ]

    def _Initialize ( self ):
        """Initialization."""
        self.mu                    = 0.0
        self.randomNumberGenerator = None
        self.sigma                 = 1.0

    def NextDeviate ( self ):
        """Return the next deviate."""
        return RandomNumberDistribution_GaussianBoxMueller ( self.randomNumberGenerator.cObject, self.mu, self.sigma )

    def NextDeviates ( self, RealArray1D deviates ):
        """Return the next deviates."""
        cdef CRandomNumberGenerator *randomNumberGenerator
        cdef CRealArray1D           *cDeviates
        cdef CInteger                i
        cdef CReal                   mu, sigma
        cDeviates             = deviates.cObject
        mu                    = self.mu
        randomNumberGenerator = self.randomNumberGenerator.cObject
        sigma                 = self.sigma
        for i from 0 <= i < deviates.extent:
            RealArray1D_SetItem ( cDeviates, i, RandomNumberDistribution_GaussianBoxMueller ( randomNumberGenerator, mu, sigma ), NULL )

    def NextStandardDeviate ( self ):
        """Return the next standard deviate."""
        return RandomNumberDistribution_GaussianBoxMueller ( self.randomNumberGenerator.cObject, 0.0, 1.0 )

    def NextStandardDeviates ( self, RealArray1D deviates ):
        """Return the next standard deviates."""
        cdef CRandomNumberGenerator *randomNumberGenerator
        cdef CRealArray1D           *cDeviates
        cdef CInteger                i
        cDeviates             = deviates.cObject
        randomNumberGenerator = self.randomNumberGenerator.cObject
        for i from 0 <= i < deviates.extent:
            RealArray1D_SetItem ( cDeviates, i, RandomNumberDistribution_GaussianBoxMueller ( randomNumberGenerator, 0.0, 1.0 ), NULL )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    @classmethod
    def WithRandomNumberGenerator ( selfClass, randomNumberGenerator, mu = 0.0, sigma = 1.0 ):
        """Constructor given a random number generator."""
        self                       = selfClass ( )
        self.mu                    = mu
        self.randomNumberGenerator = randomNumberGenerator
        self.sigma                 = sigma
        return self
