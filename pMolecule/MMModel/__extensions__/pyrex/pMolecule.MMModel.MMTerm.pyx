"""Base class for MM terms."""

from pCore                 import logFile               , \
                                  LogFileActive         , \
                                  RawObjectConstructor
from pMolecule.EnergyModel import EnergyClosurePriority

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class MMTerm:
    """The base class for MM terms.

    This class should not be used directly.
    """

    # . Public methods.
    def __dealloc__ ( self ):
        """Finalization."""
        pass

    def __init__ ( self, numberOfParameters, numberOfTerms, label = "MM Term" ):
        """Constructor."""
        pass

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def _Initialize ( self ):
        """Initialization."""
        self.label = "MM Term"

    def Energy ( self, Coordinates3 coordinates3, Coordinates3 gradients3 ):
        """Energy and gradients."""
        return 0.0

    def EnergyClosures ( self, target ):
        """Return energy closures."""
        cdef Coordinates3 coordinates3, gradients3
        def f ( ):
            coordinates3 = target.coordinates3
            gradients3   = target.scratch.Get ( "gradients3", None )
            energy       = self.Energy ( coordinates3, gradients3 )
            target.scratch.energyTerms[self.label] = energy
        return [ ( EnergyClosurePriority.IndependentEnergyTerm, f, self.label ) ]

    @staticmethod
    def MergeKeys ( parameterKeys, parameters ):
        """Merge parameter keys."""
        newKeys = set ( parameterKeys )
        if len ( newKeys ) < len ( parameterKeys ):
            # . Keys.
            newKeys = list ( newKeys )
            newKeys.sort ( )
            oldToNew = []
            for oldKey in parameterKeys:
                oldToNew.append ( newKeys.index ( oldKey ) )
            # . Parameters.
            newParameters = []
            for newKey in newKeys:
                newParameters.append ( parameters[parameterKeys.index ( newKey )] )
            # . Finish up.
            return ( True , oldToNew, newKeys      , newParameters )
        else:
            return ( False, None    , parameterKeys, parameters    )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self
