"""The base class for energy models together with energy model and closure priorities."""

from  enum  import IntEnum
from  pCore import AttributableObject , \
                   Clone              , \
                   logFile            , \
                   LogFileActive      , \
                   SummarizableObject

# . Energy model keys are unique (and the same as system keys).
#
#   New models are compatible with those of lower or equal priority and different key.
#
#   Models with priorities lower than the NullModel are compatible with everything.

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . IntEnum is preferred as it can be ordered automatically.
class EnergyClosurePriority ( IntEnum ):
    """Energy closure priorities."""
    NBInitialization      =  10 # In future use enum.auto ( ) (from Python 3.6+)
    QCInitialization      =  20
    QCIntegrals           =  30
    QCOrthogonalizer      =  40
    QCPreEnergy           =  50
    QCEnergy              =  60
    QCPostEnergy          =  70
    QCPreGradients        =  80
    QCGradients           =  90
    QCFinalization        = 100
    IndependentEnergyTerm = 110

class EnergyModelPriority ( IntEnum ):
    """Energy model priorities."""
    Restraint =   0
    NullModel =  10
    MMModel   =  20
    QCModel   =  30
    QCAddOns  =  40 # . Charge restraints, dispersion, etc.
    NBModel   =  50
    QCMMModel =  60

#===================================================================================================================================
# . State class.
#===================================================================================================================================
class EnergyModelState ( AttributableObject ):
    """The base class for energy model states."""

    _attributable = dict ( AttributableObject._attributable )
    _unpicklable  = set  ( AttributableObject._unpicklable  )
    _attributable.update ( { "target" : None } )
    _unpicklable.add     (   "target"          )

    def SummaryItems ( self ):
        """Summary items."""
        return []

    @classmethod
    def FromTarget ( selfClass, target ):
        """Constructor given a target."""
        self        = selfClass ( )
        self.target = target
        return self

#===================================================================================================================================
# . Model class.
#===================================================================================================================================
class EnergyModel ( SummarizableObject ):
    """The base class for energy models."""

    _classLabel  = "Energy Model"
    _stateName   = None
    _stateObject = None

    def BuildModel ( self, target, **options ):
        """Build the model."""
        state = None
        if self.__class__._stateObject is not None:
            state = self.__class__._stateObject.FromTarget ( target )
            setattr ( target, self.__class__._stateName, state )
        return state

    def ClearScratch ( self, scratch ):
        """Clear scratch."""
        pass

    def EnergyClosures ( self, target ):
        """Return energy closures."""
        return []

    def UnbuildModel ( self, target ):
        """Unbuild the model."""
        self.ClearScratch ( target.scratch )
        if ( self.__class__._stateName is not None ) and \
            hasattr ( target, self.__class__._stateName ):
            delattr ( target, self.__class__._stateName )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
