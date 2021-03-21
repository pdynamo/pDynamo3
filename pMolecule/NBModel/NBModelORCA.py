"""Defines a full NB model compatible with the ORCA program."""

from .NBModelFull                import NBModelFull
from .QCMMElectrostaticModelORCA import QCMMElectrostaticModelORCA
from .QCMMLennardJonesModelFull  import QCMMLennardJonesModelFull

# . Notes:
#
#   Check when assign that there is a QC model and that it is ORCA?

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class NBModelORCA ( NBModelFull ):
    """Defines a full NB model compatible with the ORCA program."""

    # . Defaults.
    _classLabel  = "ORCA NB Model"

    def QCMMModels ( self, qcModel = None, withSymmetry = False ):
        """Default companion QC/MM models for the model."""
        return { "qcmmElectrostatic" : QCMMElectrostaticModelORCA ,
                 "qcmmLennardJones"  : QCMMLennardJonesModelFull  }

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
