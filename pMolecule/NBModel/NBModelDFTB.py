"""Defines a full NB model compatible with the DFTB+ program."""

from .NBModelFull                import NBModelFull
from .QCMMElectrostaticModelDFTB import QCMMElectrostaticModelDFTB
from .QCMMLennardJonesModelFull  import QCMMLennardJonesModelFull

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class NBModelDFTB ( NBModelFull ):
    """Defines a full NB model compatible with the DFTB+ program."""

    # . Defaults.
    _classLabel = "DFTB NB Model"

    def QCMMModels ( self, qcModel = None, withSymmetry = False ):
        """Default companion QC/MM models for the model."""
        return { "qcmmElectrostatic" : QCMMElectrostaticModelDFTB ,
                 "qcmmLennardJones"  : QCMMLennardJonesModelFull  }

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
