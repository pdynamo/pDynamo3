"""Defines a full NB model compatible with the PySCF program."""

from  pMolecule.NBModel           import NBModelFull               , \
                                         QCMMLennardJonesModelFull
from  pScientific                 import Units
from .QCMMElectrostaticModelPySCF import QCMMElectrostaticModelPySCF

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class NBModelPySCF ( NBModelFull ):
    """Defines a full NB model compatible with the PySCF program."""

    _classLabel  = "PySCF NB Model"

    def QCMMModels ( self, qcModel = None, withSymmetry = False ):
        """Default companion QC/MM models for the model."""
        return { "qcmmElectrostatic" : QCMMElectrostaticModelPySCF ,
                 "qcmmLennardJones"  : QCMMLennardJonesModelFull   }

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
