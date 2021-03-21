"""PySCF QC model."""

#===================================================================================================================================
# . Modules contributed by Guilherme M. Arantes (USP, 2019).
#
# . Notes:
#
# . Only mean-field methods {'RHF', 'RKS', 'UHF, 'UKS'} are currently supported for pure QC and QCMM energies and gradients.
# . For QCMM gradients, use NBModel.NBModelPySCF.
# 
# . Post-HF methods in PySCF may be applied to the converged target.qcModel.mf (either pure QC or QCMM). 
#
#===================================================================================================================================

from .NBModelPySCF                import NBModelPySCF
from .QCMMElectrostaticModelPySCF import QCMMElectrostaticModelPySCF
from .QCModelPySCF                import QCModelPySCF
