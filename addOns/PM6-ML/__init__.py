"""Delta-ML (Machine Learning) correction to QC PM6 energy and gradients."""

#===================================================================================================================================
# . Module contributed by Jan Řezáč (Czech Academy of Sciences, 2025).
#
# . Instalation notes:
#
#     . conda environment with python 3.11:
#         $ conda install -c conda-forge python=3.11.8
#         $ conda install -c conda-forge simple-dftd3 dftd3-python
#         $ conda install -c conda-forge torchmd-net
#     . run pDynamo installation script and source the environment
#
# . Obtain ML model weights at https://github.com/Honza-R/mopac-ml/blob/main/models/PM6-ML_correction_seed8_best.ckpt
#
# . Example: examples/addOns/PM6-ML/DeltaML-PM6.py
#
# 
# . References:
#
#   Guilherme Menegon Arantes, Jan Řezáč (2025)
#   Benchmark of approximate quantum chemical and machine learning potentials for biochemical proton transfer reactions
#   DOI: 10.26434/chemrxiv-2025-xt41p
#
#   Martin Nováček, Jan Řezáč (2025)
#   PM6-ML: The Synergy of Semiempirical Quantum Chemistry and Machine Learning Transformed into a Practical Computational Method
#   DOI: 10.1021/acs.jctc.4c01330
#   J. Chem. Theo. Comput. 2025, 21, 678
#
#===================================================================================================================================

from .QCDeltaMLModel             import QCDeltaMLModel                                    , \
                                        QCDeltaMLModelState

