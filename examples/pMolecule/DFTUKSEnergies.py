"""DFT UKS tests."""

from pMolecule.QCModel import QCModelDFT
from QCTestSystems     import GetRadicalMoleculeSystems , \
                              RunQCTestSet

#===================================================================================================================================
# . Script.
#===================================================================================================================================
_QCModelOptions =  { "LDA:DZVP" : { "fitBasis" : "dgauss-a1-dftjfit" , "functional" : "lda" , "orbitalBasis" : "dgauss-dzvp" } }
testSystems     = GetRadicalMoleculeSystems ( convergerKeywords = { "densityTolerance" : 1.0e-10, "maximumIterations" : 250 } ,
                                              qcModelClass      = QCModelDFT      ,
                                              qcModelOptions    = _QCModelOptions )
RunQCTestSet ( testSystems                                     ,
               "DFT UKS"                                       ,
               maximumEnergyAtoms      = 10                    ,
               maximumEnergyTests      = 20                    ,
               maximumGradientAtoms    = 5                     ,
               maximumGradientTests    = 5                     ,
               referenceDataFileName   = None, #"DFTUKSEnergies.yaml" ,
               testGradients           = True                  )
