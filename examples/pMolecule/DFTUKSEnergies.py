"""DFT UKS tests."""

from pMolecule.QCModel import QCModelDFT
from QCTestSystems     import GetRadicalMoleculeSystems , \
                              RunQCTestSet

#===================================================================================================================================
# . Script.
#===================================================================================================================================
_QCModelOptions =  { "LDA:321G"  : { "fitBasis" : "demon" , "functional" : "lda" , "orbitalBasis" : "321g" } }
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
