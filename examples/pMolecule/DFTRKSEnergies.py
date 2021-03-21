"""DFT RKS tests."""

from pMolecule.QCModel import QCModelDFT
from QCTestSystems     import GetClosedShellMoleculeSystems , \
                              RunQCTestSet

#===================================================================================================================================
# . Script.
#===================================================================================================================================
_QCModelOptions =  { "LDA:321G"  : { "fitBasis" : "demon"  , "functional" : "lda" , "orbitalBasis" : "321g" } ,
                     "BLYP:321G" : { "fitBasis" : "demon"  , "functional" : "blyp", "orbitalBasis" : "321g" } ,
                     "LDA:SVP"   : { "fitBasis" : "weigend", "functional" : "lda" , "orbitalBasis" : "svp"  } ,
                     "BLYP:SVP"  : { "fitBasis" : "weigend", "functional" : "blyp", "orbitalBasis" : "svp"  } }
testSystems     = GetClosedShellMoleculeSystems ( convergerKeywords = { "densityTolerance" : 1.0e-10, "maximumIterations" : 250 } ,
                                                  qcModelClass      = QCModelDFT      ,
                                                  qcModelOptions    = _QCModelOptions )
RunQCTestSet ( testSystems                                     ,
               "DFT RKS"                                       ,
               maximumEnergyAtoms      = 9                     ,
               maximumEnergyTests      = 20                    ,
               maximumGradientAtoms    = 5                     ,
               maximumGradientTests    = 0                     ,
               referenceDataFileName   = None, #"DFTRKSEnergies.yaml" ,
               testGradients           = True                  )
