"""DFT RKS tests."""

from pMolecule.QCModel import QCModelDFT
from QCTestSystems     import GetClosedShellMoleculeSystems , \
                              RunQCTestSet

#===================================================================================================================================
# . Script.
#===================================================================================================================================
_QCModelOptions =  { "LDA:DZVP"   : { "fitBasis" : "dgauss-a1-dftjfit" , "functional" : "lda" , "orbitalBasis" : "dgauss-dzvp" } ,
                     "BLYP:DZVP"  : { "fitBasis" : "dgauss-a1-dftjfit" , "functional" : "blyp", "orbitalBasis" : "dgauss-dzvp" } ,
                     "LDA:SV(P)"  : { "fitBasis" : "def2-sv(p)-rifit"  , "functional" : "lda" , "orbitalBasis" : "def2-sv(p)"  } ,
                     "BLYP:SV(P)" : { "fitBasis" : "def2-sv(p)-rifit"  , "functional" : "blyp", "orbitalBasis" : "def2-sv(p)"  } }
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
