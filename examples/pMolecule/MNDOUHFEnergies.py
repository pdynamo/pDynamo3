"""MNDO UHF tests."""

from pMolecule.QCModel import QCModelMNDO
from QCTestSystems     import GetRadicalMoleculeSystems , \
                              RunQCTestSet

#===================================================================================================================================
# . Script.
#===================================================================================================================================
testSystems = GetRadicalMoleculeSystems ( convergerKeywords = { "densityTolerance" : 1.0e-12, "maximumIterations" : 250 } ,
                                          qcModelClass      = QCModelMNDO ,
                                          qcModelOptions    = { "AM1" : { "hamiltonian" : "am1" } } )
RunQCTestSet ( testSystems                                      ,
               "MNDO UHF"                                       ,
               maximumEnergyAtoms      = 100                    ,
               maximumEnergyTests      = 1000                   ,
               maximumGradientAtoms    = 100                    ,
               maximumGradientTests    = 20                     ,
               referenceDataFileName   = "MNDOUHFEnergies.yaml" ,
               testGradients           = True                   )
