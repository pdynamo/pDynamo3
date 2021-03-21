"""MNDO RHF tests."""

from pMolecule.QCModel import QCModelMNDO
from QCTestSystems     import GetClosedShellMoleculeSystems , \
                              RunQCTestSet

#===================================================================================================================================
# . Script.
#===================================================================================================================================
testSystems = GetClosedShellMoleculeSystems ( convergerKeywords = { "densityTolerance" : 1.0e-12, "maximumIterations" : 250 } ,
                                              qcModelClass      = QCModelMNDO ,
                                              qcModelOptions    = { "AM1"      : { "hamiltonian" : "am1"      } ,
                                                                    "MNDO"     : { "hamiltonian" : "mndo"     } ,
                                                                    "PDDGMNDO" : { "hamiltonian" : "pddgmndo" } ,
                                                                    "PDDGPM3"  : { "hamiltonian" : "pddgpm3"  } ,
                                                                    "PM3"      : { "hamiltonian" : "pm3"      } ,
                                                                    "PM6"      : { "hamiltonian" : "pm6"      } ,
                                                                    "RM1"      : { "hamiltonian" : "rm1"      } } )
RunQCTestSet ( testSystems                                      ,
               "MNDO RHF"                                       ,
               maximumEnergyAtoms      = 30                     ,
               maximumEnergyTests      = 1000                   ,
               maximumGradientAtoms    = 10                     ,
               maximumGradientTests    = 20                     ,
               referenceDataFileName   = "MNDORHFEnergies.yaml" ,
               testGradients           = True                   )
