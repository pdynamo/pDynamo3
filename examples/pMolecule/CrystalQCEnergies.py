"""Crystal tests - QC."""

from CrystalMMEnergies import RunCrystalTest
from pMolecule.NBModel import QCMMElectrostaticModelDensityCutOffMNDO , \
                              QCMMLennardJonesModelCutOff             , \
                              QCQCElectrostaticModelMultipoleCutOff   , \
                              QCQCLennardJonesModelCutOff
from pMolecule.QCModel import DIISSCFConverger                        , \
                              QCModelMNDO

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
_QCModel    = QCModelMNDO.WithOptions ( converger = DIISSCFConverger.WithOptions ( densityTolerance = 1.0e-12, maximumIterations = 250 ) )
_QCQCModels = { "qcmmElectrostatic" : QCMMElectrostaticModelDensityCutOffMNDO.WithDefaults ( ) ,
                "qcmmLennardJones"  : QCMMLennardJonesModelCutOff.WithDefaults             ( ) ,
                "qcqcElectrostatic" : QCQCElectrostaticModelMultipoleCutOff.WithDefaults   ( ) ,
                "qcqcLennardJones"  : QCQCLennardJonesModelCutOff.WithDefaults             ( ) }

#===================================================================================================================================
# . Script.
#===================================================================================================================================
RunCrystalTest ( checkEnergies    = False       ,
                 dataSetTag       = "QC"        ,
                 doQCMM           = False       ,
                 doQCQC           = True        ,
                 geometryOptimize = False       ,
                 qcModel          = _QCModel    ,
                 qcmmModels       = _QCQCModels )
