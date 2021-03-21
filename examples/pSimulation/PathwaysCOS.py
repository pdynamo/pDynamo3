"""Pathways tests."""

import math, os, os.path

from Definitions      import outPath
from pBabel           import ExportTrajectory
from pCore            import CPUTime                                  , \
                             logFile                                  , \
                             TestScriptExit_Fail
from pSimulation      import ChainOfStatesOptimizePath_SystemGeometry , \
                             Duplicate                                , \
                             SGOFProcessPoolFactory
from PathTestSystems  import PathTestReportsSummary                   , \
                             PathTestSystems

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Options.
_DoLong           = False
_DoShort          = True
_NumberOfImages   = 11
_PoolFactory      = None # SGOFProcessPoolFactory ( maximumProcesses = 6, poolType = "Multiprocessing" )
_OptimizerOptions = { "fixedTerminalImages"                        : True         ,
                      "forceOneSingleImageOptimization"            : False        ,
                      "forceSingleImageOptimizations"              : False        ,
                      "forceSplineRedistributionCheckPerIteration" : False        ,
                      "freezeRMSGradientTolerance"                 :    0.0       ,
                      "logFrequency"                               :    1         ,
                      "maximumIterations"                          :  5000        ,
                      "optimizer"                                  : None         ,
                      "poolFactory"                                : _PoolFactory ,
                      "rmsGradientTolerance"                       :    0.1       ,
                      "splineRedistributionTolerance"              :    1.5       ,
                      "springForceConstant"                        :    0.0       ,
                      "useSplineRedistribution"                    : True         }

# . Pathways directory.
_pathways = "pathways"

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Initialization.
isOK    = True
outPath = os.path.join ( outPath, _pathways )
if not os.path.exists ( outPath ): os.mkdir ( outPath )

# . Loop over the tests.
reports = {}
for testSystem in PathTestSystems:

    # . Check whether to skip this system.
    if ( testSystem.isLong and ( not _DoLong ) ) or ( ( not testSystem.isLong ) and ( not _DoShort ) ): continue

    # . Get the system.
    system = testSystem.GetSystem ( )
    system.Summary ( )

    # . Create a starting trajectory.
    trajectoryPath = os.path.join ( outPath, "COS_" + testSystem.tag + ".ptGeo" )
    testSystem.GetGrowingStringInitialPath ( _NumberOfImages, trajectoryPath )

    # . Pathway.
    cpuTime        = CPUTime ( )
    options        = dict ( _OptimizerOptions )
    options["log"] = logFile
    trajectory     = ExportTrajectory ( trajectoryPath, system, append = True )
    report         = ChainOfStatesOptimizePath_SystemGeometry ( system, trajectory, **options )
    report["Time"] = cpuTime.Current ( )
    reports[testSystem.tag] = report

    # . Convert the trajectory.
    amberPath = trajectoryPath[0:-4] + ".mdcrd"
    Duplicate ( trajectoryPath, amberPath, system )

# . Footer.
numberConverged = PathTestReportsSummary ( reports )
logFile.Footer ( )
if numberConverged != len ( reports ): TestScriptExit_Fail ( )
