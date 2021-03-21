"""Test centering for the cutoff NB model."""

import glob, math, os, os.path

from Definitions               import dataPath                        , \
                                      outPath
from pBabel                    import ExportTrajectory                , \
                                      ImportSystem                    , \
                                      ImportTrajectory
from pCore                     import Clone                           , \
                                      logFile                         , \
                                      TestScriptExit_Fail
from pMolecule.MMModel         import MMModelOPLS
from pMolecule.NBModel         import NBModelCutOff
from pScientific.RandomNumbers import NormalDeviateGenerator          , \
                                      RandomNumberGenerator
from pSimulation               import LangevinDynamics_SystemGeometry

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Options.
_Destination = "NBModelCutOffCentering"
_NLog        =  100
_NSave       =  100
_NSteps      = 1000
_Tolerance   = 0.001

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Paths.
dataPath = os.path.join ( dataPath, "mol2"       )
outPath  = os.path.join ( outPath , _Destination )
if not os.path.exists ( outPath ): os.mkdir ( outPath )

# . NB models.
nbModelN = NBModelCutOff.WithDefaults ( ) ; nbModelN.useCentering = False
nbModelC = NBModelCutOff.WithDefaults ( ) ; nbModelC.useCentering = True

# . Set up the system.
system = ImportSystem ( os.path.join ( dataPath, "waterBox.mol2" ) )
system.DefineMMModel ( MMModelOPLS.WithParameterSet ( "bookSmallExamples" ) )
system.DefineNBModel ( Clone ( nbModelC ) )
system.Summary ( )
system.Energy  ( )

# . Do a short dynamics.
# . Define a normal deviate generator in a given state.
normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( RandomNumberGenerator.WithSeed ( 614108 ) )

# . Dynamics.
trajectoryPath = os.path.join ( outPath, "waterBox_C_1ps.dcd" )
trajectory     = ExportTrajectory ( trajectoryPath, system )
LangevinDynamics_SystemGeometry ( system                            ,
                                  collisionFrequency     =     25.0 ,
                                  logFrequency           =    _NLog ,
                                  normalDeviateGenerator = normalDeviateGenerator ,
                                  steps                  =  _NSteps ,
                                  temperature            =    300.0 ,
                                  timeStep               =    0.001 ,
                                  trajectories = [ ( trajectory, _NSave ) ] )

# . Calculate trajectory energies with different NB models.
energies = []
for nbModel in ( nbModelC, nbModelN ):
    system.DefineNBModel ( nbModel )
    system.nbModel.Summary ( )
    trajectory = ImportTrajectory ( trajectoryPath, system )
    trajectory.ReadHeader ( )
    e = []
    while trajectory.RestoreOwnerData ( ):
        e.append ( system.Energy ( log = None ) )
    trajectory.Close ( )
    energies.append ( e )
    system.nbModel.StatisticsSummary ( system )

# . Check deviations.
( e0, e1 ) = energies
maximumDeviation = 0.0
for i in range ( len ( e0 ) ):
    maximumDeviation = max ( maximumDeviation, math.fabs ( e0[i] - e1[i] ) )

# . Summary of results.
logFile.Paragraph ( "Maximum deviation = {:.5f}".format ( maximumDeviation ) )

# . Footer.
logFile.Footer ( )
if maximumDeviation > _Tolerance: TestScriptExit_Fail ( )
