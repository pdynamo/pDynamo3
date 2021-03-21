"""Test for vibrational dynamics."""

import glob, math, os

from Definitions               import dataPathM                       , \
                                      outPath
from pBabel                    import ExportTrajectory                , \
                                      ImportSystem                    , \
                                      ImportTrajectory
from pCore                     import Clone                           , \
                                      logFile                         , \
                                      Selection                       , \
                                      TestScriptExit_Fail
from pMolecule.MMModel         import MMModelOPLS
from pMolecule.NBModel         import NBModelFull
from pMolecule.QCModel         import DIISSCFConverger                , \
                                      QCModelMNDO
from pScientific.Arrays        import Array
from pScientific.RandomNumbers import NormalDeviateGenerator          , \
                                      RandomNumberGenerator
from pSimulation               import LangevinDynamics_SystemGeometry , \
                                      LBFGSMinimize_SystemGeometry    , \
                                      NormalModes_SystemGeometry      , \
                                      QuasiHarmonic_SystemGeometry

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Dynamics options.
_CollisionFrequency =     1.0
_LogFrequency       =   10000
_NSteps0            =    1000
_NSteps1            =  100000
_Seed               =  156451
_SaveFrequency      =     100
_Temperature        =  1000.0

# . Molecules.
_MoleculeLabels = ( "hydrogenFluoride", "water", "benzene", "cyclohexane" )
_QCModels       = ( "hydrogenFluoride", )

# . Optimization options.
_Iterations   = 1000
_Tolerance    = 1.0e-05

# . Paths.
_Destination = "vibrationalDynamics"

# . Tolerances.
_RMSAbsoluteErrorTolerance = 1.0e-03

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Paths.
dataPath = os.path.join ( dataPathM, "mol"        )
outPath  = os.path.join ( outPath  , _Destination )
if not os.path.exists ( outPath ): os.mkdir ( outPath )

# . Initialization.
converger      = DIISSCFConverger.WithOptions ( densityTolerance = 1.0e-10 )
numberFailures = 0

# . Loop over molecules.
for moleculeLabel in _MoleculeLabels:

    # . Get the system.
    system = ImportSystem ( os.path.join ( dataPath, moleculeLabel + ".mol" ) )
    if moleculeLabel in _QCModels:
        system.DefineQCModel ( QCModelMNDO.WithOptions ( converger = converger ) )
    else:
        system.DefineMMModel ( MMModelOPLS.WithParameterSet ( "protein" ) )
        system.DefineNBModel ( NBModelFull.WithDefaults ( )               )
    system.Summary ( )
    system.Energy  ( )

    # . Minimize well.
    LBFGSMinimize_SystemGeometry ( system,
                                   logFrequency         = _LogFrequency ,
                                   maximumIterations    =   _Iterations ,
                                   rmsGradientTolerance =    _Tolerance )

    # . Normal mode analysis.
    nmState = NormalModes_SystemGeometry ( system )

    # . Do a dynamics simulation: equilibration and then data collection.
    normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( RandomNumberGenerator.WithSeed ( _Seed ) )
    LangevinDynamics_SystemGeometry ( system                                          ,
                                      collisionFrequency     = _CollisionFrequency    ,
                                      logFrequency           = _LogFrequency          ,
                                      normalDeviateGenerator = normalDeviateGenerator ,
                                      steps                  = _NSteps0               ,
                                      temperature            = _Temperature           ,
                                      timeStep               = 0.001                  )
    reference3 = Clone ( system.coordinates3 )
    trajectory = ExportTrajectory ( os.path.join ( outPath, moleculeLabel + ".mdcrd" ), system )
    LangevinDynamics_SystemGeometry ( system                                          ,
                                      collisionFrequency     = _CollisionFrequency    ,
                                      logFrequency           = _LogFrequency          ,
                                      normalDeviateGenerator = normalDeviateGenerator ,
                                      steps                  = _NSteps1               ,
                                      temperature            = _Temperature           ,
                                      timeStep               = 0.001                  ,
                                      trajectories           = [ ( trajectory, _SaveFrequency ) ] )

    # . Check RMSs.
    masses = Array.FromIterable ( [ atom.mass for atom in system.atoms ] )
    rms0   = system.coordinates3.RootMeanSquareDeviation ( reference3, weights = masses )
    system.coordinates3.Superimpose ( reference3, weights = masses )
    rms1 = system.coordinates3.RootMeanSquareDeviation ( reference3, weights = masses )
    if ( math.fabs ( rms1 - rms0 ) >= _RMSAbsoluteErrorTolerance ): numberFailures += 1

    # . Do a quasi-harmonic analysis.
    qhState = QuasiHarmonic_SystemGeometry ( system, temperature = _Temperature, trajectories = [ os.path.join ( outPath, moleculeLabel + ".mdcrd" ) ] )

# . Footer.
logFile.Footer ( )
if ( numberFailures != 0 ): TestScriptExit_Fail ( )
