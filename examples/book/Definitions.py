"""Definitions needed by the examples."""

import glob, math, os, os.path

from pBabel                    import ExportSystem                                 , \
                                      ExportTrajectory                             , \
                                      ImportCoordinates3                           , \
                                      ImportSystem                                 , \
                                      ImportTrajectory                             , \
                                      SMILESReader
from pCore                     import Align                                        , \
                                      Clone                                        , \
                                      logFile                                      , \
                                      Selection                                    , \
                                      TestScript_InputDataPath                     , \
                                      TestScript_OutputDataPath                    , \
                                      XHTMLLogFileWriter
from pMolecule                 import RestraintDistance                            , \
                                      RestraintEnergyModel                         , \
                                      RestraintModel                               , \
                                      RestraintMultipleTether                      , \
                                      System                                       , \
                                      SystemGeometryObjectiveFunction
from pMolecule.MMModel         import MMModelOPLS
from pMolecule.NBModel         import NBModelCutOff                                , \
                                      NBModelFull                                  , \
                                      NBModelMonteCarlo                            , \
                                      NBModelORCA                                  , \
                                      PairwiseInteractionABFS
from pMolecule.QCModel         import QCModelMNDO
from pScientific               import Constants                                    , \
                                      Units
from pScientific.Arrays        import Array                                        , \
                                      ArrayPrint2D
from pScientific.Geometry3     import Coordinates3                                 , \
                                      PairListGenerator
from pScientific.RandomNumbers import NormalDeviateGenerator                       , \
                                      RandomNumberGenerator
from pScientific.Statistics    import Statistics
from pScientific.Symmetry      import CrystalSystemCubic                           , \
                                      PeriodicBoundaryConditions                   , \
                                      SymmetryParameters
from pSimulation               import BakerSaddleOptimize_SystemGeometry           , \
                                      BuildCubicSolventBox                         , \
                                      BuildHydrogenCoordinates3FromConnectivity    , \
                                      BuildSolventBox                              , \
                                      ChainOfStatesOptimizePath_SystemGeometry     , \
                                      ConjugateGradientMinimize_SystemGeometry     , \
                                      GrowingStringInitialPath                     , \
                                      LeapFrogDynamics_SystemGeometry              , \
                                      MergeByAtom                                  , \
                                      ModifyOption                                 , \
                                      MonteCarlo_IsolateInteractionEnergy          , \
                                      MonteCarlo_ScaleIsolateInteractionParameters , \
                                      MonteCarlo_SystemGeometry                    , \
                                      NormalModes_SystemGeometry                   , \
                                      NormalModesTrajectory_SystemGeometry         , \
                                      PruneByAtom                                  , \
                                      RadialDistributionFunction                   , \
                                      SelfDiffusionFunction                        , \
                                      SolventCubicBoxDimensions                    , \
                                      SolvateSystemBySuperposition                 , \
                                      SteepestDescentPath_SystemGeometry           , \
                                      SystemDensity                                , \
                                      ThermodynamicsRRHO_SystemGeometry            , \
                                      VelocityVerletDynamics_SystemGeometry        , \
                                      WHAM_ConjugateGradientMinimize

# . Local name.
_name  = "book"

# . The input data paths.
dataPath = TestScript_InputDataPath ( _name )
molPath  = os.path.join ( dataPath, "mol" )
pdbPath  = os.path.join ( dataPath, "pdb" )
xyzPath  = os.path.join ( dataPath, "xyz" )

# . Input and output data paths.
pklPath = os.path.join ( dataPath, "pkl" )
if not os.path.exists ( pklPath ): os.mkdir ( pklPath )

# . The output data path.
scratchPath = TestScript_OutputDataPath ( _name )
