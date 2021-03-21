"""Some definitions."""

import glob, math, os.path

from pBabel                    import ExportSystem                              , \
                                      ExportTrajectory                          , \
                                      ImportCoordinates3                        , \
                                      ImportSystem                              , \
                                      ImportTrajectory
from pCore                     import Clone                                     , \
                                      logFile                                   , \
                                      Pickle                                    , \
                                      Selection                                 , \
                                      TestScript_InputDataPath                  , \
                                      TestScript_OutputDataPath                 , \
                                      Unpickle
from pMolecule                 import AtomSelection                             , \
                                      RestraintDistance                         , \
                                      RestraintEnergyModel                      , \
                                      RestraintModel                            , \
                                      RestraintMultipleTether
from pMolecule.MMModel         import MMModelOPLS
from pMolecule.NBModel         import NBModelCutOff
from pScientific.Arrays        import Array
from pScientific.RandomNumbers import NormalDeviateGenerator                    , \
                                      RandomNumberGenerator
from pScientific.Statistics    import Statistics
from pScientific.Symmetry      import CrystalSystemCubic                        , \
                                      SymmetryParameters
from pSimulation               import AddCounterIons                            , \
                                      BuildHydrogenCoordinates3FromConnectivity , \
                                      BuildSolventBox                           , \
                                      ConjugateGradientMinimize_SystemGeometry  , \
                                      DetermineSolvationParameters              , \
                                      LangevinDynamics_SystemGeometry           , \
                                      SolvateSystemBySuperposition              , \
                                      SystemDensity

# . Local name.
_name = "aSmallProtein"

# . The input and output data paths.
dataPath = TestScript_InputDataPath  ( _name )
outPath  = TestScript_OutputDataPath ( _name )
