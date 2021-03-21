"""Some definitions."""

import glob, math, os.path

from pBabel                    import ExportSystem                             , \
                                      ExportTrajectory                         , \
                                      ImportCoordinates3                       , \
                                      ImportSystem                             , \
                                      ImportTrajectory                         , \
                                      SystemRestraintTrajectoryDataHandler
from pCore                     import Align                                    , \
                                      logFile                                  , \
                                      Pickle                                   , \
                                      TestScript_InputDataPath                 , \
                                      TestScript_OutputDataPath                , \
                                      Unpickle
from pMolecule                 import RestraintDihedral                        , \
                                      RestraintEnergyModel                     , \
                                      RestraintModel
from pMolecule.MMModel         import MMModelOPLS
from pMolecule.NBModel         import NBModelFull
from pScientific.RandomNumbers import NormalDeviateGenerator                   , \
                                      RandomNumberGenerator
from pScientific.Statistics    import RegularHistogram                         , \
                                      RegularHistogramDimension                , \
                                      Statistics
from pSimulation               import ConjugateGradientMinimize_SystemGeometry , \
                                      LangevinDynamics_SystemGeometry          , \
                                      WHAM_Bootstrapping                       , \
                                      WHAM_ConjugateGradientMinimize

# . Local name.
_name = "bAlaPhiPsiPMF"

# . The input and output data paths.
dataPath = TestScript_InputDataPath  ( _name )
outPath  = TestScript_OutputDataPath ( _name )

# . Phi/psi atom indices.
phiAtomIndices = ( 4, 6,  8, 14 )
psiAtomIndices = ( 6, 8, 14, 16 )
