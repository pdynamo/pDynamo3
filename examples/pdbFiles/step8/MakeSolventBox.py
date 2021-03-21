"""Make an orthorhombic solvent box of the appropriate size."""

import os, os.path, sys
sys.path.append ( os.path.join ( os.path.dirname ( os.path.realpath ( __file__ ) ), "..", "data" ) )

from Definitions               import dataPath                                 , \
                                      outPath
from pBabel                    import ImportSystem
from pCore                     import logFile                                  , \
                                      Pickle                                   , \
                                      Unpickle
from pMolecule.MMModel         import MMModelOPLS
from pMolecule.NBModel         import NBModelCutOff
from pScientific.Symmetry      import CrystalSystemOrthorhombic                , \
                                      SymmetryParameters
from pScientific.RandomNumbers import NormalDeviateGenerator                   , \
                                      RandomNumberGenerator
from pSimulation               import BuildSolventBox                          , \
                                      ConjugateGradientMinimize_SystemGeometry , \
                                      LangevinDynamics_SystemGeometry          , \
                                      SolventMoleculeNumber                    , \
                                      SystemDensity                            , \
                                      SystemExtents

# . Header.
logFile.Header ( )

# . Parameters.
_Density =  1000.0 # . Density of water (kg m^-3).
_Refine  =  True
_Steps   = 10000

# . Box sizes.
_XBox = 75.0
_YBox = 60.0
_ZBox = 60.0

# . Define the solvent MM and NB models.
mmModel = MMModelOPLS.WithParameterSet ( "bookSmallExamples" )
nbModel = NBModelCutOff.WithDefaults ( )

# . Define the solvent molecule.
molecule = ImportSystem ( os.path.join ( dataPath, "water.mol" ) )
molecule.Summary ( )

# . Create a symmetry parameters instance with the correct dimensions.
symmetryParameters = SymmetryParameters ( )
symmetryParameters.SetCrystalParameters ( _XBox, _YBox, _ZBox, 90.0, 90.0, 90.0 )

# . Create the basic solvent box.
solvent = BuildSolventBox ( CrystalSystemOrthorhombic ( ), symmetryParameters, molecule, _Density )
solvent.label = "Water Box"
solvent.DefineMMModel ( mmModel )
solvent.DefineNBModel ( nbModel )
solvent.Summary ( )
solvent.Energy  ( )

# . Refine the system using minimization and dynamics.
if _Refine:

    # . Minimization.
    ConjugateGradientMinimize_SystemGeometry ( solvent                     ,
                                               maximumIterations    =  200 ,
                                               logFrequency         =    1 ,
                                               rmsGradientTolerance = 10.0 )

    # . Define a normal deviate generator in a given state.
    normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( RandomNumberGenerator.WithSeed ( 714717 ) )

    # . Dynamics.
    LangevinDynamics_SystemGeometry ( solvent                         ,
                                      collisionFrequency     =   25.0 ,
                                      logFrequency           =    100 ,
                                      normalDeviateGenerator = normalDeviateGenerator ,
                                      steps                  = _Steps ,
                                      temperature            =  300.0 ,
                                      timeStep               =  0.001 )
    # . Final energy.
    solvent.Energy  ( )

# . Calculate and print the final density.
logFile.Paragraph ( "Solvent density = {:.2f} kg m^-3.".format ( SystemDensity ( solvent ) ) )

# . Save the water box.
Pickle ( os.path.join ( outPath, "waterBox.pkl" ), solvent )

# . Footer.
logFile.Footer ( )
