"""Example 13."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Define the molecule and its energy model.
molecule = ImportSystem ( os.path.join ( molPath, "cyclohexane_chair.mol" ) )
molecule.DefineMMModel ( MMModelOPLS.WithParameterSet ( "bookSmallExamples" ) )
molecule.DefineNBModel ( NBModelFull.WithDefaults ( ) )
molecule.Summary ( )

# . Assign the reactant and product coordinates.
reactants = ImportCoordinates3 ( os.path.join ( xyzPath, "cyclohexane_chair.xyz"     ) )
products  = ImportCoordinates3 ( os.path.join ( xyzPath, "cyclohexane_twistboat.xyz" ) )

# . Get an initial path.
trajectoryPath = os.path.join ( scratchPath, "cyclohexane_sawpath.ptGeo" )
GrowingStringInitialPath ( molecule, 11, reactants, products, trajectoryPath )

# . Optimization.
trajectory = ExportTrajectory ( trajectoryPath, molecule, append = True )
ChainOfStatesOptimizePath_SystemGeometry ( molecule                   ,
                                           trajectory                 ,
                                           logFrequency         =   1 ,
                                           maximumIterations    = 100 ,
                                           rmsGradientTolerance = 0.1 )

# . Footer.
logFile.Footer ( )
