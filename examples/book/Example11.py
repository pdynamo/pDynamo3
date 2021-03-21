"""Example 11."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Define the molecule and its energy model.
molecule = ImportSystem ( os.path.join ( molPath, "cyclohexane_chair.mol" ) )
molecule.DefineMMModel ( MMModelOPLS.WithParameterSet ( "bookSmallExamples" ) )
molecule.DefineNBModel ( NBModelFull.WithDefaults ( ) )
molecule.Summary ( )

# . Determine the starting energy.
eStart = molecule.Energy ( )

# . Optimization.
BakerSaddleOptimize_SystemGeometry ( molecule                      ,
                                     logFrequency         =      1 ,
                                     maximumIterations    =    250 ,
                                     rmsGradientTolerance = 1.0e-3 )

# . Determine the final energy.
eStop = molecule.Energy ( )

# . Print the energy change.
logFile.Paragraph ( "Energy change after search = {:20.4f}".format ( eStop - eStart ) )

# . Save the coordinates.
molecule.label = "Cyclohexane saddle conformation - OPLS-AA optimized."
ExportSystem ( os.path.join ( scratchPath, "cyclohexane_saddle.xyz" ), molecule )

# . Footer.
logFile.Footer ( )
