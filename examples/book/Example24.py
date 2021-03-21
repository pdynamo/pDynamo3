"""Example 24."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Read the molecule definition.
molecule = ImportSystem ( os.path.join ( scratchPath, "bala_example23.pkl" ) )

# . Get the list of trajectory file names.
fileNames = glob.glob ( os.path.join ( scratchPath, "bala_window*.ptRes" ) )
fileNames.sort ( )

# . Calculate the PMF.
state = WHAM_ConjugateGradientMinimize ( fileNames                      ,
                                         bins                 = [ 100 ] ,
                                         logFrequency         =      1  ,
                                         maximumIterations    =   1000  ,
                                         rmsGradientTolerance = 1.0e-3  ,
                                         temperature          = 300.0   )

# . Write the PMF to a file.
histogram = state["Histogram"]
pmf       = state["PMF"      ]
histogram.ToTextFileWithData ( os.path.join ( scratchPath, "bala_pmf.dat" ), [ pmf ], format = "{:20.3f} {:20.3f}\n" )

# . Footer.
logFile.Footer ( )
