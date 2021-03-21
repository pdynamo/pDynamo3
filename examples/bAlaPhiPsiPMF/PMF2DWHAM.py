"""2-D PMFs - histogram the generated data."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Get the system.
system = Unpickle ( os.path.join ( outPath, "bAla.pkl" ) )
system.Summary ( )

# . Get the list of trajectory file names.
fileNames = glob.glob ( os.path.join ( outPath, "bAla_phi_*_psi_*.ptRes" ) )
fileNames.sort ( )

# . Calculate the PMF.
state = WHAM_ConjugateGradientMinimize ( fileNames                           ,
                                         bins                 = [ 100, 100 ] ,
                                         logFrequency         =          1   ,
                                         maximumIterations    =       1000   ,
                                         rmsGradientTolerance =     1.0e-3   ,
                                         temperature          =     300.0    )

# . Write the PMF to a file.
histogram = state["Histogram"]
pmf       = state["PMF"      ]
histogram.ToTextFileWithData ( os.path.join ( outPath, "phiPsiPMF.dat" ), [ pmf ], format = "{:20.3f} {:20.3f} {:20.3f}\n" )

# . Footer.
logFile.Footer ( )
