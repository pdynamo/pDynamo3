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

# . Histogram the trajectory data.
handler   = SystemRestraintTrajectoryDataHandler.FromTrajectoryPaths ( fileNames )
histogram = handler.HistogramData ( [ 180, 180 ] )
counts    = [ float ( count ) for count    in histogram.counts ]
histogram.ToTextFileWithData ( os.path.join ( outPath, "phiPsiValues.dat" ), [ counts ], format = "{:20.3f} {:20.3f} {:20.3f}\n" )

# . Footer.
logFile.Footer ( )
