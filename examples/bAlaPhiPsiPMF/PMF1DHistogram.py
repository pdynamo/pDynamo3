"""1-D PMFs - histogram the generated data."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Loop over angles.
for ( tag, indices ) in ( ( "Phi", phiAtomIndices ), ( "Psi", psiAtomIndices ) ):

    tagL = tag.lower ( )

    # . Get the system.
    system = Unpickle ( os.path.join ( outPath, "bAla.pkl" ) )
    system.Summary ( )

    # . Get the list of trajectory file names.
    fileNames = [ fileName for fileName in glob.glob ( os.path.join ( outPath, "bAla_" + tagL + "_*.ptRes" ) ) if ( fileName.count ( "_" ) == 2 ) ]
    fileNames.sort ( )

    # . Histogram the trajectory data.
    handler   = SystemRestraintTrajectoryDataHandler.FromTrajectoryPaths ( fileNames )
    histogram = handler.HistogramData ( [ 360 ] )
    counts    = [ float ( count ) for count    in histogram.counts ]
    histogram.ToTextFileWithData ( os.path.join ( outPath, tagL + "Values.dat" ), [ counts ], format = "{:20.3f} {:20.3f}\n" )

# . Footer.
logFile.Footer ( )
