"""Example 17."""

# . Note that a simple average will not work if the phi or psi angles cross the +/-180 boundary.

from Definitions import *

# . Header.
logFile.Header ( )

# . Read the molecule definition.
molecule = ImportSystem ( os.path.join ( molPath, "bala_c7eq.mol" ) )
molecule.Summary ( )

# . Define the trajectory.
trajectory = ImportTrajectory ( os.path.join ( scratchPath, "bala_c7eq.ptGeo" ), molecule )
trajectory.ReadHeader ( )

# . Loop over the frames in the trajectory.
phi = []
psi = []
while trajectory.RestoreOwnerData ( ):
    phi.append ( molecule.coordinates3.Dihedral ( 4, 6,  8, 14 ) )
    psi.append ( molecule.coordinates3.Dihedral ( 6, 8, 14, 16 ) )

# . Finish up.
trajectory.ReadFooter ( )
trajectory.Close ( )

# . Set up the statistics calculation.
phiStatistics = Statistics ( phi )
psiStatistics = Statistics ( psi )

# . Output the results.
table = logFile.GetTable ( columns = [ 20, 20, 20 ] )
table.Start   ( )
table.Title   ( "Phi/Psi Angles" )
table.Heading ( "Frame" )
table.Heading ( "Phi"   )
table.Heading ( "Psi"   )
for ( i, ( h, s ) ) in enumerate ( zip ( phi, psi ) ):
    table.Entry ( "{:d}"  .format ( i ) )
    table.Entry ( "{:.2f}".format ( h ) )
    table.Entry ( "{:.2f}".format ( s ) )
table.Entry ( "Mean:",               align = Align.Left )
table.Entry ( "{:.2f}".format ( phiStatistics.mean ) )
table.Entry ( "{:.2f}".format ( psiStatistics.mean ) )
table.Entry ( "Standard Deviation:", align = Align.Left )
table.Entry ( "{:.2f}".format ( phiStatistics.standardDeviation ) )
table.Entry ( "{:.2f}".format ( psiStatistics.standardDeviation ) )
table.Stop ( )

# . Footer.
logFile.Footer ( )
