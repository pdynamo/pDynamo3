"""Example 3."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Read in a system.
molecule = ImportSystem ( os.path.join ( xyzPath, "bala_c7eq.xyz" ) )
molecule.Summary ( )

# . Define a table for the results.
table = logFile.GetTable ( columns = [ 15, 15 ] )
table.Start ( )
table.Title ( "Bond Analysis" )
table.Heading ( "Safety Factor" )
table.Heading ( "Bonds Found"   )

# . Loop over the buffer sizes.
for i in range ( 21 ):
    safety = 0.1 * float ( i )
    molecule.BondsFromCoordinates3 ( safety = safety )
    table.Entry ( "{:4.1f}".format ( safety ) )
    table.Entry ( "{:d}".format ( len ( molecule.connectivity.bonds ) ) )

# . Finish up the table.
table.Stop ( )

# . Generate the bonds with the default safety factor.
molecule.BondsFromCoordinates3 ( safety = 0.5 )
molecule.Summary ( )

# . Print the bonds.
table = logFile.GetTable ( columns = 4 * [ 5, 5, 10 ] )
table.Start ( )
table.Title ( "Bond Lengths (Angstroms)" )
for ( i, j ) in molecule.connectivity.bondIndices:
    table.Entry ( "{:d}".format ( i ) )
    table.Entry ( "{:d}".format ( j ) )
    table.Entry ( "{:6.3f}".format ( molecule.coordinates3.Distance ( i, j ) ) )
table.Stop ( )

# . Print the angles.
table = logFile.GetTable ( columns = 3 * [ 5, 5, 5, 10 ] )
table.Start ( )
table.Title ( "Angles (Degrees)" )
for ( i, j, k ) in molecule.connectivity.angleIndices:
    table.Entry ( "{:d}".format ( i ) )
    table.Entry ( "{:d}".format ( j ) )
    table.Entry ( "{:d}".format ( k ) )
    table.Entry ( "{:6.1f}".format ( molecule.coordinates3.Angle ( i, j, k ) ) )
table.Stop ( )

# . Print the dihedrals.
table = logFile.GetTable ( columns = 4 * [ 5, 5, 5, 5, 10 ] )
table.Start ( )
table.Title ( "Dihedrals (Degrees)" )
for ( i, j, k, l ) in molecule.connectivity.dihedralIndices:
    table.Entry ( "{:d}".format ( i ) )
    table.Entry ( "{:d}".format ( j ) )
    table.Entry ( "{:d}".format ( k ) )
    table.Entry ( "{:d}".format ( l ) )
    table.Entry ( "{:6.1f}".format ( molecule.coordinates3.Dihedral ( i, j, k, l ) ) )
table.Stop ( )

# . Footer.
logFile.Footer ( )
