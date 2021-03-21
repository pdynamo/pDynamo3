"""Example 7."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Define the list of structures.
xyzFiles = [ "bala_alpha.xyz", "bala_c5.xyz", "bala_c7ax.xyz", "bala_c7eq.xyz" ]

# . Define the MM and NB models.
mmModel = MMModelOPLS.WithParameterSet ( "bookSmallExamples" )
nbModel = NBModelFull.WithDefaults ( )

# . Generate the molecule.
molecule = ImportSystem ( os.path.join ( molPath, "bala_c7eq.mol" ) )
molecule.DefineMMModel ( mmModel )
molecule.DefineNBModel ( nbModel )
molecule.Summary ( )

# . Loop over the structures in the xyz files.
results = []
for xyzFile in xyzFiles:
    molecule.coordinates3 = ImportCoordinates3 ( os.path.join ( xyzPath, xyzFile ) )
    energy = molecule.Energy       ( )
    dipole = molecule.DipoleMoment ( )
    results.append ( ( xyzFile[5:-4], energy, dipole.Norm2 ( ) ) )

# . Output the results.
table = logFile.GetTable ( columns = [ 20, 20, 20 ] )
table.Start  ( )
table.Title  ( "Energy Model Results for bALA" )
table.Heading ( "Conformation" )
table.Heading ( "Energy" )
table.Heading ( "Dipole" )
for ( label, energy, dipole ) in results:
    table.Entry ( label, align = Align.Left )
    table.Entry ( "{:.1f}".format ( energy ) )
    table.Entry ( "{:.3f}".format ( dipole ) )
table.Stop ( )

# . Footer.
logFile.Footer ( )
