"""Example 5."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Define the energy models.
energyModels = [ QCModelMNDO.WithOptions ( hamiltonian = "am1"  ) ,
                 QCModelMNDO.WithOptions ( hamiltonian = "mndo" ) ,
                 QCModelMNDO.WithOptions ( hamiltonian = "pm3"  ) ]

# . Get the filename.
fileName = os.path.join ( xyzPath, "water.xyz" )

# . Loop over the energy models.
results = []
for model in energyModels:
    molecule = ImportSystem ( fileName )
    molecule.DefineQCModel ( model )
    molecule.Summary ( )
    energy  = molecule.Energy        ( )
    charges = molecule.AtomicCharges ( )
    dipole  = molecule.DipoleMoment  ( )
    results.append ( ( model.hamiltonian.upper ( ), energy, charges, dipole.Norm2 ( ) ) )

# . Output the results.
table = logFile.GetTable ( columns = [ 10, 20, 20, 20, 20, 20 ] )
table.Start  ( )
table.Title  ( "Energy Model Results for Water" )
table.Heading ( "Model"  )
table.Heading ( "Energy" )
table.Heading ( "Charges", columnSpan = 3 )
table.Heading ( "Dipole" )
for ( label, energy, charges, dipole ) in results:
    table.Entry ( label, align = Align.Left )
    table.Entry ( "{:.1f}".format ( energy ) )
    for charge in charges: table.Entry ( "{:.3f}".format ( charge ) )
    table.Entry ( "{:.3f}".format ( dipole ) )
table.Stop ( )

# . Footer.
logFile.Footer ( )
