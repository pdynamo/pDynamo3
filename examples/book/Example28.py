"""Example 28."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Define the MM and NB models.
mmModel = MMModelOPLS.WithParameterSet ( "bookSmallExamples" )
nbModel = NBModelMonteCarlo.WithDefaults ( )

# . Define the solute molecule.
solute = ImportSystem ( os.path.join ( molPath, "methane.mol" ) )
solute.Summary ( )

# . Define the solvent box.
solvent = ImportSystem ( os.path.join ( scratchPath, "water216_cubicBox_mc.pkl" ) )

# . Create the solvated system.
solution       = SolvateSystemBySuperposition ( solute, solvent )
solution.label = "Methane in Water."
solution.DefineMMModel ( mmModel )
solution.DefineNBModel ( nbModel )
solution.Summary ( )

# . Do a Monte Carlo calculation to equilibrate the system.
MonteCarlo_SystemGeometry ( solution        ,
                            blocks =     10 ,
                            moves  = 100000 )

# . Find out how many waters are left.
nWaters = len ( solution.connectivity.isolates ) - 1

# . Save the system.
ExportSystem ( os.path.join ( scratchPath, "ch4_water{:d}_cubicBox_mc.pkl".format ( nWaters ) ), solution )

# . Footer.
logFile.Footer ( )
