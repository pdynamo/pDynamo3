"""Example 27."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Define some parameters.
_Linear       =          6
_MoleculeName =    "water"
_Molecules    = _Linear**3

# . Define the MM and NB models.
mmModel = MMModelOPLS.WithParameterSet ( "bookSmallExamples" )
nbModel = NBModelMonteCarlo.WithDefaults ( )

# . Define the solvent molecule.
molecule = ImportSystem ( os.path.join ( molPath, _MoleculeName + ".mol" ) )
molecule.Summary ( )

# . Build the cubic system.
solvent = BuildCubicSolventBox ( molecule, _Molecules )
solvent.label = "Water Box"
solvent.DefineMMModel ( mmModel )
solvent.DefineNBModel ( nbModel )
solvent.Summary ( )

# . Do Monte Carlo simulations to equilibrate the system.
MonteCarlo_SystemGeometry ( solvent           ,
                            blocks   =     10 ,
                            moves    = 100000 ,
                            pressure = 1000.0 )
MonteCarlo_SystemGeometry ( solvent           ,
                            blocks   =     10 ,
                            moves    = 100000 )

# . Calculate and print the final density.
mass    = sum ( [ atom.mass for atom in solvent.atoms ] )
volume  = solvent.symmetryParameters.volume
density = ( mass / volume ) * ( Units.Mass_AMU_To_Kg * 1.0e+30 )
logFile.Paragraph ( "Solvent density = {:.2f} kg m^-3.".format ( density ) )

# . Save the system.
ExportSystem ( os.path.join ( scratchPath, "{:s}{:d}_cubicBox_mc.pkl".format ( _MoleculeName, _Molecules ) ), solvent )

# . Footer.
logFile.Footer ( )
