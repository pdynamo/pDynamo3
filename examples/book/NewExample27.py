"""Example 27.

Updated version using different box builder and Monte Carlo refinement.
"""

from Definitions import *

# . Header.
logFile.Header ( )

# . Define some parameters.
_Density      = 1000.0 # . Density of water (kg m^-3).
_MoleculeName = "water"
_Molecules    = 216

# . Define the solvent MM and NB models.
mmModel = MMModelOPLS.WithParameterSet ( "bookSmallExamples" )
nbModel = NBModelMonteCarlo.WithDefaults ( )

# . Define the solvent molecule.
molecule = ImportSystem ( os.path.join ( molPath, _MoleculeName + ".mol" ) )
molecule.Summary ( )

# . Get the dimensions of a cubic box with the required number of molecules.
a = SolventCubicBoxDimensions ( molecule, _Molecules, _Density )

# . Create a symmetryParameters instance with the correct dimensions.
symmetryParameters = SymmetryParameters ( )
symmetryParameters.SetCrystalParameters ( a, a, a, 90.0, 90.0, 90.0 )

# . Create the basic solvent box.
solvent = BuildSolventBox ( CrystalSystemCubic ( ), symmetryParameters, molecule, _Density )
solvent.label = "Water Box"
solvent.DefineMMModel ( mmModel )
solvent.DefineNBModel ( nbModel )
solvent.Summary ( )

# . Do Monte Carlo simulations to equilibrate the system.
MonteCarlo_SystemGeometry ( solvent                  ,
                            blocks          =     50 ,
                            log             =   None ,
                            moves           =   1000 ,
                            volumeFrequency =      0 )
MonteCarlo_SystemGeometry ( solvent                  ,
                            blocks          =     20 ,
                            moves           = 100000 )

# . Calculate and print the final density.
logFile.Paragraph ( "Solvent density = {:.2f} kg m^-3.".format ( SystemDensity ( solvent ) ) )

# . Save the system.
ExportSystem ( os.path.join ( scratchPath, "{:s}{:d}_cubicBox_mc.pkl".format ( _MoleculeName, _Molecules ) ), solvent )

# . Footer.
logFile.Footer ( )
