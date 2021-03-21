"""Example 19."""

from Definitions import *

# . Header.
logFile.Header ( )

# . Define the cut-offs.
_CutOffs = [ ( float ( i ), float ( i ) + 4.0 ) for i in range ( 1, 41 ) ]

# . Set up the system.
molecule = ImportSystem ( os.path.join ( pdbPath, "crambin.pdb" ), useComponentLibrary = True )
molecule.DefineMMModel ( MMModelOPLS.WithParameterSet ( "protein" ) )
molecule.DefineNBModel ( NBModelFull.WithDefaults ( ) )
molecule.Summary ( )

# . Build missing coordinates.
rng = RandomNumberGenerator ( )
rng.SetSeed ( 861199 )
BuildHydrogenCoordinates3FromConnectivity ( molecule, randomNumberGenerator = rng )
# . Stop if not all coordinates are defined.
if molecule.coordinates3.numberUndefined > 0:
    logFile.Paragraph ( "The system has undefined coordinates." )
# . Continue.
else:
    # . Get the energy with no cut-off.
    eF = molecule.Energy ( log = None )

    # . Get the energies with different cut-off values.
    table = logFile.GetTable ( columns = [ 20, 20 ] )
    table.Start   ( )
    table.Title   ( "Cut-off/Full Energy Difference" )
    table.Heading ( "Inner Cut-off" )
    table.Heading ( "Difference"    )
    for ( cut, cutB ) in _CutOffs:
        generator           = PairListGenerator.WithOptions       ( cutOff      = cutB )
        pairwiseInteraction = PairwiseInteractionABFS.WithOptions ( innerCutOff = cut, outerCutOff = cutB )
        nbModel = NBModelCutOff.WithOptions ( generator = generator, pairwiseInteraction = pairwiseInteraction )
        molecule.DefineNBModel ( nbModel )
        eT = molecule.Energy ( log = None )
        table.Entry ( "{:.1f}".format ( cut     ) )
        table.Entry ( "{:.4f}".format ( eT - eF ) )
    table.Stop ( )

# . Footer.
logFile.Footer ( )
