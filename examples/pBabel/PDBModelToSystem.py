"""Read PDB model files and convert them to systems using atom data from original PDB files."""

import glob, os, os.path

from Definitions               import dataPath            , \
                                      outPath
from pBabel                    import ExportSystem        , \
                                      PDBFileReader       , \
                                      PDBModel
from pCore                     import logFile             , \
                                      TestScriptExit_Fail
from pMolecule.MMModel         import MMModelOPLS
from pMolecule.NBModel         import NBModelCutOff
from pScientific.RandomNumbers import RandomNumberGenerator
from pSimulation               import BuildHydrogenCoordinates3FromConnectivity

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Initialization.
numberErrors = 0

# . Paths.
modelPath = os.path.join ( dataPath, "pdbModel" )
pdbPath   = os.path.join ( dataPath, "pdb"      )
outPath   = os.path.join ( outPath , "xyz"      )

# . Set up the output directory.
if not os.path.exists ( outPath ): os.mkdir ( outPath )
outFiles = glob.glob ( os.path.join ( outPath, "*.xyz" ) )
for outFile in outFiles: os.remove ( outFile )

# . Get the files to process.
modelFiles = glob.glob ( os.path.join ( modelPath, "*.model" ) )

# . Get the model file names.
pdbNames = set ( )
for modelFile in modelFiles:
    ( head, tail ) = os.path.split ( modelFile )
    pdbNames.add ( tail[0:-6] )
pdbNames = list ( pdbNames )
pdbNames.sort ( )

# . Loop over the files.
for pdbName in pdbNames:

    # . Check file names.
    modelFile = os.path.join ( modelPath, pdbName + ".model" )
    pdbFile   = os.path.join ( pdbPath,   pdbName + ".pdb"   )
    if os.path.exists ( modelFile ) and os.path.exists ( pdbFile ):

        # . Get the model and its raw counterpart.
        model1   = PDBModel.FromModelFile       ( modelFile )
        rawModel = PDBFileReader.PathToPDBModel ( pdbFile   )

        # . Route 1 - make an atomic model.
        model1.Summary ( )
        try:

            # . Make the atomic model.
            model1.MakeAtomicModelFromComponentLibrary ( )
            model1.ExtractAtomData ( rawModel )
            model1.Summary ( )

            # . Make a system.
            system1 = model1.MakeSystem ( )
            system1.Summary ( )

            # . Add energy models.
            system1.DefineMMModel ( MMModelOPLS.WithParameterSet ( "protein" ) )
            system1.DefineNBModel ( NBModelCutOff.WithDefaults ( ) )

            # . Build as many undefined coordinates as possible.
            if system1.coordinates3.numberUndefined > 0:
                rng = RandomNumberGenerator.WithSeed ( 117513 )
                BuildHydrogenCoordinates3FromConnectivity ( system1, randomNumberGenerator = rng )

            # . Calculate an energy if all coordinates have been defined.
            if system1.coordinates3.numberUndefined <= 0:
                system1.Energy ( doGradients = True )

        # . Error.
        except Exception as e:
            numberErrors += 1
            logFile.Paragraph ( "Error occurred> " +  e.args[0] )

        # . Route 2 - extract atoms.
        model2  = PDBModel.FromModelFile ( modelFile )
        model2.ExtractAtoms ( rawModel )
        model2.Summary ( )
        system2 = model2.MakeSystem ( )
        system2.Summary ( )

        # . Output the xyz file if there are no undefined coordinates.
        n = system2.coordinates3.numberUndefined
        if n > 0 : logFile.Paragraph ( "System has {:d} undefined coordinates.".format ( n ) )
        else     : ExportSystem ( os.path.join ( outPath, pdbName + ".xyz" ), system2 )

        # . Separator.
        logFile.Separator ( )

# . Footer.
logFile.Footer ( )
if numberErrors > 0: TestScriptExit_Fail ( )

