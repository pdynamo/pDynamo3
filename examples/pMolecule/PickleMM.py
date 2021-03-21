"""Test pickling and unpickling."""

import glob, math, os, os.path

from collections       import defaultdict
from Definitions       import dataPath                   , \
                              outPath
from pBabel            import ImportSystem
from pCore             import logFile                    , \
                              Pickle                     , \
                              TestScriptExit_Fail        , \
                              Unpickle                   , \
                              YAMLMappingFile_FromObject , \
                              YAMLMappingFile_ToObject   , \
                              YAMLPickle                 , \
                              YAMLUnpickle
from pMolecule         import System
from pMolecule.MMModel import MMModelCHARMM
from pMolecule.NBModel import NBModelCutOff
from pSimulation       import BuildHydrogenCoordinates3FromConnectivity

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
_Destination = "picklingMM"
_Tolerance   = 0.1

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Paths.
dataPath = os.path.join ( dataPath, "pdb"        )
outPath  = os.path.join ( outPath , _Destination )
if not os.path.exists ( outPath ): os.mkdir ( outPath )

# . Get all files.
pdbFiles = glob.glob ( os.path.join ( dataPath, "*.pdb" ) )
pdbFiles.sort ( )

# . Read all PDB files.
failures  = defaultdict ( int )
successes = 0
tests     = 0
for pdbFile in pdbFiles:

    ( head, tail ) = os.path.split ( pdbFile )
    tag = tail[0:-4]

    logFile.Paragraph ( "Processing " + pdbFile + ":" )
    system = ImportSystem ( pdbFile, useComponentLibrary = True )
    BuildHydrogenCoordinates3FromConnectivity ( system )

    # . Setup.
    if system.coordinates3.numberUndefined > 0: raise ValueError ( "There are unbuilt hydrogen coordinates." )
    system.DefineMMModel ( MMModelCHARMM.WithParameterSet ( "c36a2" ) )
    system.DefineNBModel ( NBModelCutOff.WithDefaults ( ) )
    system.Summary ( )
    referenceEnergy = system.Energy ( doGradients = True )
    system.nbModel.StatisticsSummary ( system )

    # . Pickling.
    pklFile   = os.path.join ( outPath, tag +         ".pkl"  )
    yamlFile  = os.path.join ( outPath, tag +         ".yaml" )
    yamlFileM = os.path.join ( outPath, tag + "_mapping.yaml" )
    try:
        Pickle ( pklFile, system )
        successes += 1
    except Exception as e:
        failures["Pickle Failures"] += 1
        logFile.Paragraph ( "Error occurred> " +  repr ( e ) )
    try:
        YAMLPickle ( yamlFile, system )
        successes += 1
    except Exception as e:
        failures["YAML Pickle Failures"] += 1
        logFile.Paragraph ( "Error occurred> " +  repr ( e ) )
    try:
        YAMLMappingFile_FromObject ( yamlFileM, system.label, system )
        successes += 1
    except Exception as e:
        failures["YAML Mapping Pickle Failures"] += 1
        logFile.Paragraph ( "Error occurred> " +  repr ( e ) )

    # . Unpickling.
    try:
        pklSystem = Unpickle ( pklFile )
        pklSystem.label += " (Pickled)"
        pklSystem.Summary ( )
        energy = pklSystem.Energy  ( doGradients = True )
        pklSystem.nbModel.StatisticsSummary ( pklSystem )
        if math.fabs ( referenceEnergy - energy ) > _Tolerance: raise ValueError ( "Energy mismatch." )
        successes += 1
    except Exception as e:
        failures["Unpickle Failures"] += 1
        logFile.Paragraph ( "Error occurred> " +  repr ( e ) )
    try:
        yamlSystem = YAMLUnpickle ( yamlFile )
        yamlSystem.label += " (YAML Pickled)"
        yamlSystem.Summary ( )
        energy = yamlSystem.Energy  ( doGradients = True )
        yamlSystem.nbModel.StatisticsSummary ( yamlSystem )
        if math.fabs ( referenceEnergy - energy ) > _Tolerance: raise ValueError ( "Energy mismatch." )
        successes += 1
    except Exception as e:
        failures["YAML Unpickle Failures"] += 1
        logFile.Paragraph ( "Error occurred> " +  repr ( e ) )
    try:
        yamlSystemM = YAMLMappingFile_ToObject ( yamlFileM, System )
        yamlSystemM.label += " (YAML Mapping Pickled)"
        yamlSystemM.Summary ( )
        energy = yamlSystemM.Energy  ( doGradients = True )
        yamlSystemM.nbModel.StatisticsSummary ( yamlSystemM )
        if math.fabs ( referenceEnergy - energy ) > _Tolerance: raise ValueError ( "Energy mismatch." )
        successes += 1
    except Exception as e:
        failures["YAML Mapping Unpickle Failures"] += 1
        logFile.Paragraph ( "ParagraphError occurred> " +  repr ( e ) )

# . Summary of results.
f = sum ( failures.values ( ) )
n = 6 * len ( pdbFiles )
if f == 0:
    logFile.Paragraph ( "All {:d} pickling and unpickling tests were successful.".format ( n ) )
else:
    items = [ ( "Total Tests"     , "{:d}".format ( n         ) ) ,
              ( "Successes"       , "{:d}".format ( successes ) ) ]
    for key in sorted ( failures.keys ( ) ):
        f = failures[key]
        if f > 0: items.append ( ( key, "{:d}".format ( f ) ) )
    logFile.SummaryOfItems ( items, order = False, title = "Pickling and Unpickling Tests" )

# . Footer.
logFile.Footer ( )
if f != 0: TestScriptExit_Fail ( )
