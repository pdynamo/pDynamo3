"""Test for reading Gromacs files."""

import os

from Definitions       import dataPath                        , \
                              outPath                         , \
                              referenceDataPath
from pBabel            import ExportSystem                    , \
                              GromacsDefinitionsFileReader    , \
                              GromacsParameterFileReader      , \
                              ImportCoordinates3              , \
                              ImportSystem
from pCore             import logFile                         , \
                              LogFileActive                   , \
                              TestDataSet                     , \
                              TestScriptExit_Fail             , \
                              YAMLMappingFile_ToObject
from pMolecule         import SystemGeometryObjectiveFunction
from pMolecule.NBModel import NBModelCutOff                   , \
                              NBModelFull  
from pScientific       import Units
from pSimulation       import PrintComponentData

# . Energy values should be the same as those obtained with Gromacs 4.5 to within 1.0e-3 kJ/mole for covalent terms and within 0.1 kJ/mole for LJ/elect. terms
# . They differ from Charmm and Amber program values because of parameter manipulation and (back)conversions.

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Molecule names.
_SystemLabels = [ "ava", "1atp_peptide" ] 

# . Force field names.
_ForceFields  = [ "CHARMM", "AMBER" ]

# . Options.
_MaximumAtoms = 100

# . NB Model options (unnecessary).
#ABFS_options = { "innerCutOff" :  8.0 ,
#                 "outerCutOff" : 12.0 ,
#                 "listCutOff"  : 14.0 }

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Initialization.
observed         = {}
gromacsReference = True

# . We have two references: One for energies obtained with Gromacs; another for pDynamo energies
if gromacsReference: referenceDataPath = os.path.join ( referenceDataPath, "GromacsTopCrdRead_gromacsValues.yaml" )
else:                referenceDataPath = os.path.join ( referenceDataPath, "GromacsTopCrdRead_pDynamoValues.yaml" )

# . Output setup.
dataPath        = os.path.join ( dataPath, "gromacs" )
outPath         = os.path.join ( outPath , "gromacs" )
if not os.path.exists ( outPath ): os.mkdir ( outPath )

# . Generate systems.
for label in _SystemLabels:

    # . Force field.
    for ff in _ForceFields:

        # . Header.
        logFile.Heading ( "System \"{:s}\" with force field \"{:s}\"".format ( label, ff.lower ( ) ), includeBlankLine = True, pageWidth = 120 )

        # . Get the parameters.
        fileName            = os.path.join                                ( dataPath, label + "_" + ff )
        parameters          = GromacsParameterFileReader.PathToParameters ( fileName + ".top" )
        system              = GromacsDefinitionsFileReader.PathToSystem   ( fileName + ".top", parameters = parameters )
        system.coordinates3 = ImportCoordinates3                          ( fileName + ".gro" )
        system.label        = label
        if system.symmetryParameters is not None: system.DefineNBModel ( NBModelCutOff.WithDefaults ( ) )
        else                                    : system.DefineNBModel ( NBModelFull.WithDefaults   ( ) )
        system.Summary ( )
        energy = system.Energy ( )
        logFile.Paragraph ( "Energy (kcal/mole) = {:.4f}".format ( energy / Units.Energy_Kilocalories_Per_Mole_To_Kilojoules_Per_Mole ) )

        # . Get the dictionary of energies.
        localObserved = dict ( system.scratch.energyTerms )

        # . U-B term in Gromacs includes harmonic angle contribution, but not in pDynamo. Their sum should be equal, though.
        if gromacsReference and ( ff == "CHARMM" ):
            localObserved["Harmonic Angle + U-B"] = localObserved["Harmonic Angle"] + localObserved["Urey-Bradley"]
            del localObserved["Harmonic Angle"]
            del localObserved["Urey-Bradley"]

        # . Test gradients.
        if len ( system.atoms ) <= _MaximumAtoms:
            of = SystemGeometryObjectiveFunction.FromSystem ( system )
            localObserved["Gradient Error"] = of.TestGradients ( )

        # . Write PDB file 
        ExportSystem ( os.path.join ( outPath, label + ".pdb" ), system )
        PrintComponentData ( system.sequence, doFrequencies = False )

        # . Accumulate current data
        dataLabel = label + "-" + ff
        observed[dataLabel] = localObserved

# . Verify the observed data against the reference data.
referenceData = YAMLMappingFile_ToObject ( referenceDataPath, TestDataSet )
results       = referenceData.VerifyAgainst ( observed )
isOK          = results.WasSuccessful ( )
results.Summary ( fullSummary = True )

# . Footer.
logFile.Footer ( )
if not isOK: TestScriptExit_Fail ( )
