"""Test pickling and unpickling."""

import glob, math, os, os.path

from collections          import defaultdict
from Definitions          import dataPath                   , \
                                 outPath
from pBabel               import ImportSystem
from pCore                import logFile                    , \
                                 Pickle                     , \
                                 Selection                  , \
                                 TestScriptExit_Fail        , \
                                 Unpickle                   , \
                                 YAMLMappingFile_FromObject , \
                                 YAMLMappingFile_ToObject   , \
                                 YAMLPickle                 , \
                                 YAMLUnpickle
from pMolecule            import AtomSelection              , \
                                 System
from pMolecule.MMModel    import MMModelCHARMM
from pMolecule.NBModel    import NBModelCutOff
from pMolecule.QCModel    import QCModelMNDO
from pScientific.Symmetry import CrystalSystemCubic         , \
                                 PeriodicBoundaryConditions

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Scalars.
_BoxSize     = 28.0
_Destination = "picklingQCMM"
_Tolerance   = 0.1

# . Atom names.
_Tags = ( "CB",  "CG",  "CD1",  "CD2",  "CE1",  "CE2",  "CZ",  "OH",  "HB2",  "HB3",  "HD1",  "HD2",  "HE1",  "HE2",  "HH" )

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Initialization.
failures  = defaultdict ( int )
successes = 0
tests     = 0

# . Paths.
dataPath = os.path.join ( dataPath, "pdb"        )
outPath  = os.path.join ( outPath , _Destination )
if not os.path.exists ( outPath ): os.mkdir ( outPath )

# . Models.
mmModel = MMModelCHARMM.WithParameterSet ( "c36a2" )
nbModel = NBModelCutOff.WithDefaults ( )
qcModel = QCModelMNDO.WithDefaults   ( )

# . Get the file.
pdbFile = os.path.join ( dataPath, "2E4E_folded_solvated.pdb" )
( head, tail ) = os.path.split ( pdbFile )
tag = tail[0:-4]

logFile.Paragraph ( "Processing " + pdbFile + ":" )
system = ImportSystem ( pdbFile, useComponentLibrary = True )

# . Free atoms.
freeAtoms = AtomSelection.FromAtomPattern ( system, "A:*:*" )

# . QC selection.
indices = set ( )
for atomTag in _Tags:
    indices.add ( system.sequence.AtomIndex ( "A:TYR.2:" + atomTag ) )
tyrosine = Selection ( indices )

# . Setup.
system.symmetry           = PeriodicBoundaryConditions.WithCrystalSystem ( CrystalSystemCubic ( ) )
system.symmetryParameters = system.symmetry.MakeSymmetryParameters ( a = _BoxSize )
system.freeAtoms          = freeAtoms
system.DefineMMModel   ( mmModel )
system.DefineQCModel   ( qcModel, qcSelection = tyrosine )
system.DefineNBModel   ( nbModel )
system.Summary ( )
referenceEnergy = system.Energy  ( doGradients = True )
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
    logFile.Paragraph ( "Error occurred> " +  repr ( e ) )

# . Summary of results.
f = sum ( failures.values ( ) )
n = 6
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
