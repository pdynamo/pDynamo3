"""Test for calculating MM, QC and QC/MM water dimer binding energies."""

import glob, math, os, os.path

from Definitions       import dataPath
from pBabel            import ImportSystem
from pCore             import logFile             , \
                              LogFileActive       , \
                              Selection           , \
                              TestScriptExit_Fail
from pMolecule.MMModel import MMModelOPLS
from pMolecule.NBModel import NBModelFull
from pMolecule.QCModel import QCModelMNDO
from pSimulation       import PruneByAtom

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Model names.
_Models = ( "MM", "QC" )

# . Bounds on binding energies (kJ mol^-1).
_LowerBound = -35.0
_UpperBound =   5.0

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def MonomerEnergies ( dimer, selection, nbModel, qcModel, log = logFile ):
    """Determine the monomer energies."""
    if LogFileActive ( log ): log.Heading ( "Monomer Calculation", includeBlankLine = True )
    e = {}
    # . Get the monomer.
    monomer       = PruneByAtom ( dimer, selection )
    monomer.label = "Water Monomer"
    monomer.Summary ( )
    # . MM model.
    monomer.DefineNBModel ( nbModel )
    e["MM"] = monomer.Energy ( )
    # . QC model.
    monomer.DefineQCModel ( qcModel )
    e["QC"] = monomer.Energy ( )
    return e

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Output setup.
dataPath = os.path.join ( dataPath, "mol" )

# . Define the energy models.
mmModel = MMModelOPLS.WithParameterSet ( "bookSmallExamples" )
nbModel = NBModelFull.WithDefaults ( )
qcModel = QCModelMNDO.WithDefaults ( )

# . Define the dimer with an MM model.
dimer = ImportSystem ( os.path.join ( dataPath, "waterDimer_cs.mol" ) )
dimer.DefineMMModel ( mmModel )
dimer.Summary ( )

# . Define the monomer selections.
selection1 = Selection.FromIterable ( range ( 0, 3 ) )
selection2 = Selection.FromIterable ( range ( 3, 6 ) )

# . Get the monomer energies.
e1 = MonomerEnergies ( dimer, selection1, nbModel, qcModel )
e2 = MonomerEnergies ( dimer, selection2, nbModel, qcModel )

# . Get the binding energies.
e12 = {}
for model1 in _Models:
    for model2 in _Models:
        key = model1 + " " + model2
        logFile.Heading ( model1 + "/" + model2 + " Dimer Calculation", includeBlankLine = True )
        # . Define the energy model.
        if   key == "QC QC": dimer.DefineQCModel ( qcModel )
        elif key == "QC MM": dimer.DefineQCModel ( qcModel, qcSelection = selection1 )
        elif key == "MM QC": dimer.DefineQCModel ( qcModel, qcSelection = selection2 )
        else:                dimer.DefineQCModel ( None )
        if "MM" in key: dimer.DefineNBModel ( nbModel )
        dimer.Summary ( )
        # . Store the results.
        e12[key] = dimer.Energy ( ) - e1[model1] - e2[model2]

# . Output the results.
keys = list ( e12.keys ( ) )
keys.sort ( )
table = logFile.GetTable ( columns = [ 20, 20, 20 ] )
table.Start  ( )
table.Title  ( "Water Dimer Binding Energies" )
table.Heading ( "Monomer 1" )
table.Heading ( "Monomer 2" )
table.Heading ( "Binding Energy" )
for key in keys:
    ( model1, model2 ) = key.split ( )
    table.Entry ( model1 )
    table.Entry ( model2 )
    table.Entry ( "{:.1f}".format ( e12[key] ) )
table.Stop ( )

# . Footer.
isOK = True
for e in e12.values ( ):
    if ( e < _LowerBound ) or ( e > _UpperBound ):
        isOK = False
        break
logFile.Footer ( )
if not isOK: TestScriptExit_Fail ( )
