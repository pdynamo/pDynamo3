"""Test SQL atom selection."""

import glob, os

from Definitions       import dataPath
from pBabel            import ImportSystem
from pCore             import logFile             , \
                              TestScriptExit_Fail
from pMolecule         import AtomSelection       , \
                              AtomSelectionError  , \
                              SQLAtomSelector
from pMolecule.MMModel import MMModelOPLS
from pMolecule.NBModel import NBModelCutOff
from pSimulation       import BuildHydrogenCoordinates3FromConnectivity

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Paths.
dataPath = os.path.join ( dataPath, "pdb" )

# . Get all files.
pdbFiles = glob.glob ( os.path.join ( dataPath, "*.pdb" ) )
pdbFiles.sort ( )

# . Read all PDB files.
numberFailed = 0
for pdbFile in pdbFiles:

    logFile.Paragraph ( "Processing " + pdbFile + ":" )
    system = ImportSystem ( pdbFile, useComponentLibrary = True )
    BuildHydrogenCoordinates3FromConnectivity ( system )
    try:

        # . Setup.
        if system.coordinates3.numberUndefined > 0: raise
        system.DefineMMModel ( MMModelOPLS.WithParameterSet ( "protein" ) )
        system.DefineNBModel ( NBModelCutOff.WithDefaults ( ) )
        system.Summary ( )
        system.Energy  ( doGradients = True )

        # . Selection.
        selector              = SQLAtomSelector.WithSystem ( system )

        # . Standard selections.
        aromatics             = selector.aromatics
        backbone              = selector.backbone
        boundaryAtoms         = selector.boundaryAtoms
        counterions           = selector.counterions
        heavyAtoms            = selector.heavyAtoms
        hydrogens             = selector.hydrogens
        mmAtoms               = selector.mmAtoms
        polymerAtoms1         = selector.linearPolymerAtoms
        protein               = selector.protein
        qcAtoms               = selector.qcAtoms
        ringAtoms             = selector.ringAtoms
        water                 = selector.water

        # . Where selections.
        nearOrigin1           = selector.Where ( "X*X + Y*Y + Z*Z < 25.0" )
        positive              = selector.Where ( "Charge > 0.0" )
        threonines1           = selector.Where ( "Path LIKE '%:THR.%:%'" )
        threonines2           = selector.Where ( "ResNam='THR'" )

        # . Atom selection methods.
        nearOrigin2           = AtomSelection.Within           ( system, nearOrigin1 , 5.0 )
        nearOrigin2           = AtomSelection.ByComponent      ( system, nearOrigin2 )
        polymerAtoms2         = AtomSelection.ByLinearPolymer  ( system, threonines1 )
        neighbors             = AtomSelection.ByBondedNeighbor ( system, threonines1 , iterations = 3 )

        # . Selection operators.
        # . Remember no upper bound here so relies on maximum index in nearOrigin2.
        complementNearOrigin2 = ~ nearOrigin2
        null                  = ( nearOrigin2 & complementNearOrigin2 )
        total1                = ( nearOrigin2 | complementNearOrigin2 )
        total2                = ( nearOrigin2 ^ complementNearOrigin2 )

        # . Basic checks.
        n    = len ( system.atoms )
        isOK = ( len ( boundaryAtoms ) == 0 ) and \
               ( len ( qcAtoms       ) == 0 ) and \
               ( len ( mmAtoms       ) == n ) and \
               ( len ( heavyAtoms    )  + len ( hydrogens     ) == n ) and \
               ( len ( polymerAtoms1 ) == len ( polymerAtoms2 )      ) and \
               ( len ( counterions   )  + len ( protein       ) + len ( water ) == n ) and \
               ( len ( threonines1   ) == len ( threonines2   )      ) and \
               ( len ( null          ) == 0 ) and \
               ( len ( total1        ) == len ( total2        )      ) and \
               ( len ( total1        ) == max ( total1        ) +  1 )
        if not isOK: raise AtomSelectionError ( "Atom selection error." )

    except Exception as e:
        numberFailed += 1
        logFile.Paragraph ( "Error occurred> " +  e.args[0] )

# . Footer.
logFile.SummaryOfItems ( [ ( "Successes", "{:d}".format ( len ( pdbFiles ) - numberFailed ) )   ,
                           ( "Failures" , "{:d}".format (                    numberFailed ) ) ] ,
                           order = False, title = "SQL Atom Selection Tests" )
logFile.Footer ( )
if numberFailed != 0: TestScriptExit_Fail ( )
