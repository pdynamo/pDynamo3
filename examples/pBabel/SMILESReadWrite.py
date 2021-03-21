"""Test to read and write SMILES."""

# . The SMILES module needs work.
# . No attempt is made to make sure everything is OK.

import glob, os, os.path

from Definitions import dataPath
from pBabel      import SMILESReaderError       , \
                        SMILESReader            , \
                        SMILESWriter
from pCore       import Align                   , \
                        logFile                 , \
                        LogFileActive           , \
                        TestScriptExit_Fail     , \
                        YAMLUnpickle            , \
                        YAMLPickleFileExtension

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Initialization.
isOK                   = True
errorsIn               = 0
errorsOut              = 0
smilesErrorsIn         = 0
successes              = 0
unknownSmilesErrorsIn  = 0

# . Define the SMILES to process and those which are known to give (and should give) problems.
smilesPath  = os.path.join ( dataPath, "smiles" )
badCases    = YAMLUnpickle ( os.path.join ( smilesPath, "badCases{:s}".format ( YAMLPickleFileExtension ) ) )
paths       = glob.glob ( os.path.join ( smilesPath, "*{:s}".format ( YAMLPickleFileExtension ) ) )
smiles      = {}
for path in paths:
    data = YAMLUnpickle ( path )
    smiles.update ( data )

# . Loop over the examples.
results = []
logFile.Separator ( )
for key in sorted ( smiles.keys ( ) ):
    oldSmiles = smiles[key]
    isOK      = False
    logFile.Paragraph ( "*** Processing {:s} {:s} ***".format ( key, oldSmiles ) )
    try:
        system = SMILESReader.StringToSystem ( oldSmiles )
        isOK   = True
    except SMILESReaderError as e:
        if key in badCases:
            smilesErrorsIn += 1
            tag = "Expected"
        else:
            unknownSmilesErrorsIn += 1
            tag = "Unexpected"
        logFile.Paragraph ( "{:s} SMILES input error> ".format ( tag ) + e.args[0] )
    except Exception as f:
        errorsIn += 1
        logFile.Paragraph ( "Unexpected input error occurred> " + f.args[0] )
        logFile.Paragraph ( "                                 " + oldSmiles )

    if not isOK: continue

    system.label = key
    system.Summary ( )
    formula = system.atoms.FormulaString ( )

    try:
        newSmiles = SMILESWriter.StringFromSystem ( system )#, removeExplicitHydrogens = False, useKekule = True )
        results.append ( ( key, oldSmiles, newSmiles, formula ) )
        successes += 1
    except Exception as g:
        errorsOut += 1
        logFile.Paragraph ( "Unexpected output error occurred> " + g.args[0] )
        logFile.Paragraph ( "                                  " + oldSmiles )
    logFile.Separator ( )

# . Write out the old and new SMILES.
if LogFileActive ( logFile ):
    kLength = 20
    sLength = 20
    for ( key, oldSmiles, newSmiles, _ ) in results:
        kLength = max ( kLength, len ( key ) )
        sLength = max ( sLength, len ( oldSmiles ), len ( newSmiles ) )

    # . Define a table for the results.
    table = logFile.GetTable ( columns = [ kLength + 2, 8, sLength + 2 ] )
    table.Start ( )
    table.Title ( "SMILES Comparison" )
    table.Heading ( "Name"   , columnSpan = 2 )
    table.Heading ( "SMILES" )
    for ( key, oldSmiles, newSmiles, formula ) in results:
        table.Entry ( key      , align = Align.Left )
        if oldSmiles == newSmiles:
            table.Entry ( "In/Out"  , align = Align.Left )
        else:
            table.Entry ( "In"      , align = Align.Left )
            table.Entry ( oldSmiles , align = Align.Left )
            table.Entry ( "" )
            table.Entry ( "Out"     , align = Align.Left )
        table.Entry ( newSmiles , align = Align.Left )
        table.Entry ( "" )
        table.Entry ( "Formula" , align = Align.Left )
        table.Entry ( formula   , align = Align.Left )
    table.Stop ( )

    # . Summary.
    items = [ ( "Total Tests"                    , "{:d}".format ( len ( smiles )        ) ) ,
              ( "Successes"                      , "{:d}".format ( successes             ) ) ,
              ( "Expected SMILES Input Errors"   , "{:d}".format ( smilesErrorsIn        ) ) ,
              ( "Unexpected SMILES Input Errors" , "{:d}".format ( unknownSmilesErrorsIn ) ) ,
              ( "Other Input Errors"             , "{:d}".format ( errorsIn              ) ) ,
              ( "Output Errors"                  , "{:d}".format ( errorsOut             ) ) ]
    logFile.SummaryOfItems ( items, order = False, title = "Summary of SMILES Tests" )

# . Footer.
logFile.Footer ( )
if ( not isOK ) or ( unknownSmilesErrorsIn > 0 ) or ( ( errorsIn + errorsOut ) != 0 ): TestScriptExit_Fail ( )
