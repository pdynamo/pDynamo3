"""Testing for point groups."""

import glob, math, os

from Definitions          import structuresPath
from pBabel               import ImportSystem
from pCore                import Align               , \
                                 logFile             , \
                                 LogFileActive       , \
                                 TestScriptExit_Fail
from pScientific.Arrays   import Array
from pScientific.Symmetry import Find3DGraphPointGroup

# . This test is quite rapid apart from the icosahedral point groups (I, Ih)
# . which take > 97% of the time.

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Options.
_doLong     = False
_doPrinting = True

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Initialization.
log = logFile
if not _doPrinting: log = None

# . Paths.
paths = glob.glob ( os.path.join ( structuresPath, "pointGroupExamples", "xyz", "*.xyz" ) )
paths.sort ( )

# . Calculation.
results   = []
successes = 0
for path in paths:

    ( head, tail ) = os.path.split ( path )
    inGroup = tail[0:-4].split ( "_" )[0]

    # . Skip icosahedral cases for short jobs.
    if ( not _doLong ) and ( ( inGroup == "I" ) or ( tail[0:-4] == "Ih_c" ) ): continue

    system = ImportSystem ( path )
    system.Summary ( log = log )

    masses       = Array.FromIterable ( [ atom.mass for atom in system.atoms ] )
    numbers      = [ atom.atomicNumber for atom in system.atoms ]
    localResults = Find3DGraphPointGroup ( numbers                                ,
                                           system.coordinates3                    ,
                                           doCharacterSymmetryOperations = False  ,
                                           log                           = log    ,
                                           weights                       = masses )
    pointGroup = localResults.get ( "Point Group", None )
    if pointGroup is None: outGroup = None
    else:                  outGroup = pointGroup.label

    if outGroup is None:
        if LogFileActive ( log ):
            log.Paragraph ( "No point group found for " + system.label + "." )
        outGroup = "-"
    else:
        if LogFileActive ( log ):
            log.Paragraph ( "Point group " + outGroup + " for " + system.label + "." )
        if inGroup == outGroup: successes += 1

    results.append ( ( inGroup, outGroup, system.label ) )

# . Printing.
if LogFileActive ( log ):
    n     = max ( [ len ( label ) for ( _, _, label ) in results ] ) + 2
    table = logFile.GetTable ( columns = [ 2, 12, 12, n ] )
    table.Start  ( )
    table.Title  ( "Summary of Results" )
    table.Heading ( "OK"       )
    table.Heading ( "Expected" )
    table.Heading ( "Obtained" )
    table.Heading ( "System"   )
    for ( inGroup, outGroup, label ) in results:
        if inGroup == outGroup: tag = ""
        else:                   tag = "*"
        table.Entry ( tag     , align = Align.Left )
        table.Entry ( inGroup , align = Align.Left )
        table.Entry ( outGroup, align = Align.Left )
        table.Entry ( label )
    table.Stop ( )
    log.Paragraph ( "Number of successes = {:d} (out of {:d}).".format ( successes, len ( results ) ) )
    if successes == len ( results ): log.Paragraph ( "Test Succeeded." )
    else:                            log.Paragraph ( "Test Failed."    )

# . Footer.
logFile.Footer ( )
if successes != len ( results ): TestScriptExit_Fail ( )
