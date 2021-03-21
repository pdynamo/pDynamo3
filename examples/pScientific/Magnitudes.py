"""Magnitude tests."""

# . Notes:
#
#   before, after, bounds were originally added directly to enum using decorator.
#   This was complicated and these quantities are not so difficult to calculate
#   directly when needed.
#
#   Additionally, comparison was added using __eq__, __ne__, etc. but again this
#   proved unnecessary. Incidentally, explicit definition of __eq__ requires that
#   __hash__ also be provided explicitly if the enums are to be used as keys.
#

import math, random

from pCore       import logFile             , \
                        TestScriptExit_Fail
from pScientific import Magnitude           , \
                        Magnitude_Adjust

# . Options.
_Powers    =  25
_Tolerance = 1.0e-12 # . Maximum that can be expected for large numbers is ~1.0e-13.
_Tries     = 100

# . Header.
logFile.Header ( )

# . Magnitude names.
names = [ magnitude.name for magnitude in Magnitude.__members__.values ( ) ]

# . Summary of enum.
# . Gather data.
data  = [ ( "-", "-", "-", "{:s}".format ( repr ( 0.0 ) ) ) ]
for ( i, ( name, magnitude ) ) in enumerate ( Magnitude.__members__.items ( ) ):
    if i <= 0                :  before = "-"
    else:                       before = names[i-1]
    if i >= len ( names ) - 1: after   = "-"
    else:                      after   = names[i+1]
    bound = "{:s}".format ( repr ( math.pow ( 10.0, magnitude.power ) ) )
    data.append ( ( name, before, after, bound ) )
data.append ( ( "-", "-", "-", "{:s}".format ( repr ( math.inf ) ) ) )

# . Printing.
table = logFile.GetTable ( columns = [ 20, 20, 20, 20 ] )
table.Start   ( )
table.Title   ( "Magnitudes"  )
table.Heading ( "Short Name"  )
table.Heading ( "Predecessor" )
table.Heading ( "Successor"   )
table.Heading ( "Value"       )
for items in data:
    for item in items: table.Entry ( item )
table.Stop ( )

# . Test magnitude adjustment.
failed = 0
powers = list ( range ( -_Powers, _Powers+1 ) )
table  = logFile.GetTable ( columns = [ 22, 3, 22, 3, 20 ] )
table.Start   ( )
table.Title   ( "Magnitude Adjustment"  )
table.Heading ( "Before"  , columnSpan = 2 )
table.Heading ( "After"   , columnSpan = 2 )
table.Heading ( "% Error" )
for t in range ( _Tries ):
    p  = random.choice  ( powers )
    m0 = Magnitude[random.choice ( names )]
    v0 = random.uniform ( math.pow ( 10.0, p ), math.pow ( 10.0, p+1 ) )
    if random.choice ( [ "+", "-" ] ) == "-": v0 *= -1.0
    ( v1, m1 ) = Magnitude_Adjust ( v0, magnitude = m0, powersOfThree = True )#False )
    n0 = v0 * math.pow ( 10.0, m0.power )
    n1 = v1 * math.pow ( 10.0, m1.power )
    d  = math.fabs ( n0 - n1 ) / math.fabs ( n0 )
    if d > _Tolerance: failed += 1
    table.Entry ( "{:.10g}".format ( v0 ) ) ; table.Entry ( m0.symbol )
    table.Entry ( "{:.10g}".format ( v1 ) ) ; table.Entry ( m1.symbol )
    table.Entry ( "{:.4g}".format  ( d ) )
table.Stop ( )

# . Footer.
if failed == 0: logFile.Paragraph ( "Test Succeeded." )
else:           logFile.Paragraph ( "Test Failed."    )
logFile.Footer ( )
if failed > 0: TestScriptExit_Fail ( )
