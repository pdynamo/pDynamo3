"""Tests for real N-D arrays."""

import math, os, os.path, random

from pCore import logFile             , \
                  Selection           , \
                  TestScriptExit_Fail

# . Header.
logFile.Header ( )

# . Options.
_MaximumPopulationSize = 1000
_MaximumSelectionSize  = 500
_NumberOfSelections    = 50

# . Create the selection/set pairs.
population = list ( range ( 1, _MaximumPopulationSize+1 ) )
selections = []
for i in range ( _NumberOfSelections ):
    k = random.randint ( 1, _MaximumSelectionSize )
    s = set ( random.sample ( population, k ) )
    selections.append ( ( Selection.FromIterable ( s ), s ) )

# . Run the tests.
nComplement          = 0
nDifference          = 0
nIntersection        = 0
nSymmetricDifference = 0
nUnion               = 0
# . First loop over selections.
for i in range ( _NumberOfSelections ):
    ( iSelection, iSet ) = selections[i]
    # . Complement.
    test1 = ~ iSelection
    test2 = set ( range ( max ( iSelection ) +1 ) ) - iSet
    if len ( set ( test1 ) - test2 ) != 0: nComplement += 1
    # . Second loop over selections.
    for j in range ( i+1 ):
        ( jSelection, jSet ) = selections[j]
        # . Difference.
        test1 = iSelection - jSelection
        test2 = iSet       - jSet
        if len ( set ( test1 ) - test2 ) != 0: nDifference += 1
        # . Intersection.
        test1 = iSelection & jSelection
        test2 = iSet       & jSet
        if len ( set ( test1 ) - test2 ) != 0: nIntersection += 1
        # . Symmetric difference.
        test1 = iSelection ^ jSelection
        test2 = iSet       ^ jSet
        if len ( set ( test1 ) - test2 ) != 0: nSymmetricDifference += 1
        # . Union.
        test1 = iSelection | jSelection
        test2 = iSet       | jSet
        if len ( set ( test1 ) - test2 ) != 0: nUnion += 1

# . Finish up.
numberTests = _NumberOfSelections + 2 * ( _NumberOfSelections * ( _NumberOfSelections + 1 ) )
isOK        = ( ( nComplement + nDifference + nIntersection + nSymmetricDifference + nUnion ) <= 0 )
if isOK:
    logFile.Paragraph ( "All {:d} tests passed successfully.".format ( numberTests ) )
else:
    items = [ ( "Number Of Tests"                , "{:d}".format ( numberTests          ) ) ,
              ( "Complement Failures"            , "{:d}".format ( nComplement          ) ) ,
              ( "Difference Failures"            , "{:d}".format ( nDifference          ) ) ,
              ( "Intersection Failures"          , "{:d}".format ( nIntersection        ) ) ,
              ( "Symmetric Difference Failures"  , "{:d}".format ( nSymmetricDifference ) ) ,
              ( "Union Failures"                 , "{:d}".format ( nUnion               ) ) ]
    logFile.SummaryOfItems ( items, order = False, title = "Selection Test Summary" )

# . Footer.
logFile.Footer ( )
if not isOK: TestScriptExit_Fail ( )
