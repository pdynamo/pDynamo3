"""Testing for relevant cycles."""

import glob, itertools, math, os, os.path

from collections       import defaultdict
from Definitions       import dataPath
from pCore             import Align                    , \
                              logFile                  , \
                              TestScript_InputDataPath , \
                              TestScriptExit_Fail      , \
                              YAMLUnpickle
from pScientific.Graph import BiconnectedComponents    , \
                              ConnectedComponents      , \
                              Edge                     , \
                              Graph                    , \
                              HanserAllCycles          , \
                              Node                     , \
                              PatonMinimalCycleBasis   , \
                              VismaraRelevantCycles

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Options.
_MaximumCyclomaticNumber =    20   # . For all cycles perception.
_MaximumDegree           = 10000
_TestFigueras            =  True
_TestHanser              =  True
_TestPaton               =  True
_TestRelevant            =  True

if _TestFigueras:
    from pScientific.Graph import FiguerasRings

# . Utility method.
def _CycleFrequencyString ( cycles ):
    lengths = defaultdict ( int )
    for cycle in cycles: lengths[len(cycle)] += 1
    items = [ "{:d} x {:d}".format ( lengths[key], key ) for key in sorted ( lengths.keys ( ) ) ]
    return "{:s}".format ( " ; ".join ( items ) )

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Get data.
data             = YAMLUnpickle ( os.path.join ( dataPath, "BiconnectedGraphExamples.yaml" ) )
connectivities   = data["Connectivities" ]
hanserCycles     = data["Hanser Cycles"  ]
relevantCycles   = data["Relevant Cycles"]

# . Initialization.
hanserFailures = []
successes      = 0
totalTests     = 0

# . Loop over examples.
for key in sorted ( connectivities.keys ( ) ):
    #if key != "Figueras 1": continue
    connectivity = connectivities[key]

    # . Create graph.
    graph = Graph ( )

    # . Nodes.
    nodes = {}
    for i in range ( len ( connectivity ) ):
        node = Node ( )
        nodes[i] = node
        graph.AddNode ( node )

    # . Edges.
    for ( i, connections ) in enumerate ( connectivity ):
        for j in connections:
            if j > i: graph.AddEdge ( Edge.WithNodes ( nodes[i], nodes[j] ) )

    # . Various graph statistics.
    connectedComponents   =   ConnectedComponents ( graph )
    biconnectedComponents = BiconnectedComponents ( graph )

    # . Output.
    CN      = len ( graph.edges ) - len ( graph.nodes ) + len ( connectedComponents )
    details = []
    items = [ ( "Number of nodes"             , "{:d}".format ( len ( graph.nodes           ) ) ) ,
              ( "Number of edges"             , "{:d}".format ( len ( graph.edges           ) ) ) ,
              ( "Connected components"        , "{:d}".format ( len ( connectedComponents   ) ) ) ,
              ( "Biconnected components"      , "{:d}".format ( len ( biconnectedComponents ) ) ) ,
              ( "Cyclomatic number"           , "{:d}".format ( CN        ) ) ,
              ( "Upper limit to cycle number" , "{:d}".format ( 2**CN - 1 ) ) ]

    # . Figueras SSSR - for old time's sake.
    if _TestFigueras:
        table = {}
        for ( i, connections ) in enumerate ( connectivity ):
            if len ( connections ) > 0: table[i] = connections
        cycles = FiguerasRings ( table )
        if len ( cycles ) > 0:
            details.append ( ( "Figueras SSSR", _CycleFrequencyString ( cycles ) ) )
            items.append   ( ( "Figueras SSSR", "{:d}".format ( len ( cycles ) ) ) )

    # . All cycles.
    if _TestHanser and ( CN <= _MaximumCyclomaticNumber ):
        ( cycleSets, isOK ) = HanserAllCycles ( graph, biconnectedComponents = biconnectedComponents, maximumDegree = _MaximumDegree )
        if len ( cycleSets ) > 0:
            cycles = list ( itertools.chain.from_iterable ( cycleSets ) )
            string = _CycleFrequencyString ( cycles )
            if isOK:
                details.append ( ( "Hanser all cycles", string ) )
                items.append   ( ( "Hanser all cycles", "{:d}".format ( len ( cycles ) ) ) )
            else:
                details.append ( ( "Hanser all cycles (partial)", string ) )
                items.append   ( ( "Hanser all cycles (partial)", "{:d}".format ( len ( cycles ) ) ) )
            reference = hanserCycles.get ( key, None )
            if reference is not None:
                if len ( cycles ) == reference: successes += 1
                else: hanserFailures.append ( ( key, len ( cycles ), reference ) )
                totalTests += 1

    # . Paton MCB.
    if _TestPaton:
        cycleSets = PatonMinimalCycleBasis ( graph )
        if len ( cycleSets ) > 0:
            cycles = list ( itertools.chain.from_iterable ( cycleSets ) )
            details.append ( ( "Paton MCB cycles", _CycleFrequencyString ( cycles ) ) )
            items.append   ( ( "Paton MCB cycles", "{:d}".format ( len ( cycles ) ) ) )
            totalTests += 1
            if len ( cycles ) == CN: successes += 1

    # . Vismara relevant cycles.
    if _TestRelevant:
        cycles = VismaraRelevantCycles ( graph, biconnectedComponents = biconnectedComponents )
        if len ( cycles ) > 0:
            details.append ( ( "Vismara relevant cycles", _CycleFrequencyString ( cycles ) ) )
            items.append   ( ( "Vismara relevant cycles", "{:d}".format ( len ( cycles ) ) ) )
            reference = relevantCycles.get ( key, None )
            if reference is not None:
                vismaraLengths   = sorted ( [ len ( cycle ) for cycle in cycles ] )
                referenceLengths = []
                for ( length, frequency ) in reference.items ( ):
                    referenceLengths.extend ( frequency * [ length ] )
                referenceLengths.sort ( )
                if referenceLengths == vismaraLengths: successes += 1
                totalTests += 1

    # . Summary and cycles.
    logFile.SummaryOfItems ( items, order = False, title = "Graph Example - {:s}".format ( key ) )
    table = logFile.GetTable ( columns = [ max ( [ len ( k ) for ( k, _ ) in details ] ) + 2, max ( [ len ( v ) for ( _, v ) in details ] ) ] )
    table.Start  ( )
    table.Title  ( "Cycle Decompositions" )
    for ( key, value ) in details:
        table.Entry ( key  , align = Align.Left )
        table.Entry ( value, align = Align.Left )
    table.Stop ( )

# . Failures.
if len ( hanserFailures ) > 0:
    n     = max ( [ len ( key ) for ( key, _, _ ) in hanserFailures ] ) + 1
    table = logFile.GetTable ( columns = [ n, 12, 12 ] )
    table.Start  ( )
    table.Title  ( "Hanser Failures" )
    table.Heading ( "Name"      )
    table.Heading ( "Cycles"    )
    table.Heading ( "Reference" )
    for ( key, cycles, reference ) in hanserFailures:
        table.Entry ( key, align = Align.Left )
        table.Entry ( "{:d}".format ( cycles    ) )
        table.Entry ( "{:d}".format ( reference ) )
    table.Stop ( )

# . Final printing.
logFile.Paragraph ( "There were {:d} successful tests out of {:d} total tests on {:d} connectivities.".format ( successes, totalTests, len ( connectivities) ) )
if successes == totalTests: logFile.Paragraph ( "Test passed." )
else:                       logFile.Paragraph ( "Test failed." )

# . Footer.
logFile.Footer ( )
if successes != totalTests: TestScriptExit_Fail ( )
