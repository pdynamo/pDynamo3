"""Maximum matches in graphs."""

import glob, math, os, os.path, sys

from pCore             import Align                    , \
                              logFile                  , \
                              TestScript_InputDataPath , \
                              TestScriptExit_Fail      , \
                              YAMLUnpickle
from pScientific.Graph import BiconnectedComponents    , \
                              ConnectedComponents      , \
                              Edge                     , \
                              EdmondsMaximumMatching   , \
                              Graph                    , \
                              Node

# . Header.
logFile.Header ( )

# . Get data.
dataPath       = os.path.join ( os.getenv ( "PDYNAMO3_HOME" ), "examples", "pScientific", "data" )
data           = YAMLUnpickle ( os.path.join ( dataPath, "BiconnectedGraphExamples.yaml" ) )
connectivities = data["Connectivities" ]

# . Initialization.
oddGraphs      = 0
perfectMatches = 0
totalTests     = len ( connectivities )

# . Loop over examples.
for key in sorted ( connectivities.keys ( ) ):

    # . Create graph.
    connectivity = connectivities[key]
    graph        = Graph ( )
    nodes        = {}
    for i in range ( len ( connectivity ) ):
        node = Node ( )
        nodes[i] = node
        graph.AddNode ( node )
    for ( i, connections ) in enumerate ( connectivity ):
        for j in connections:
            if j > i: graph.AddEdge ( Edge.WithNodes ( nodes[i], nodes[j] ) )

    # . Various graph statistics.
    connectedComponents   =   ConnectedComponents ( graph )
    biconnectedComponents = BiconnectedComponents ( graph )

    # . Output.
    CN    = len ( graph.edges ) - len ( graph.nodes ) + len ( connectedComponents )
    items = [ ( "Number of nodes"        , "{:d}".format ( len ( graph.nodes           ) ) ) ,
              ( "Number of edges"        , "{:d}".format ( len ( graph.edges           ) ) ) ,
              ( "Connected components"   , "{:d}".format ( len ( connectedComponents   ) ) ) ,
              ( "Biconnected components" , "{:d}".format ( len ( biconnectedComponents ) ) ) ,
              ( "Cyclomatic number"      , "{:d}".format ( CN ) ) ]
    logFile.SummaryOfItems ( items, order = False, title = "Graph Example - {:s}".format ( key ) )

    # . Matching.
    result = EdmondsMaximumMatching ( graph )
    nMatch = len ( result      )
    nNodes = len ( graph.nodes )
    if nNodes % 2 != 0: oddGraphs += 1
    nodeIndex = { node : order for ( order, node ) in enumerate ( graph.nodes ) }
    pairs     = set ( )
    for ( node, match ) in result.items ( ):
        i = nodeIndex[node ]
        j = nodeIndex[match]
        pairs.add ( ( min ( i, j ), max ( i, j ) ) )
    pairs = sorted ( pairs )
    if nMatch == nNodes:
        title = "Perfect match of length {:d}".format ( nMatch )
        perfectMatches += 1
    else:
        title = "Maximum match of length {:d} / {:d}".format ( nMatch, nNodes )
    table = logFile.GetTable ( columns = 5 * [ 4, 1, 4 ] )
    table.Start  ( )
    table.Title  ( title )
    for ( i, j ) in pairs:
        table.Entry ( "{:d}".format ( i ), align = Align.Right )
        table.Entry ( "-" )
        table.Entry ( "{:d}".format ( j ), align = Align.Left  )
    table.Stop ( )

# . Final printing.
isOK  = ( ( oddGraphs + perfectMatches ) == totalTests )
items = [ ( "Total tests"     , "{:d}".format ( totalTests     ) ) ,
          ( "Odd graphs"      , "{:d}".format ( oddGraphs      ) ) ,
          ( "Perfect matches" , "{:d}".format ( perfectMatches ) ) ]
logFile.SummaryOfItems ( items, order = False, title = "Final Results" )
if isOK: logFile.Paragraph ( "Test passed." )
else:    logFile.Paragraph ( "Test failed." )

# . Footer.
logFile.Footer ( )
if not isOK: TestScriptExit_Fail ( )
