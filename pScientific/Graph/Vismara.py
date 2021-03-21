"""Calculate the relevant cycles of an undirected graph."""

from  pCore                 import Clone
from .BiconnectedComponents import BiconnectedComponents
from .Dijkstra              import DijkstraSingleSource
from .Edge                  import Edge
from .GaussElimination      import GaussElimination
from .GraphStatus           import GraphError

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def VismaraRelevantCycles ( graph, biconnectedComponents = None ):
    """Calculate the relevant cycles of an undirected graph."""

    # . Initialization.
    cycles = []

    # . Get the biconnected components if necessary.
    if biconnectedComponents is None:
        biconnectedComponents = BiconnectedComponents ( graph )

    # . Determine the cycles for each component.
    for component in biconnectedComponents:

        # . Get a graph corresponding to the component.
        subgraph = graph.MakeSubgraph ( component )

        # . Get prototypes.
        prototypes = _CalculatePrototypes ( subgraph )

        # . Calculate relevant cycles.
        relevants = _CalculateRelevantCycles ( prototypes )

        # . Add into cycles.
        cycles.extend ( relevants )

    # . Finish up.
    return cycles

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def _AllDirectedPaths ( graph, start, stop, currentPath = None ):
    """Find all directed paths between nodes start and stop."""
    # . Initialization.
    if currentPath is None: myCurrentPath = []
    else:                   myCurrentPath = currentPath
    myCurrentPath.append ( start )
    results = []
    # . All paths found.
    if start is stop:
        results.append ( myCurrentPath )
    # . Check for all nodes that have directed edges from start.
    else:
        for node in graph.nodes:
            if _NodeIsReachableFrom ( graph, node, start ):
                extras = _AllDirectedPaths ( graph, node, stop, currentPath = list ( myCurrentPath ) )
                results.extend ( extras )            
    # . Finish up.
    return results

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def _CalculatePrototypes ( graph ):
    """Calculate prototypes."""

    # . Initialization.
    prototypes = []

    # . Loop over all nodes.
    for ( rIndex, r ) in enumerate ( graph.nodes ):

        # . Skip the zero graph.
        if rIndex == 0: continue

        # . Make a subgraph containing all nodes of lesser order.
        reduced = graph.MakeSubgraph ( graph.nodes[0:rIndex+1] )

        # . Shortest paths from r to all other nodes in the reduced graph.
        ( distances, paths ) = DijkstraSingleSource ( reduced, r )

        # . Get Vr.
        Vr = []
        for node in reduced.nodes[0:-1]:
            if len ( paths.get ( node, [] ) ) > 0: Vr.append ( node )
        Vr.append ( r )

        if len ( Vr ) <= 1: continue

        # . Make a subgraph with reachable nodes.
        digraph = graph.MakeSubgraph ( Vr, induced = False )

        # . Remove r from Vr.
        Vr.pop ( -1 )

        # . Loop over all nodes from Vr.
        for y in Vr:

            # . Initialization.
            s = []
            yDistance = distances[y]
            yPath     = paths[y]

            # . For all neighbors to y that occur in Vr.
            for z in reduced.adjacentNodes.get ( y, [] ):
                if z is r:
                    digraph.AddEdge ( Edge.WithNodes ( y, r ) )
                elif z in Vr:

                    # . Initialization.
                    zDistance = distances[z]
                    zPath     = paths[z]

                    # . z belongs to a shortest path from r to y.
                    if ( zDistance + 1 ) == yDistance:
                        s.append ( z )
                        try: digraph.AddEdge ( Edge.WithNodes ( y, z ) )
                        except: pass

                    # . Odd cycle - P(r,z) + (z,y) + P(y,r).
                    # . id should be sufficient here for ordering as long as nodes constant.
                    elif ( id ( z ) < id ( y ) ) and ( zDistance != yDistance + 1 ) and _PathIntersectsAt ( yPath, zPath, r ):
                        newPath = zPath + yPath[::-1][:-1]
                        prototypes.append ( ( len ( newPath ), newPath, ( r, y, z, len ( yPath ), len ( zPath ) ), digraph, _MakeBooleanEdgeVector ( graph, newPath ) ) )

            # . Loop over nodes in s.
            n = len ( s )
            for ( pIndex, p ) in enumerate ( s ):
                pPath = paths[p]
                for q in s[pIndex+1:]:
                    qPath = paths[q]
                    # . Even cycle - P(r,p) + (p,y) + (y,q) + P(q,r).
                    if _PathIntersectsAt ( pPath, qPath, r ):
                        newPath = pPath + [ y ] + qPath[::-1][:-1]
                        prototypes.append ( ( len ( newPath ), newPath, ( r, p, q, y, len ( pPath ), len ( qPath ) ), digraph, _MakeBooleanEdgeVector ( graph, newPath ) ) )

    # . Sort the prototypes - in order of length.
    prototypes.sort ( )

    # . Process the prototypes.
    # . Get the size of the cycle basis space.
    cyclomatic = len ( graph.edges ) - len ( graph.nodes ) + 1

    # . Initialization.
    equalCycles   = []
    lesserCycles  = []
    length        =  0
    newPrototypes = []

    # . Loop over prototypes.
    for prototype in prototypes:

        # . Get the cycle.
        cycle           = prototype[1]
        cycleEdgeVector = prototype[-1]
        newLength       = len ( cycle )

        # . Check for a new length.
        if newLength != length:

            # . Still missing cycles in the MCB.
            if cyclomatic > 0:
                lesserCycles.extend ( equalCycles )
                length      = newLength
                equalCycles = []
            # . All cycles for MCB found so finish.
            else:
                break

        # . Check for linear dependence wrt shorter cycles.
        testCycles = Clone ( lesserCycles )
        testCycles.append ( Clone ( cycleEdgeVector ) )
        rank  = len ( testCycles )
        rankg = GaussElimination ( testCycles, 0, 0 )

        # . Cycle is independent.
        if ( rank == 1 ) or ( rank == rankg ):

            # . Save the prototype.
            newPrototypes.append ( prototype )

            # . Check for independence versus cycles of equal length.
            testCycles = Clone ( equalCycles )
            testCycles.append ( Clone ( cycleEdgeVector ) )
            rank  = len ( testCycles )
            rankg = GaussElimination ( testCycles, 0, 0 )

            # . Cycle independent.
            if ( rank == 1 ) or ( rank == rankg ):
                equalCycles.append ( cycleEdgeVector )
                cyclomatic -= 1

#    print "\nNumber of pruned prototypes = ", len ( prototypes ) - len ( newPrototypes )

    # . Sort the prototypes - in order of length.
    newPrototypes.sort ( )
    return newPrototypes

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def _CalculateRelevantCycles ( prototypes ):
    """Calculate the relevant cycles given the prototypes."""

    # . Initialization.
    relevants = []

    # . Loop over prototypes.
    for prototype in prototypes:

        # . Get the cycle.
        cycle  = prototype[1]
        length = len ( cycle )

        # . There are no other cycles of length 3 and 4 in the prototype family as all nodes are adjacent.
        if length < 5:
            relevants.append ( cycle )

        # . Calculate relevant cycles from the prototype.
        else:

            # . Get node indices.
            nodes   = prototype[2]
            r       = nodes[0]
            p       = nodes[1]
            q       = nodes[2]
            digraph = prototype[3]
            isEven  = ( len ( nodes ) > 5 )
            pLength = nodes[-2]
            qLength = nodes[-1]

            # . All paths P(p,r) and P(q,r).
            pPaths = _AllDirectedPaths ( digraph, p, r )
            qPaths = _AllDirectedPaths ( digraph, q, r )
            for ( pathLength, paths ) in ( ( pLength, pPaths ), ( qLength, qPaths ) ):
                for path in paths:
                    xnodes = set ( path )
                    if len ( xnodes ) != len ( path ): raise GraphError ( "\nInvalid directed path: {:d} {:d} {:d}.".format ( len ( xnodes ), len ( path ), pathLength ) )

            # . Even cycles.
            if isEven:
                x = nodes[3]
                for pPath in pPaths:
                    for qPath in qPaths:
                        relevants.append ( pPath[::-1] + [ x ] + qPath[:-1] )
            # . Odd cycles.
            else:
                for pPath in pPaths:
                    for qPath in qPaths:
                        relevants.append ( pPath[::-1] + qPath[:-1] )

#    print "Number relevants - number prototypes = ", len ( relevants ) - len ( prototypes )

    #. Finish up.
    relevants.sort ( key = len )
    return relevants

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def _MakeBooleanEdgeVector ( graph, path ):
    """Return a list of booleans in graph-edge order indicating which edges occur in the path."""
    vector = [ False for i in range ( len ( graph.edges ) ) ]
    for ( i, tail ) in enumerate ( path ):
        if i == len ( path ) - 1: head = path[0]
        else:                     head = path[i+1]
        found = False
        for ( e, edge ) in enumerate ( graph.edges ):
            if ( ( edge.node1 is head ) and ( edge.node2 is tail ) ) or ( ( edge.node2 is head ) and ( edge.node1 is tail ) ):
                found     = True
                vector[e] = True
                break
        if not found: raise GraphError ( "No edges connect two nodes." )
    return vector

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def _NodeIsReachableFrom ( graph, targetNode, sourceNode ):
    """Is the target node reachable from source node?"""
    isReachable = False
    for edge in graph.edges:
        if ( sourceNode is edge.node1 ) and ( targetNode is edge.node2 ):
            isReachable = True
            break
    return isReachable

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def _PathIntersectsAt ( aPath, bPath, node ):
    """Check to see if two paths intersect at node."""
    common = set ( aPath ).intersection ( set ( bPath ) )
    return ( ( len ( common ) == 1 ) and ( common.pop ( ) is node ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
