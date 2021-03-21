"""Solve the single source shortest path problem for a graph with non-negative edge weights."""

import heapq # . A priority queue.

from .GraphStatus import GraphError

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def DijkstraShortestPaths ( graph, source, reverse = False, target = None ):
    """Returns shortest paths in a weighted graph."""
    # . Check edge weights.
    if graph.MinimumEdgeWeight ( ) < 0: raise GraphError ( "Graph has negative edge weights." )
    # . Edges to use.
    if reverse: outEdges = graph.adjacentEdges # inEdges
    else:       outEdges = graph.adjacentEdges # outEdges
    # . Initialization.
    predecessorEdges = { source : None }
    predecessorNodes = { source : None }
    weights          = {}
    priorityQueue    = []
    seen             = { source : 0.0 }
    heapq.heappush ( priorityQueue, ( 0.0, source ) )
    # . Loop until the queue is empty.
    while len ( priorityQueue ) > 0:
        ( weight, node ) = heapq.heappop ( priorityQueue )
        if node in weights: continue
        weights[node]    = weight
        if node is target: break
        for edge in outEdges.get ( node, [] ):
            other     = edge.Opposite ( node )
            newWeight = weight + edge.weight
            if other in weights:
                if newWeight < weights[other]:
                    raise GraphError ( "Logic error: shorter path found." )
            elif ( other not in seen ) or ( newWeight < seen[other] ):
                heapq.heappush ( priorityQueue, ( newWeight, other ) )
                predecessorEdges[other] = edge
                predecessorNodes[other] = node
                seen            [other] = newWeight
    # . Build paths.
    paths = {}
    for node in weights:
        paths[node] = PathHead ( headNode = node )
    for node in weights:
        paths[node].UpdateTail ( headEdge = predecessorEdges[node], tail = paths.get ( predecessorNodes[node], None ) )
    for node in weights:
        paths[node].CalculateWeight ( )
    # . Finish up.
    if target is None: result = paths
    else:              result = paths[target]
    return result

def DijkstraSingleSource ( graph, source, cutOff = None, target = None ):
    """Returns shortest paths and lengths in a weighted graph."""
    # . Initialization.
    distances = { source : 0 }
    paths     = { source : [ source ] }
    # . Do nothing if there are no edges or the source and target are identical.
    if ( len ( graph.edges ) > 0 ) and ( source is not target ):
        # . Check edge weights.
        if graph.MinimumEdgeWeight ( ) < 0: raise GraphError ( "Graph has negative edge weights." )
        # . Initialization.
        distances = {} # . Needed to avoid immediate exit.
        queue     = []
        seen      = { source : 0 }
        heapq.heappush ( queue, ( 0, source ) )
        # . Loop until the queue is empty.
        while len ( queue ) > 0:
            ( distance, node ) = heapq.heappop ( queue )
            if node in distances: continue
            distances[node] = distance
            if node is target: break
            for edge in graph.adjacentEdges.get ( node, [] ):
                other    = edge.Opposite ( node )
                distance = distances[node] + edge.weight
                if cutOff is not None:
                    if distance > cutOff: continue
                if other in distances:
                    if distance < distances[other]:
                        raise GraphError ( "Logic error: shorter path found." )
                elif ( other not in seen ) or ( distance < seen[other] ):
                    heapq.heappush ( queue, ( distance, other ) )
                    paths[other] = paths[node] + [other]
                    seen [other] = distance
    # . Finish up.
    return ( distances, paths )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
