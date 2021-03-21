"""Shortest paths using Bellman-Ford algorithm."""

from .GraphStatus import GraphError
from .PathHead    import PathHead

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def BellmanFordShortestPaths ( graph, source ):
    """Calculate the shortest path distance between the source node and all other nodes in the graph using Bellman-Ford's algorithm."""
    # . Initialization.
    paths = {}
    paths[source] = PathHead ( headNode = source )

    # . Iterate and relax.
    for i in range ( len ( graph.nodes ) - 1 ):
        for edge in graph.edges:
            for ( source, target ) in edge.SourceTargetPairs ( ):
                sourcePath = paths.get ( source, None )
                if sourcePath is not None:
                    targetPath = paths.get ( target, None )
                    if targetPath is None:
                        paths[target] = PathHead ( headEdge = edge, headNode = target, tail = sourcePath )
                    else:
                        if targetPath.weight > sourcePath.weight + edge.weight:
                            targetPath.UpdateTail ( headEdge = edge, tail = sourcePath )

    # . Detect negative weight cycles.
    for edge in graph.edges:
        for ( source, target ) in edge.SourceTargetPairs ( ):
            if paths[target].weight > ( paths[source].weight + edge.weight ):
                raise GraphError ( "Graph has negative weight cycle on edge {:s}.".format ( edge ) )

    # . Finish up.
    return paths

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
