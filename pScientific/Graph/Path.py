"""Path classes and functions."""

"""Definitions:

walk        - any consecutive sequence of nodes and edges.
closed walk - first and last nodes are the same.
open walk   - first and last nodes are different.

trail       - a walk in which all edges are distinct.

path        - a walk in which no vertices occur twice.
cycle       - a closed path.

"""

from  pCore       import AttributableObject
from .GraphStatus import GraphError

# . Be careful with comparisons as they can be done by length, by weight or by node/edge identity.

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class Path ( AttributableObject ):
    """An open or closed path in a graph."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "edges"             : list ,
                             "nodes"             : list ,
                             "weight"            :  0.0 ,
                             "_weightComparison" : True } )

# . Length = number of edges or nodes?

    def __contains__ ( self, item ):
        return ( item in self.edges ) or ( item in self.nodes )

    def __len__ ( self ):
        """Length."""
        return len ( self.edges )

    def AddEdge ( self, edge ):
        """Add an edge."""
        if   edge.node1 is self.nodes[-1]: newNode = edge.node2
        elif edge.node2 is self.nodes[-1]: newNode = edge.node1
        else: newNode = None
        if newNode is None: raise GraphError ( "Edge does not connect to last node." )
        else:               self.AddNode ( newNode )
        self.edges.append ( edge )
        self.weight += edge.weight

    def AddNode ( self, node ):
        """Add a node."""
        if node not in self.nodes: self.nodes.append ( node )
        else: raise GraphError ( "Node already in path." )

    def Duplicate ( self ):
        """Duplicate the current path but retaining node and edge identity."""
        copy = Path ( )
        copy.edges  = list ( self.edges )
        copy.nodes  = list ( self.nodes )
        copy.weight = self.weight
        return copy

#???

    def IsClosed ( self ):
        """Is the path closed (i.e. a cycle)."""
        return ( len ( self.nodes ) > 2 ) and ( self.nodes[0] is self.nodes[-1] )

    # . Properties.
    @property
    def length ( self ): return len ( self )

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def EdgeVectorToPath ( graph, vector ):
    """Return a list of booleans in graph-edge order indicating which edges occur in the path."""
    # . Initialization.
    path = []
    # . Find the connection map.
    connections = {}
    for ( e, isPresent ) in enumerate ( vector ):
        if isPresent:
            node1 = graph.edges[e].node1
            node2 = graph.edges[e].node2
            for ( node, other ) in ( ( node1, node2 ), ( node2, node1 ) ):
                local = connections.get ( node, set ( ) )
                local.add ( other )
                connections[node] = local
    # . Check the consistency of the map.
    endNode            = None
    maximumConnections = 0
    numberConnections  = 0
    numberOnes         = 0
    for ( node, values ) in connections.items ( ):
        n = len ( values )
        if n == 1:
            endNode     = node
            numberOnes += 1
        maximumConnections = max ( maximumConnections, n )
        numberConnections += n
    # . No nodes.
    if len ( connections ) == 0:
        pass
    # . A closed or open path.
    elif ( maximumConnections == 2 ) and ( numberOnes in ( 0, 2 ) ):
        # . Create the path.
        if numberOnes == 0: current = list ( connections.keys ( ) )[0]
        else:               current = endNode
        while True:
            path.append ( current )
            local = connections[current]
            if len ( local ) == 0: break
            next  = local.pop ( )
            connections[next].remove ( current )
            current = next
    # . A problem.
    else: raise GraphError ( "Cannot create a valid path from the edge vector." )
    # . Finish up.
    return path

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def PathToEdgeVector ( graph, path, closePath = False ):
    """Return a list of booleans in graph-edge order indicating which edges occur in the path."""
    # . Check whether to close the path.
    extraNode = []
    if closePath and ( len ( path ) > 0 ) and ( path[0] is not path[-1] ): extraNode = [ path[0] ]
    # . Create the vector.
    vector = [ False for i in range ( len ( graph.edges ) ) ]
    tail   = path[0]
    for ( i, head ) in enumerate ( ( path + extraNode )[1:] ):
        found = False
        for ( e, edge ) in enumerate ( graph.edges ):
            if ( ( edge.node1 is head ) and ( edge.node2 is tail ) ) or ( ( edge.node2 is head ) and ( edge.node1 is tail ) ):
                found     = True
                vector[e] = True
                break
        if not found: raise GraphError ( "A path contains two nodes without an accompanying edge." )
        tail = head
    return vector

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
