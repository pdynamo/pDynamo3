"""Enumerate all the cycles of an undirected graph."""

# . From Hanser T, Jauffret P, Kaufmann G. J Chem Inf Comp Sci 36, 1146-1152, 1996.
# . With hints from John May's thesis 2014.

# . PathGraph and PathEdge are classes constructed to work with the algorithm. They are incomplete for general use.

from  collections           import defaultdict
from  itertools             import combinations
from .BiconnectedComponents import BiconnectedComponents
from .Edge                  import Edge
from .Graph                 import Graph
from .GraphStatus           import GraphError

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
_DefaultMaximumCycleSize = None # . Maximum cycle size.
_DefaultMaximumDegree    = 1000 # . Maximum reduced node degree.

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PathEdge ( Edge ):
    """A path edge."""

    _attributable = dict ( Edge._attributable )
    _attributable.update ( { "nodes" : list , # . The intermediate nodes along the path.
                             "_path" : None } )

    def __len__ ( self ):
        return ( len ( self.nodes ) + 2 )

    def IsDisjoint ( self, other ):
        """Are the paths disjoint apart from endpoints."""
        return ( len ( set ( self.nodes + other.nodes ) ) == ( len ( self.nodes ) + len ( other.nodes ) ) )

    @property
    def path ( self ):
        """Return the path."""
        if self._path is None:
            self._path = [ self.node1 ] + self.nodes + [ self.node2 ]
        return self._path

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PathGraph ( Graph ):
    """A path graph."""

    _attributable = dict ( Graph._attributable )
    _attributable.update (  { "degree"           : 0                        ,
                              "maximumCycleSize" : _DefaultMaximumCycleSize ,
                              "maximumDegree"    : _DefaultMaximumDegree    ,
                              "ordering"         : None                     } )

    def AddEdge ( self, edge ):
        """Add an edge to the graph."""
        self.edges.append ( edge )
        self.adjacentNodes[edge.node1].add ( edge.node2 ) # . order n1 < n2 so only n1 does indexing.
        self.adjacentEdges[edge.node1].add ( edge       )

    @classmethod
    def FromSubGraph ( selfClass, graph, nodes, maximumCycleSize = _DefaultMaximumCycleSize ,
                                                maximumDegree    = _DefaultMaximumDegree    ):
        """Constructor from a graph and a subset of nodes."""
        self = selfClass ( )
        if maximumCycleSize is None: self.maximumCycleSize = max ( len ( nodes ) + 1, 3 )
        else:                        self.maximumCycleSize = maximumCycleSize
        self.maximumDegree    = maximumDegree
        if len ( nodes ) > 0:
            # . Add nodes by degree (low to high).
            degrees = defaultdict ( int )
            edges   = [ edge for edge in graph.edges if ( edge.node1 in nodes ) and ( edge.node2 in nodes ) ]
            for edge in edges:
                degrees[edge.node1] += 1
                degrees[edge.node2] += 1 
            work       = sorted ( [ ( degrees[node], node ) for node in nodes ] )
            ordering   = { node : order for ( order, ( _, node ) ) in enumerate ( work ) }
            self.nodes = [ node for ( _, node ) in work ]
            for edge in edges:
                n1 = edge.node1
                n2 = edge.node2
                if ordering[n1] < ordering[n2]: self.AddEdge ( PathEdge.WithNodes ( n1, n2 ) )
                else:                           self.AddEdge ( PathEdge.WithNodes ( n2, n1 ) )
            self.ordering = ordering
        return self

    def Reduce ( self ):
        """Reduce the graph by removing all the nodes."""
        # . Loop over nodes - low to high priority.
        cycles = []
        while len ( self.nodes ) > 0:
            node        = self.nodes.pop ( 0 )
            edges       = self.adjacentEdges.pop ( node, [] )
            degree      = len ( edges )
            self.degree = max ( self.degree, degree ) # . Information only.
            if degree <= self.maximumDegree:
                # . Find edge order - use order here as path edges cannot be compared directly (as they may have the same node2).
                lEdges    = list ( edges )
                edgeOrder = sorted ( [ ( self.ordering[edge.node2], order ) for ( order, edge ) in enumerate ( lEdges ) ] )
                edges     = [ lEdges[order] for ( _, order ) in edgeOrder ]
                # . Loop over pairs of edges emanating from the node with e1 < e2.
                for i in range ( degree - 1 ):
                    edge1  = edges[i]
                    limit  = self.maximumCycleSize + 1 - len ( edge1 )
                    n1     = edge1.node2
                    nodes1 = edge1.nodes[::-1] + [ node ] # . Reversed.
                    for j in range ( i+1, degree ):
                        edge2 = edges[j]
                        # . Accept the new path if the intermediate nodes in the edge paths are unique, and the new path is not too long.
                        if edge1.IsDisjoint ( edge2 ) and ( len ( edge2 ) <= limit ):
                            n2 = edge2.node2
                            if n1 is n2: cycles.append ( [ n1 ] + nodes1 + edge2.nodes ) # . Cycle not closed.
                            else: self.AddEdge ( PathEdge.WithNodes ( n1, n2, nodes = nodes1 + edge2.nodes ) )
        # . Finish up.
        isOK = ( self.degree <= self.maximumDegree ) #  . A flag for incomplete searching.
        self.Clear ( )
        return ( cycles, isOK )

#===================================================================================================================================
# . Function.
#===================================================================================================================================
# . Could use checks here for simple cases (e.g. no branching).
def HanserAllCycles ( graph, biconnectedComponents = None                     ,
                             maximumCycleSize      = _DefaultMaximumCycleSize ,
                             maximumDegree         = _DefaultMaximumDegree    ):
    """Calculate the relevant cycles of an undirected graph."""
    cycleSets = []
    isOK      = True
    if biconnectedComponents is None:
        biconnectedComponents = BiconnectedComponents ( graph )
    for component in biconnectedComponents:
        pathGraph           = PathGraph.FromSubGraph ( graph, component, maximumCycleSize = maximumCycleSize, maximumDegree = maximumDegree )
        ( cycles, localOK ) = pathGraph.Reduce ( )
        #print ( "\nHAC> Maximum Reduced Node Degree = {:d}.\n".format ( pathGraph.degree ) )
        if len ( cycles ) > 0: cycleSets.append ( cycles )
        isOK = isOK and localOK
    return ( cycleSets, isOK )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
