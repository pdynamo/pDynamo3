"""Graph class."""

# . A class for representing undirected graphs.
# . Directed and multigraphs might need modification.

from  collections import defaultdict
from  pCore       import AttributableObject
from .GraphStatus import GraphError
from .Path        import Path

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class Graph ( AttributableObject ):
    """A basic graph class."""
 
    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "adjacentEdges" : None ,
                             "adjacentNodes" : None ,
                             "edges"         : None ,
                             "nodes"         : None } )

    def __len__ ( self ):
        """The number of nodes."""
        if self.nodes is None: return 0
        else:                  return len ( self.nodes )

    def _Initialize ( self ):
        """Initialization."""
        super ( Graph, self )._Initialize ( )
        self.Clear ( )

    def AddEdge ( self, edge ):
        """Add an edge to the graph."""
        if self.HasEdge ( edge ):
            raise GraphError ( "Edge already present in graph." )
        else:
            self.edges.append ( edge )
            self.adjacentNodes[edge.node1].add ( edge.node2 )
            self.adjacentNodes[edge.node2].add ( edge.node1 )
            self.adjacentEdges[edge.node1].add ( edge       )
            self.adjacentEdges[edge.node2].add ( edge       )

    def AddEdges ( self, edges ):
        """Add edges to the graph."""
        for edge in edges: self.AddEdge ( edge )

    def AddNode ( self, node ):
        """Add a node to the graph."""
        if node in self.adjacentNodes:
            raise GraphError ( "Node already present in graph." )
        else:
            self.nodes.append ( node )

    def AddNodes ( self, nodes ):
        """Add nodes to the graph."""
        for node in nodes: self.AddNode ( node )

    def Clear ( self ):
        """Clear all edges and nodes."""
        self.ClearEdges ( )
        self.nodes = []

    def ClearEdges ( self ):
        """Clear all edges."""
        self.adjacentEdges = defaultdict ( set )
        self.adjacentNodes = defaultdict ( set )
        self.edges         = []

    def GetEdge ( self, node1, node2 ):
        """Get the edge between two nodes."""
        edge = None
        for toTry in self.adjacentEdges[node1]:
            if toTry.Opposite ( node1 ) is node2:
                edge = toTry
                break
        return edge

    # . Check by end points not by edge identity.
    def HasEdge ( self, edge ):
        """Is the edge in the graph?"""
        return ( ( edge.node2 in self.adjacentNodes[edge.node1] ) and \
                 ( edge.node1 in self.adjacentNodes[edge.node2] ) )

    def MakePathFromNodes ( self, nodes ):
        """Generate a path from a sequence of nodes."""
        path = Path ( )
        for node in nodes:
            path.AddNode ( node )
        if len ( nodes ) > 1:
            tail = nodes[0]
            for head in nodes[1:]:
                found = False
                for edge in self.adjacentEdges[tail]:
                    if ( ( edge.node1 is tail ) and ( edge.node2 is head ) ) or \
                       ( ( edge.node2 is tail ) and ( edge.node1 is head ) ) :
                       found = True
                       path.AddEdge ( edge )
                       break
                if not found: raise GraphError ( "Unable to find edge between two putative path nodes." )
                tail = head
        return path

    def MakeSubgraph ( self, nodes, induced = True ):
        """Generate a subgraph with the selected nodes."""
        pruned = self.__class__ ( )
        if len ( nodes ) > 0:
            for node in nodes:
                if node in self.nodes: pruned.AddNode ( node )
                else: raise GraphError ( "Node not present in graph." )
            if induced:
                    for edge in self.edges:
                        if ( edge.node1 in nodes ) and ( edge.node2 in nodes ):
                            pruned.AddEdge ( edge )
        return pruned

    def MinimumEdgeWeight ( self ):
        """Return the minimum edge weight."""
        if len ( self.edges ) > 0: return min ( [ edge.weight for edge in self.edges ] )
        else:                      return 0.0

    def RemoveEdge ( self, edge ):
        """Remove an edge from the graph."""
        try:    self.edges.remove ( edge )
        except: raise GraphError ( "Edge not in graph." )
        iNode = edge.node1
        jNode = edge.node2
        self.adjacentEdges[iNode].remove ( edge  )
        self.adjacentEdges[jNode].remove ( edge  )
        self.adjacentNodes[iNode].remove ( jNode )
        self.adjacentNodes[jNode].remove ( iNode )

    def RemoveEdges ( self, edges ):
        """Remove edges from the graph."""
        for edge in edges: self.RemoveEdge ( edge )

    def RemoveNode ( self, node ):
        """Remove a node from the graph."""
        try:    self.nodes.remove ( node )
        except: raise GraphError ( "Node not in graph." )
        edges = self.adjacentEdges.pop ( node, set ( ) )
        for edge in edges:
            self.edges.remove ( edge )
            other = edge.Opposite ( node )
            self.adjacentEdges[other].remove ( edge )
        for other in self.adjacentNodes.pop ( node, set ( ) ):
            self.adjacentNodes[other].remove ( node )
        return edges

    def RemoveNodes ( self, nodes ):
        """Remove nodes from the graph."""
        edges = set ( )
        for node in nodes:
            edges.update ( self.RemoveNode ( node ) )
        return edges

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass

