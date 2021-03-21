"""Graph edge."""

from  pCore       import AttributableObject
from .GraphStatus import GraphError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class Edge ( AttributableObject ):
    """A graph edge."""
 
    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "node1"  : None ,
                             "node2"  : None ,
                             "weight" :  1.0 } )

    def Opposite ( self, node ):
        """Return the other node for the edge."""
        if   node is self.node1: return self.node2
        elif node is self.node2: return self.node1
        else: raise GraphError ( "Node not associated with edge." )

    def SourceTargetPairs ( self ):
        """Return a tuple of source/target pairs."""
        return ( ( self.node1, self.node2 ), ( self.node2, self.node1 ) )

    @classmethod
    def WithNodes ( selfClass, node1, node2, **options ):
        """Constructor with nodes and other options."""
        options = dict ( options )
        options["node1"] = node1
        options["node2"] = node2
        return selfClass.WithOptions ( **options )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
