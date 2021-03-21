"""A class to represent the head of a path.""" # . walk?

from  pCore       import AttributableObject
from .GraphStatus import GraphError
from .Path        import Path

# . Be careful with comparisons as they can be done by length, by weight or by node/edge identity.

# . Check for same node occurring multiple times?

# . Should really be called WalkHead as no checks made for duplicate nodes?

#===================================================================================================================================
# . Function.
#===================================================================================================================================
class PathHead ( AttributableObject ):
    """A path head class."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "headEdge" : None ,
                             "headNode" : None ,
                             "index"    :   -1 ,
                             "tail"     : None ,
                             "weight"   :  0.0 } )

    def __len__ ( self ):
        """Length of the path."""
        length   = 1
        headEdge = self.headEdge
        tail     = self.tail
        while ( headEdge is not None ) and ( tail is not None ):
            length  += 1
            headEdge = tail.headEdge
            tail     = tail.tail
        return length

    def _CheckOptions ( self ):
        """Check options."""
        super ( PathHead, self )._CheckOptions ( )
        self.VerifyTail ( )

    def CalculateWeight ( self ):
        """Get the weight."""
        weight   = 0.0
        headEdge = self.headEdge
        tail     = self.tail
        while ( headEdge is not None ) and ( tail is not None ):
            weight  += headEdge.weight
            headEdge = tail.headEdge
            tail     = tail.tail
        self.weight = weight

    def Edges ( self ):
        """Return the edges in the path (tail->head)."""
        edges = []
        edge  = self.headEdge
        tail  = self.tail
        while ( edge is not None ) and ( tail is not None ):
            edges.append ( edge )
            edge = tail.headEdge
            tail = tail.tail
        edges.reverse ( )
        return edges

    def IsIdenticalTo ( self, other ):
        """Are paths equal?"""
        result = False
        if self is other:
            result = True
        elif len ( self ) == len ( other ):
            result    = True
            selfTail  = self
            otherTail = other
            while ( selfTail is not None ) and ( otherTail is not None ):
                if ( selfTail.headEdge is otherTail.headEdge ) and ( selfTail.headNode is otherTail.headNode ):
                    selfTail  = selfTail.tail
                    otherTail = otherTail.tail
                else:
                    result = False
                    break
        return result

    def IsSimplePath ( self ):
        """Check to see if the path is simple."""
        nodes = set ( [ self.headNode ] )
        tail  = self.tail
        while tail is not None:
            nodes.add ( tail.headNode )
            tail = tail.tail
        return ( len ( nodes ) == len ( self ) )

    def Nodes ( self ):
        """Return the nodes in the path (tail->head)."""
        nodes = [ self.headNode ]
        tail  = self.tail
        while tail is not None:
            nodes.append ( tail.headNode )
            tail = tail.tail
        nodes.reverse ( )
        return nodes

    def PushHead ( self, edge ):
        """Return a new path with |edge| added."""
        newHeadNode = edge.Opposite ( self.headNode )
#        if newHeadNode is None: raise GraphError ( "Edge does not connect to head node." )
        newPath = PathHead ( headEdge = edge, headNode = newHeadNode, tail = self )
        return newPath

    def Split ( self, splitHeadNode ):
        """Split the path at the head node of a path element."""
        splitEdge = None
        splitNode = None
        splitTail = None
        if self.headNode is splitHeadNode:
            splitTail = self
        else:
            edge = self.headEdge
            node = self.headNode
            tail = self.tail
            while tail is not None:
                if tail.headNode is splitHeadNode:
                    splitEdge = edge
                    splitNode = node
                    splitTail = tail
                    break
                edge = tail.headEdge
                node = tail.headNode
                tail = tail.tail
        return ( splitNode, splitEdge, splitTail )

# . Convert to walk?
    def ToPath ( self ):
        """Convert to a path."""
        path     = Path ( source = self.headNode )
        headEdge = self.headEdge
        tail     = self.tail
        while ( headEdge is not None ) and ( tail is not None ):
            path.Prepend ( headEdge )
            headEdge = tail.headEdge
            tail     = tail.tail
        return path

    def UpdateTail ( self, headEdge = None, tail = None ):
        """Update the tail."""
        self.headEdge = headEdge
        self.tail     = tail
        self.VerifyTail ( )

    def VerifyTail ( self ):
        """Verify the tail."""
        isOK = True
        if self.headEdge is None:
            isOK = ( self.tail is None )
        else:
            isOK = ( self.tail is not None ) and ( self.tail.headNode is self.headEdge.Opposite ( self.headNode ) )
        if not isOK: raise GraphError ( "Invalid tail for path." )
        self.CalculateWeight ( )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
