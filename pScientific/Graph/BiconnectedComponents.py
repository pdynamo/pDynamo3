"""Calculate the biconnected components of an undirected graph."""

# . Non-recursive algorithm taken from networkX.

from itertools import chain

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
# . By default dyads are removed.
def BiconnectedComponents ( graph, removeDyads = True ):
    """Returns the biconnected components of a graph."""
    components = [ set ( chain.from_iterable ( component ) ) for component in _DepthFirstSearch ( graph, components = True ) ]
    if removeDyads: return [ component for component in components if len ( component ) > 2 ]
    else:           return components

def _DepthFirstSearch ( graph, components = True ):
    """Depth-first search algorithm to generate articulation points and biconnected components."""
    visited = set ( )
    for start in graph.nodes:
        if start in visited: continue
        discovery     = { start: 0 }  # . Time of first discovery of node during search.
        low           = { start: 0 }
        root_children = 0
        visited.add ( start )
        edge_stack = []
        stack = [ ( start, start, iter ( graph.adjacentNodes[start] ) ) ]
        while len ( stack ) > 0:
            ( grandparent, parent, children ) = stack[-1]
            try:
                child = next ( children )
                if grandparent is child: continue
                if child in visited:
                    if discovery[child] <= discovery[parent]:  # . Back edge.
                        low[parent] = min ( low[parent], discovery[child] )
                        if components: edge_stack.append ( ( parent, child ) )
                else:
                    low[child] = discovery[child] = len ( discovery )
                    visited.add ( child )
                    stack.append ( ( parent, child, iter ( graph.adjacentNodes[child] ) ) )
                    if components: edge_stack.append ( ( parent, child ) )
            except StopIteration:
                stack.pop ( )
                if len ( stack ) > 1:
                    if low[parent] >= discovery[grandparent]:
                        if components:
                            ind = edge_stack.index ( ( grandparent, parent ) )
                            #print ( "Y1>", edge_stack[ind:], len ( edge_stack[ind:] ) )
                            yield edge_stack[ind:]
                            edge_stack = edge_stack[:ind]
                        else:
                            yield grandparent
                    low[grandparent] = min ( low[parent], low[grandparent] )
                elif len ( stack ) > 0:  # . Length 1 so grandparent is root.
                    root_children += 1
                    if components:
                        ind = edge_stack.index ( ( grandparent, parent ) )
                        #print ( "Y2>", edge_stack[ind:], len ( edge_stack[ind:] ) )
                        yield edge_stack[ind:]
        if not components:
            # . Root node is articulation point if it has more than 1 child.
            if root_children > 1:
                yield start

##===================================================================================================================================
## . Original algorithm.
##===================================================================================================================================
#
# . The original algorithm is recursive and will reach Python's recursion limit, even for relatively small graphs.
#
#def BiconnectedComponents ( graph ):
#    """Calculate the biconnected components of an undirected graph."""
#    # . Find the biconnected components.
#    components  = []
#    count1      =  0
#    count2      =  0
#    current     = []
#    depth       = {}
#    low         = {}
#    predecessor = {}
#    for node in graph.nodes:
#        if node not in depth:
#            count1 += 1
#            depth[node] = count1
#            current.append ( node )
#            ( count1, count2 ) = _DepthFirstSearch ( graph, node, components, depth, low, predecessor, current, count1, count2 )
#    # . Sort into components - largest to smallest.
#    components.sort ( key = len, reverse = True )
#    return components
#
#def _DepthFirstSearch ( graph, node, components, depth, low, predecessor, current, count1, count2 ):
#    """Depth first search."""
#    low[node] = depth[node]
#    for opposite in graph.adjacentNodes.get ( node, [] ):
#        if opposite not in depth:
#            count1 += 1
#            depth[opposite] = count1
#            current.append ( opposite )
#            predecessor[opposite] = node
#            ( count1, count2 ) = _DepthFirstSearch ( graph, opposite, components, depth, low, predecessor, current, count1, count2 )
#            low[node] = min ( low[node], low[opposite] )
#        else:
#            low[node] = min ( low[node], depth[opposite] )
#    if ( node in predecessor ) and ( low[node] == depth[predecessor[node]] ):
#        nodes = set ( )
#        while True:
#            next = current.pop ( -1 )
#            for opposite in graph.adjacentNodes.get ( next, [] ):
#                if depth[next] > depth[opposite]:
#                    nodes.add ( next     )
#                    nodes.add ( opposite )
#            if next is node: break
#        if len ( nodes ) > 2:
#            components.append ( nodes )
#            count2 += 1
#    return ( count1, count2 )
#
#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass

