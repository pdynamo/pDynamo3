"""Determine a minimal cycle basis for an undirected graph."""

# . From Paton K. Comm ACM 12, 514-518, 1969.

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def PatonMinimalCycleBasis ( graph ):
    """Determine a minimal cycle basis for an undirected graph."""
    cycleSets = []
    nodes     = set ( graph.nodes )
    root      = graph.nodes[0]
    while len ( nodes ) > 0:
        if root is None: root = nodes.pop ( )
        cycles      = []
        stack       = [ root ]
        predecessor = { root: root    }
        used        = { root: set ( ) }
        while len ( stack ) > 0:
            current = stack.pop ( )
            for next in graph.adjacentNodes[current]:
                if next not in used:
                    predecessor[next] =         current
                    used       [next] = set ( [ current ] )
                    stack.append ( next )
                elif next is current:
                    cycles.append ( [ current ] )
                elif next not in used[current]:
                    cycle = [ next, current ]
                    last  = predecessor[current]
                    while last not in used[next]:
                        cycle.append ( last )
                        last = predecessor[last]
                    cycle.append  ( last  )
                    cycles.append ( cycle )
                    used[next].add ( current )
        nodes -= set ( predecessor )
        root   = None
        if len ( cycles ) > 0: cycleSets.append ( cycles )
    return cycleSets

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass

