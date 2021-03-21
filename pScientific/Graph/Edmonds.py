"""Find maximum cardinality matchings in general undirected graphs."""

#
# . Modified for pDynamo graph format from "CardinalityMatching.py" in the PADS library by D. Eppstein, UC Irvine.
# . The algorithm is Edmonds' blossom-contraction algorithm, as described in Galil's 1986 Computing Surveys paper.
#

#===================================================================================================================================
# . Utility class and function.
#===================================================================================================================================
class UnionFind:
    """Union-find data structure.

    Each unionFind instance X maintains a family of disjoint sets of
    hashable objects, supporting the following two methods:

    - X[item] returns a name for the set containing the given item.
      Each set is named by an arbitrarily-chosen one of its members; as
      long as the set remains unchanged it will keep the same name. If
      the item is not yet part of a set in X, a new singleton set is
      created for it.

    - X.union(item1, item2, ...) merges the sets containing each item
      into a single larger set.  If any item is not yet part of a set
      in X, it is added to X as one of the members of the merged set.
    """

    def __init__(self):
        """Create a new empty union-find structure."""
        self.weights = {}
        self.parents = {}

    def __getitem__(self, object):
        """Find and return the name of the set containing the object."""

        # check for previously unknown object
        if object not in self.parents:
            self.parents[object] = object
            self.weights[object] = 1
            return object

        # find path of objects leading to the root
        path = [object]
        root = self.parents[object]
        while root != path[-1]:
            path.append(root)
            root = self.parents[root]

        # compress the path and return
        for ancestor in path:
            self.parents[ancestor] = root
        return root
        
    def __iter__(self):
        """Iterate through all items ever found or unioned by this structure."""
        return iter(self.parents)

    def union(self, *objects):
        """Find the sets containing the objects and merge them all."""
        roots = [self[x] for x in objects]
        heaviest = max([(self.weights[r],r) for r in roots])[1]
        for r in roots:
            if r != heaviest:
                self.weights[heaviest] += self.weights[r]
                self.parents[r] = heaviest

def ArbitraryItem ( S ):
    """
    Select an arbitrary item from set or sequence S.
    Avoids bugs caused by directly calling iter(S).next() and
    mysteriously terminating loops in callers' code when S is empty.
    """
    try:                  return next ( iter ( S ) )
    except StopIteration: raise IndexError ( "No items to select." )

#===================================================================================================================================
# . Matching function.
#===================================================================================================================================
def EdmondsMaximumMatching ( G, initialMatching = None ):
    """Find a maximum cardinality matching in a graph G."""

    # Copy initial matching so we can use it nondestructively
    # and augment it greedily to reduce main loop iterations
    matching = GreedyMatching(G,initialMatching)

    def augment():
        """Search for a single augmenting path.
        Returns true if the matching size was increased, false otherwise.
        """

        # Data structures for augmenting path search:
        #
        # leader: union-find structure; the leader of a blossom is one
        # of its vertices (not necessarily topmost), and leader[v] always
        # points to the leader of the largest blossom containing v
        #
        # S: dictionary of blossoms at even levels of the structure tree.
        # Dictionary keys are names of blossoms (as returned by the union-find
        # data structure) and values are the structure tree parent of the blossom
        # (a T-node, or the top vertex if the blossom is a root of a structure tree).
        #
        # T: dictionary of vertices at odd levels of the structure tree.
        # Dictionary keys are the vertices; T[x] is a vertex with an unmatched
        # edge to x.  To find the parent in the structure tree, use leader[T[x]].
        #
        # unexplored: collection of unexplored vertices within blossoms of S
        #
        # base: if x was originally a T-vertex, but becomes part of a blossom,
        # base[t] will be the pair (v,w) at the base of the blossom, where v and t
        # are on the same side of the blossom and w is on the other side.

        leader = UnionFind()
        S = {}
        T = {}
        unexplored = []
        base = {}

        # Subroutines for augmenting path search.
        # Many of these are called only from one place, but are split out
        # as subroutines to improve modularization and readability.

        def blossom(v,w,a):
            """Create a new blossom from edge v-w with common ancestor a."""

            def findSide(v,w):
                path = [leader[v]]
                b = (v,w)   # new base for all T nodes found on the path
                while path[-1] != a:
                    tnode = S[path[-1]]
                    path.append(tnode)
                    base[tnode] = b
                    unexplored.append(tnode)
                    path.append(leader[T[tnode]])
                return path

            a = leader[a]   # sanity check
            path1,path2 = findSide(v,w), findSide(w,v)
            leader.union(*path1)
            leader.union(*path2)
            S[leader[a]] = S[a] # update structure tree

        topless = object()  # should be unequal to any graph vertex
        def alternatingPath(start, goal = topless):
            """Return sequence of vertices on alternating path from start to goal.
            The goal must be a T node along the path from the start to
            the root of the structure tree. If goal is omitted, we find
            an alternating path to the structure tree root.
            """
            path = []
            while 1:
                while start in T:
                    v, w = base[start]
                    vs = alternatingPath(v, start)
                    vs.reverse()
                    path += vs
                    start = w
                path.append(start)
                if start not in matching:
                    return path     # reached top of structure tree, done!
                tnode = matching[start]
                path.append(tnode)
                if tnode == goal:
                    return path     # finished recursive subpath
                start = T[tnode]

        def alternate(v):
            """Make v unmatched by alternating the path to the root of its structure tree."""
            path = alternatingPath(v)
            path.reverse()
            for i in range(0,len(path)-1,2):
                matching[path[i]] = path[i+1]
                matching[path[i+1]] = path[i]

        def addMatch(v, w):
            """Here with an S-S edge vw connecting vertices in different structure trees.
            Find the corresponding augmenting path and use it to augment the matching.
            """
            alternate(v)
            alternate(w)
            matching[v] = w
            matching[w] = v

        def ss(v,w):
            """Handle detection of an S-S edge in augmenting path search.
            Like augment(), returns true iff the matching size was increased.
            """

            if leader[v] == leader[w]:
                return False        # self-loop within blossom, ignore

            # parallel search up two branches of structure tree
            # until we find a common ancestor of v and w
            path1, head1 = {}, v
            path2, head2 = {}, w

            def step(path, head):
                head = leader[head]
                parent = leader[S[head]]
                if parent == head:
                    return head     # found root of structure tree
                path[head] = parent
                path[parent] = leader[T[parent]]
                return path[parent]

            while 1:
                head1 = step(path1, head1)
                head2 = step(path2, head2)

                if head1 == head2:
                    blossom(v, w, head1)
                    return False

                if leader[S[head1]] == head1 and leader[S[head2]] == head2:
                    addMatch(v, w)
                    return True

                if head1 in path2:
                    blossom(v, w, head1)
                    return False

                if head2 in path1:
                    blossom(v, w, head2)
                    return False

        # Start of main augmenting path search code.

        for v in G.nodes: # . MJF.
            if v not in matching:
                S[v] = v
                unexplored.append(v)

        current = 0     # index into unexplored, in FIFO order so we get short paths
        while current < len(unexplored):
            v = unexplored[current]
            current += 1

            for w in G.adjacentNodes[v]: # . MJF.
                if leader[w] in S:  # S-S edge: blossom or augmenting path
                    if ss(v,w):
                        return True

                elif w not in T:    # previously unexplored node, add as T-node
                    T[w] = v
                    u = matching[w]
                    if leader[u] not in S:
                        S[u] = w    # and add its match as an S-node
                        unexplored.append(u)

        return False    # ran out of graph without finding an augmenting path

    # augment the matching until it is maximum
    while augment():
        pass

    return matching

def GreedyMatching(G, initialMatching=None):
    """Near-linear-time greedy heuristic for creating high-cardinality matching.
    If there is any vertex with one unmatched neighbor, we match it.
    Otherwise, if there is a vertex with two unmatched neighbors, we contract
    it away and store the contraction on a stack for later matching.
    If neither of these two cases applies, we match an arbitrary edge.
    """

    # Copy initial matching so we can use it nondestructively
    matching = {}
    if initialMatching:
        for x in initialMatching:
            matching[x] = initialMatching[x]

    # Copy graph to new subgraph of available edges
    # Representation: nested dictionary rep->rep->pair
    # where the reps are representative vertices for merged clusters
    # and the pair is an unmatched original pair of vertices
    avail = {}
    has_edge = False
    for v in G.nodes: # . MJF.
        if v not in matching:
            avail[v] = {}
            for w in G.adjacentNodes[v]: # . MJF.
                if w not in matching:
                    avail[v][w] = (v,w)
                    has_edge = True
            if not avail[v]:
                del avail[v]
    if not has_edge:
        return matching

    # make sets of degree one and degree two vertices
    deg1 = {v for v in avail if len(avail[v]) == 1}
    deg2 = {v for v in avail if len(avail[v]) == 2}
    d2edges = []
    def updateDegree(v):
        """Cluster degree changed, update sets."""
        if v in deg1:
            deg1.remove(v)
        elif v in deg2:
            deg2.remove(v)
        if len(avail[v]) == 0:
            del avail[v]
        elif len(avail[v]) == 1:
            deg1.add(v)
        elif len(avail[v]) == 2:
            deg2.add(v)

    def addMatch(v,w):
        """Add edge connecting two given cluster reps, update avail."""
        p,q = avail[v][w]
        matching[p] = q
        matching[q] = p
        for x in avail[v].keys():
            if x != w:
                del avail[x][v]
                updateDegree(x)
        for x in avail[w].keys():
            if x != v:
                del avail[x][w]
                updateDegree(x)
        avail[v] = avail[w] = {}
        updateDegree(v)
        updateDegree(w)

    def contract(v):
        """Handle degree two vertex."""
        u,w = avail[v]  # find reps for two neighbors
        d2edges.extend([avail[v][u],avail[v][w]])
        del avail[u][v]
        del avail[w][v]
        if len(avail[u]) > len(avail[w]):
            u,w = w,u   # swap to preserve near-linear time bound
        for x in avail[u].keys():
            del avail[x][u]
            if x in avail[w]:
                updateDegree(x)
            elif x != w:
                avail[x][w] = avail[w][x] = avail[u][x]
        avail[u] = avail[v] = {}
        updateDegree(u)
        updateDegree(v)
        updateDegree(w)

    # loop adding edges or contracting deg2 clusters
    while avail:
        if deg1:
            v = ArbitraryItem(deg1)
            w = ArbitraryItem(avail[v])
            addMatch(v,w)
        elif deg2:
            v = ArbitraryItem(deg2)
            contract(v)
        else:
            v = ArbitraryItem(avail)
            w = ArbitraryItem(avail[v])
            addMatch(v,w)

    # at this point the edges listed in d2edges form a matchable tree
    # repeat the degree one part of the algorithm only on those edges
    avail = {}
    d2edges = [(u,v) for u,v in d2edges if u not in matching and v not in matching]
    for u,v in d2edges:
        avail[u] = {}
        avail[v] = {}
    for u,v in d2edges:
        avail[u][v] = avail[v][u] = (u,v)
    deg1 = {v for v in avail if len(avail[v]) == 1}
    while deg1:
        v = ArbitraryItem(deg1)
        w = ArbitraryItem(avail[v])
        addMatch(v,w)

    # . Finish up.
    return matching

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
