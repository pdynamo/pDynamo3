"""A "connectivity" is a set of atoms and bonds and related information such as angles, dihedrals, molecules and ringsets."""

from  collections           import defaultdict
from  itertools             import chain                         , \
                                   combinations
from  pCore                 import logFile                       , \
                                   LogFileActive                 , \
                                   SelfPairList
from  pScientific.Graph     import BiconnectedComponents         , \
                                   ConnectedComponents           , \
                                   Edge                          , \
                                   Graph                         , \
                                   Node
from  pScientific.Geometry3 import SelfPairList_FromCoordinates3
from .Atom                  import Atom
from .Bond                  import Bond
from .ConnectivityUtilities import DetermineAromaticity          , \
                                   MakeAtomAttributes

# . A connectivity should, in principle, store a complete and valid Kekule structure.
# . Atoms and bonds belonging to aromatic rings are flagged only.

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class Connectivity ( Graph ):
    """A class to store connectivity data."""

    # . Attributes that can be derived from atom and bond information.
    _derivable = ( "_angleIndices"    ,
                   "_angles"          ,
                   "_dihedralIndices" ,
                   "_dihedrals"       ,
                   "_isolateIndices"  ,
                   "_isolates"        ,
                   "_nodeIndices"     ,
                   "_ringSetIndices"  ,
                   "_ringSets"        )

    # . Default attributes.
    _attributable = dict ( Graph._attributable )
    _attributable.update ( { key : None for key in _derivable } )

    # . Method aliases (for ease).
    # . Care required here to ensure these methods do not rely on Connectivity attributes (just those of Graph).
    AddAtom = Graph.AddNode
    AddBond = Graph.AddEdge

    #-------------------------------------------------------------------------------------------------------------------------------
    # . Index methods.
    #-------------------------------------------------------------------------------------------------------------------------------
    def _Indices123And13AsSets ( self ):
        """Get a list of 123 and 13 indices."""
        s13 = set ( [ ( i, k ) for ( i, j, k ) in self.angleIndices ] ) # . i > k.
        s12 = set ( self.bondIndices ) # . i > j.
        return ( s13.union ( s12 ), s13.difference ( s12 ) )

    def _Indices1234And14 ( self ):
        """Make 1234 and 14 indices."""
        # . Ensure there are no 1-2 or 1-3 interactions in the 1-4 list.
        # . i > j.
        ( s123, _ ) = self._Indices123And13AsSets ( )
        s14p = set ( [ ( max ( i, l ), min ( i, l ) ) for ( i, j, k, l ) in self.dihedralIndices ] )
        return ( sorted ( s123.union ( s14p ) ), sorted ( s14p.difference ( s123 ) ) )

    def _IndicesAngle ( self ):
        """Make angle indices."""
        # . Sorted by j with i > k.
        mapping = self.nodeIndices
        indices = []
        for ( iNode, jNode, kNode ) in self.angles:
            ( i, j, k ) = ( mapping[iNode], mapping[jNode], mapping[kNode] )
            if i > k: indices.append ( ( j, i, k ) )
            else:     indices.append ( ( j, k, i ) )
        return [ ( i, j, k ) for ( j, i, k ) in sorted ( indices ) ]

    def _IndicesBond ( self ):
        """Make bond indices.""" # . Equivalent to a list of 1-2 indices.
        # . Sorted with i > j.
        mapping = self.nodeIndices
        indices = []
        for edge in self.edges:
            ( i, j ) = (  mapping[edge.node1], mapping[edge.node2] )
            if i > j: indices.append ( ( i, j ) )
            else    : indices.append ( ( j, i ) )
        return sorted ( indices )

    def _IndicesDihedral ( self ):
        """Make dihedral indices."""
        # . Sorted with j > k.
        mapping = self.nodeIndices
        indices = []
        for ( iNode, jNode, kNode, lNode ) in self.dihedrals:
            ( i, j, k, l ) = ( mapping[iNode], mapping[jNode], mapping[kNode], mapping[lNode] )
            if j > k: indices.append ( ( j, k, i, l ) )
            else:     indices.append ( ( k, j, l, i ) )
        return [ ( i, j, k, l ) for ( j, k, i, l ) in sorted ( indices ) ]

    def _IndicesIsolate ( self ):
        """Make isolate indices."""
        # . A selection container.
        bondPairList     = SelfPairList.FromIndexPairs ( list ( chain ( * list ( self.bondIndices ) ) ) )
        indices          = bondPairList.GetConnectedComponents ( upperBound = len ( self.atoms ) )
        indices.itemName = "Isolate"
        return indices

    def _IndicesNeighbor3 ( self ):
        """Get the indices of nodes with three neighbors."""
        mapping = self.nodeIndices
        indices = []
        for node in self.nodes:
            neighbors = self.adjacentNodes[node]
            if len ( neighbors ) == 3:
                indices.append ( [ mapping[node] ] + sorted ( [ mapping[n] for n in neighbors ] ) ) 
        return sorted ( indices )

    def _IndicesNode ( self ):
        """Make node indices."""
        return { node : order for ( order, node ) in enumerate ( self.nodes ) }

    def _IndicesRingSet ( self ):
        """Make ring set indices."""
        mapping = self.nodeIndices
        indices = []
        for ringSet in self.ringSets:
            indices.append ( sorted ( [ mapping[node] for node in ringSet ] ) )
        return sorted ( indices )

    #-------------------------------------------------------------------------------------------------------------------------------
    # . IC methods.
    #-------------------------------------------------------------------------------------------------------------------------------
    def _MakeAngles ( self ):
        """Make angles."""
        items = []
        for jNode in self.nodes:
            for ( iNode, kNode ) in combinations ( self.adjacentNodes[jNode], 2 ):
                items.append ( ( iNode, jNode, kNode ) )
        return items

    def _MakeDihedrals ( self ):
        """Make angles."""
        items = []
        for edge in self.edges:
            jNode = edge.node1
            kNode = edge.node2
            for iNode in self.adjacentNodes[jNode]:
                if iNode is not kNode:
                    for lNode in self.adjacentNodes[kNode]:
                        if ( lNode is not iNode ) and ( lNode is not jNode ):
                            items.append ( ( iNode, jNode, kNode, lNode ) )
        return items

    #-------------------------------------------------------------------------------------------------------------------------------
    # . Other methods.
    #-------------------------------------------------------------------------------------------------------------------------------
    def BondsFromCoordinates ( self, coordinates3, radii, safety ):
        """Estimate bonds from coordinates using a distance search."""
        if coordinates3.rows != len ( self.nodes ): raise TypeError ( "Coordinates3 argument of invalid extent." )
        self.ClearDerivableAttributes ( )
        self.ClearEdges               ( )
        pairList = SelfPairList_FromCoordinates3 ( coordinates3, radii = radii, safety = safety )
        for ( i, j ) in pairList:
            self.AddEdge ( Bond.WithNodes ( self.nodes[i], self.nodes[j] ) )

    def BondsFromIterable ( self, iterable ):
        """Add bonds from an iterable."""
        self.ClearDerivableAttributes ( )
        if isinstance ( iterable, dict ): # . Dictionary -> list.
            items = [ ( n1, n2, key ) for ( key, values ) in iterable.items ( ) for ( n1, n2 ) in values  ]
        else:
            try:    items = iter ( iterable )
            except: TypeError ( "Argument must be an iterable sequence." )
        for item in items: self.AddEdge ( Bond.FromIterable ( item, self.nodes ) )

    def Clear ( self ):
        """Clear all attributes."""
        self.ClearDerivableAttributes ( )
        super ( Connectivity, self ).Clear ( )

    def ClearDerivableAttributes ( self ):
        """Clear derivable attributes."""
        for attribute in self.__class__._derivable: self.__dict__.pop ( attribute, None )

    def CompleteConnectivity ( self ):
        """Complete the connectivity."""
        MakeAtomAttributes ( self, "connections", "hydrogens", "valence" )
        DetermineAromaticity ( self )

    @classmethod
    def FromAtoms ( selfClass, atoms, bonds = None ):
        """Constructor from an atom list and, optionally, bonds."""
        self = selfClass ( )
        self.AddNodes ( atoms )
        if bonds is not None: self.BondsFromIterable ( bonds )
        self.CompleteConnectivity ( )
        return self

    def IdentifyBoundaryAtoms ( self, selection, results ):
        """Fill the dictionary |results| boundary atoms, and their out-of-selection and in-selection partners."""
        if len ( selection ) > 0:
            nodeToIndex   = self.nodeIndices
            nodes         = [ self.nodes[s] for s in selection ]
            boundaryAtoms = set ( )
            for node in nodes: boundaryAtoms.update ( self.adjacentNodes[node] )
            boundaryAtoms = boundaryAtoms.difference ( nodes )
            if len ( boundaryAtoms ) > 0:
                for node in boundaryAtoms:
                    iNode = nodeToIndex[node]
                    data  = results.get ( iNode, [ set ( ), set ( ), set ( ) ] )
                    for n in self.adjacentNodes[node]:
                        i = nodeToIndex[n]
                        if   ( n in boundaryAtoms ): data[0].add ( i ) # . Boundary atom partners.
                        elif ( n not in nodes     ): data[1].add ( i ) # . Out-of-selection partners.
                        else:                        data[2].add ( i ) # . In-selection partners.
                    results[iNode] = data

    def Kekulize ( self ):
        """Kekulize the bonds of a connectivity."""
        for bond in self.edges:
            bond.isAromatic = False

    def Make1234And14PairLists ( self ):
        """Make 1234 and 14 pairlists."""
        ( i1234, i14 ) = self._Indices1234And14 ( )
        return ( SelfPairList.FromIndexPairs ( list ( chain ( * list ( i1234 ) ) ) ) ,
                 SelfPairList.FromIndexPairs ( list ( chain ( * list ( i14   ) ) ) ) )

    @classmethod
    def Merge ( selfClass, items, information = {} ):
        """Merging."""
        # . Atoms.
        atoms       = information.get ( "Atom Container", None )
        merged      = selfClass.FromAtoms ( atoms )
        # . Bonds.
        atomMapping = information.get ( "Atom Mapping", None )
        for ( item, mapping ) in zip ( items, atomMapping ):
            oldBonds = getattr ( item, "bonds", None )
            if oldBonds is not None:
                for old in oldBonds:
                    merged.AddEdge ( Bond.WithNodes ( mapping[old.node1], mapping[old.node2], type = old.type ) )
        merged.CompleteConnectivity ( )
        return merged

    def Prune ( self, selection, information = {} ):
        """Pruning."""
        # . Atoms.
        atoms  = information.get ( "Atom Container", None )
        pruned = self.__class__.FromAtoms ( atoms )
        # . Bonds.
        if len ( self.edges ) > 0:
            mapping = information.get ( "Atom Mapping", None )
            nodes   = [ self.nodes[s] for s in selection ]
            for old in self.edges:
                if ( old.node1 in nodes ) and ( old.node2 in nodes ):
                    pruned.AddEdge ( Bond.WithNodes ( mapping[old.node1], mapping[old.node2], type = old.type ) )
        pruned.CompleteConnectivity ( )
        return pruned

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            log.SummaryOfItems ( self.SummaryItems ( ), title = "Connectivity Summary" )

    def SummaryItems ( self ):
        """Summary items."""
        items = [ ( "Connectivity", True ) ]
        if self.atoms     is not None: items.append ( ( "Atoms"     , "{:d}".format ( len ( self.atoms     ) ) ) )
        if self.bonds     is not None: items.append ( ( "Bonds"     , "{:d}".format ( len ( self.bonds     ) ) ) )
        if self.angles    is not None: items.append ( ( "Angles"    , "{:d}".format ( len ( self.angles    ) ) ) )
        if self.dihedrals is not None: items.append ( ( "Dihedrals" , "{:d}".format ( len ( self.dihedrals ) ) ) )
        if self.isolates  is not None: items.append ( ( "Isolates"  , "{:d}".format ( len ( self.isolates  ) ) ) )
        if self.ringSets  is not None: items.append ( ( "Ring Sets" , "{:d}".format ( len ( self.ringSets  ) ) ) )
        return items

    def ToMapping ( self ):
        """Return a minimal bond mapping for serialization."""
        mapping = {}
        if len ( self.edges ) > 0:
            indices = self.nodeIndices
            bonds   = defaultdict ( list )
            for bond in self.edges:
                bonds[bond.type.label].append ( ( indices[bond.node1], indices[bond.node2] ) )
            mapping["bonds"] = bonds
        return mapping

    # . Properties - ICs.
    @property
    def angles ( self ):
        items = self.__dict__.get ( "_angles", None )
        if items is None:
            items = self._MakeAngles ( )
            self.__dict__["_angles"] = items
        return items

    @property
    def atoms ( self ): return self.nodes

    @property
    def bonds ( self ): return self.edges

    @property
    def dihedrals ( self ):
        items = self.__dict__.get ( "_dihedrals", None )
        if items is None:
            items = self._MakeDihedrals ( )
            self.__dict__["_dihedrals"] = items
        return items

    @property
    def isolates ( self ):
        items = self.__dict__.get ( "_isolates", None )
        if items is None:
            items = ConnectedComponents ( self )
            self.__dict__["_isolates"] = items
        return items

    @property
    def ringSets  ( self ):
        items = self.__dict__.get ( "_ringSets", None )
        if items is None:
            items = BiconnectedComponents ( self )
            self.__dict__["_ringSets"] = items
        return items

    # . Properties - indices.
    @property
    def angleIndices ( self ):
        items = self.__dict__.get ( "_angleIndices", None )
        if items is None:
            items = self._IndicesAngle ( )
            self.__dict__["_angleIndices"] = items
        return items

    @property
    def bondIndices ( self ):
        items = self.__dict__.get ( "_bondIndices", None )
        if items is None:
            items = self._IndicesBond ( )
            self.__dict__["_bondIndices"] = items
        return items

    @property
    def dihedralIndices ( self ):
        items = self.__dict__.get ( "_dihedralIndices", None )
        if items is None:
            items = self._IndicesDihedral ( )
            self.__dict__["_dihedralIndices"] = items
        return items

    @property
    def isolateIndices ( self ):
        items = self.__dict__.get ( "_isolateIndices", None )
        if items is None:
            items = self._IndicesIsolate ( )
            self.__dict__["_isolateIndices"] = items
        return items

    @property
    def nodeIndices ( self ):
        items = self.__dict__.get ( "_nodeIndices", None )
        if items is None:
            items = self._IndicesNode ( )
            self.__dict__["_nodeIndices"] = items
        return items

    @property
    def ringSetIndices ( self ):
        items = self.__dict__.get ( "_ringSetIndices", None )
        if items is None:
            items = self._IndicesIsolate ( )
            self.__dict__["_ringSetIndices"] = items
        return items

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
