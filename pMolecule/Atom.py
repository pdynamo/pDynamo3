"""Classes for handling atoms."""

from collections         import defaultdict
from pCore               import logFile            , \
                                LogFileActive      , \
                                Selection          , \
                                ShallowClone       , \
                                TreeLeafNode       , \
                                SummarizableObject
from  pScientific        import PeriodicTable
from  pScientific.Arrays import Array
from  pScientific.Graph  import Node

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class Atom ( Node, TreeLeafNode ):
    """An atom."""

    #
    # . Other actual or possible attributes:
    #
    #   aromaticElectrons       the number of electrons the atom donates to an aromatic system
    #   aromaticValence         the extra valence used up by belonging to an aromatic system
    #   chiralityClass          Smiles
    #   chiralityNumber         Smiles
    #   implicitHydrogens       hydrogens with no atom specification
    #   inRing                  is the atom in a ring?
    #   isotope                 isotope specification
    #   isReduced               Smiles
    #   radical                 radical specification
    #
    # . No attributes are inherited from Node.
    #
    # . Attributes/properties inherited from TreeLeafNode:
    #
    #   label
    #   parent
    #   path
    #
    _attributable = dict ( TreeLeafNode._attributable )
    _attributable.update ( { "atomicNumber"   :    -1 ,
                             "connections"    :     0 ,
                             "formalCharge"   :     0 ,
                             "geometry"       :  None ,
                             "hydrogens"      :     0 ,
                             "index"          :    -1 ,
                             "isAromatic"     : False ,
                             "oxidationState" :  None ,
                             "valence"        :     0 } )

    def __getattr__ ( self, name ):
        """Get an elemental attribute."""
        try   : return getattr ( PeriodicTable[self.__dict__.get ( "atomicNumber", -1 )], name )
        except: raise AttributeError ( "Unknown atom or elemental attribute: " + name + "." )

#    def __getstate__ ( self ): return self.__dict__ # . Needed for unspecified attributes?

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class AtomContainer ( SummarizableObject ):
    """A container class for atoms."""

    _attributable = dict ( SummarizableObject._attributable )
    _classLabel   = "Atom Container"
    _attributable.update ( { "items" : list } )

    def __getitem__ ( self, i ):
        """Get an item."""
        if isinstance ( i, int ):
            try   : return self.items[i]
            except: raise IndexError ( "Container index out of range." )
        else: raise TypeError ( "Container index must be integer." )

    def __len__ ( self ):
        """Length."""
        return len ( self.items )

    def _SetItemAttributes ( self, mapping ):
        """Set item attributes given a mapping."""
        for ( attribute, data ) in mapping.items ( ):
            if isinstance ( data, dict ):
                for ( index, value ) in data.items ( ):
                    setattr ( self.items[index], attribute, value )
            else:
                for ( index, value ) in enumerate ( data ):
                    setattr ( self.items[index], attribute, value )

    def _SetItemsFromIterable ( self, iterable ):
        """The iterable sequence must consist of Atom objects, atomic numbers or symbols.
        """
        # . Test for an iterable sequence.
        try:    items = iter ( iterable )
        except: TypeError ( "Argument initializer must be an iterable sequence." )
        # . Create the atoms by either using an atom directly or searching for an atomic number.
        atoms     = []
        unlabeled = 0
        for item in items:
            if isinstance ( item, Atom ): atom = item
            else:
                if   isinstance ( item, float          ): atomicNumber = PeriodicTable.AtomicNumberFromMass ( item )
                elif isinstance ( item, int            ): atomicNumber = item
                elif isinstance ( item, str            ): atomicNumber = PeriodicTable.AtomicNumber ( item )
                elif hasattr    ( item, "atomicNumber" ): atomicNumber = item.atomicNumber
                else: raise TypeError ( "Unrecognized atom initializer: " + str ( item ) + "." )
                atom = Atom.WithOptions ( atomicNumber = atomicNumber )
            atoms.append ( atom )
            if atom.label is None: unlabeled += 1
        self.items = atoms
        # . Make labels and reindex.
        if unlabeled > 0: self.MakeLabels ( )
        self.Reindex ( )

    def ElementDecomposition ( self, selection = None ):
        """Return a dictionary containing the indexes of the atoms of each element."""
        if selection is None: indices = range ( len ( self ) )
        else:                 indices = selection
        decomposition = defaultdict ( list )
        for i in indices: decomposition[self.items[i].atomicNumber].append ( i )
        return decomposition

    def ElementFrequencies ( self, selection = None ):
        """Return a dictionary containing the element frequencies."""
        if selection is None: indices = range ( len ( self ) )
        else:                 indices = selection
        frequencies = defaultdict ( int )
        for i in indices: frequencies[self.items[i].atomicNumber] += 1
        return frequencies

    def FormulaString ( self ):
        """Return a formula string."""
        frequencies = self.ElementFrequencies ( )
        tokens      = [ "{:s}{:d}".format ( PeriodicTable.Symbol ( n ), frequencies[n] ) for n in sorted ( frequencies.keys ( ) ) ]
        return ( "".join ( tokens ) )

    @classmethod
    def FromIterable ( selfClass, iterable, attributes = None ):
        """Constructor from iterable."""
        self = selfClass ( )
        self._SetItemsFromIterable ( iterable )
        if attributes is not None: self._SetItemAttributes ( attributes )
        return self

    @classmethod
    def FromMapping ( selfClass, mapping ):
        """Constructor from mapping."""
        atomicNumbers = mapping.pop ( "atomicNumber" )
        return selfClass.FromIterable ( atomicNumbers, attributes = mapping )

    def MakeLabels ( self, startFromOne = False, useFrequencies = False ):
        """Make default atom labels."""
        if startFromOne: i0 = 1
        else:            i0 = 0
        if useFrequencies:
            frequencies = defaultdict ( int )
            for ( i, atom ) in enumerate ( self.items ):
                atom.label = PeriodicTable.Symbol ( atom.atomicNumber, index = frequencies[atom.atomicNumber]+i0 )
                frequencies[atom.atomicNumber] += 1
        else:
            for ( i, atom ) in enumerate ( self.items ):
                atom.label = PeriodicTable.Symbol ( atom.atomicNumber, index = i+i0 )

    @classmethod
    def Merge ( selfClass, items, information = {} ):
        """Merging."""
        atomIncrements = []
        atomMapping    = []
        atoms          = []
        numberAtoms    = 0
        for item in items:
            if isinstance ( item, AtomContainer ):
                atomIncrements.append ( numberAtoms )
                mapping      = {}
                numberAtoms += len ( item )
                for atom in item:
                    newAtom       = ShallowClone ( atom )
                    mapping[atom] = newAtom
                    atoms.append ( newAtom )
                atomMapping.append ( mapping )
            else:
                raise ValueError ( "Invalid atom container in merge." )
        merged                          = selfClass ( )
        merged.__dict__["items"]        = atoms
        information["Atom Container"  ] = merged
        information["Atom Mapping"    ] = atomMapping
        information["Index Increments"] = atomIncrements
        return merged

    def Prune ( self, selection, information = {} ):
        """Pruning."""
        atomMapping = {}
        pruned      = self.__class__ ( )
        items       = []
        for s in selection:
            oldAtom              = self.items[s]
            newAtom              = ShallowClone ( oldAtom )
            atomMapping[oldAtom] = newAtom
            items.append ( newAtom )
        pruned.items = items
        information["Atom Container"] = pruned
        information["Atom Mapping"  ] = atomMapping
        return pruned

    def Reindex ( self ):
        """Reindex the container items."""
        for ( i, atom ) in enumerate ( self.items ):
            atom.index = i

    def SummaryItems ( self ):
        """Summary items."""
        frequencies = self.ElementFrequencies ( )
        atoms       = len ( self )
        hydrogens   = frequencies.pop ( 1, 0 )
        unknowns    = 0
        for ( n, f ) in frequencies.items ( ):
            if n <= 0: unknowns += f
        heavies = atoms - hydrogens - unknowns
        items = [ ( "Atoms"       , True                        ) ,
                  ( "Atoms"       , "{:d}".format ( atoms     ) ) ,
                  ( "Heavy Atoms" , "{:d}".format ( heavies   ) ) ,
                  ( "Hydrogens"   , "{:d}".format ( hydrogens ) ) ]
        if unknowns > 0: items.append ( ( "Unknowns" , "{:d}".format ( unknowns  ) ) )
        return items

    def ToMapping ( self ):
        """Return a mapping for serialization."""
        mapping      = { "atomicNumber" : [ item.atomicNumber for item in self.items ] }
        formalCharge = { index : item.formalCharge for ( index, item ) in enumerate ( self.items ) if item.formalCharge != 0 }
        if len ( formalCharge ) > 0: mapping["formalCharge"] = formalCharge
        return { "atoms" : mapping }

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
