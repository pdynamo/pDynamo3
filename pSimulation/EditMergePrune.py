"""Functions for editing, merging and pruning systems."""

#
# . Editing and edit/merging:
#   - are system-specific.
#   - requires unique atom labels within a system and unique system keys.
#     The atom labels in the final system are constructed as sKey:aLabel and so are also unique.
#   - return a new system that consists of a connectivity and coordinates only (no sequence, etc.).
#   - coordinates are not guessed for any atoms that are added.
#
# . Merging and pruning:
#   - work for any object that have Merge and Prune methods.
#   - errors are, in general, ignored. Any problems that arise are treated by simply 
#     skipping the offending attribute and omitting them from the final object.
#

from pMolecule             import Atom   , \
                                  Bond   , \
                                  System
from pScientific.Geometry3 import Coordinates3

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class EditMergePruneError ( Exception ):
    pass

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Atom keywords to keep in edited systems.
_AtomOptions = { "atomicNumber" : -1 ,
                 "formalCharge" :  0 }

# . The new label separator.
_LabelSeparator = ":"

# . Undefined coordinate value.
_XYZZero = 999999.0

#===================================================================================================================================
# . Editing utility functions.
#===================================================================================================================================
def _MakeAtomLabel ( sKey, aLabel ):
    """Make a new atom label."""
    if ( sKey is None ) or ( sKey == "" ): return aLabel
    else:                                  return ( sKey + _LabelSeparator + aLabel )

def _MakeBondKey ( iLabel, jLabel ):
    """Make a bond dictionary key."""
    return ( max ( iLabel, jLabel ), min ( iLabel, jLabel ) )

def _ProcessSystem ( sKey, system, atomsToDelete, bondsToDelete, atomsToAdd, index0 ):
    """Get the old and new atoms and old bonds of an input system."""
    # . Atom labels must be unique.
    labels = set ( [ atom.label for atom in system.atoms ] )
    if len ( labels ) != len ( system.atoms ): raise EditMergePruneError ( "Atom labels within a system must be unique." )
    # . Atoms to delete must exist.
    toRemove = set ( atomsToDelete )
    if not toRemove.issubset ( labels ): raise EditMergePruneError ( "Unknown atoms to delete." )
    # . Atoms to add must not have duplicate labels.
    toAdd = set ( [ atom.label for atom in atomsToAdd ] )
    if len ( ( labels - toRemove ) & toAdd ) > 0: raise EditMergePruneError ( "Duplicate atoms to add." )
    # . Initialization.
    atoms          = []
    bondDictionary = {}
    oldToNew       = {}
    oldToNewXYZ    = []
    index          = index0
    # . Get atoms with new labels.
    for ( i, atom ) in enumerate ( system.atoms ):
        if atom.label not in toRemove:
            options          = { key : getattr ( atom, key, default ) for ( key, default ) in _AtomOptions.items ( ) }
            options["label"] = _MakeAtomLabel ( sKey, atom.label )
            newAtom          = Atom.WithOptions ( **options )
            atoms.append ( newAtom )
            oldToNew[atom.label] = newAtom
            oldToNewXYZ.append ( ( i, index ) )
            index += 1
    if len ( toAdd ) > 0:
        for atom in atomsToAdd:
            options          = { key : getattr ( atom, key, default ) for ( key, default ) in _AtomOptions.items ( ) }
            options["label"] = _MakeAtomLabel ( sKey, atom.label )
            newAtom          = Atom.WithOptions ( **options )
            atoms.append ( newAtom )
            oldToNew[atom.label] = newAtom
            index += 1
    # . Get bonds between atoms that have not been deleted.
    bondDictionary = {}
    for bond in system.connectivity.bonds:
        aLabel1 = bond.node1.label
        aLabel2 = bond.node2.label
        if ( aLabel1 not in toRemove ) and ( aLabel2 not in toRemove ):
            iLabel = _MakeAtomLabel ( sKey, aLabel1 )
            jLabel = _MakeAtomLabel ( sKey, aLabel2 )
            bKey   = _MakeBondKey ( iLabel, jLabel )
            bondDictionary[bKey] = ( oldToNew[bond.node1.label], oldToNew[bond.node2.label], bond.type )
    # . Bonds to remove.
    if bondsToDelete is not None:
        for ( aLabel1, aLabel2 ) in bondsToDelete:
            iLabel = _MakeAtomLabel ( sKey, aLabel1 )
            jLabel = _MakeAtomLabel ( sKey, aLabel2 )
            bKey   = _MakeBondKey ( iLabel, jLabel )
            old    = bondDictionary.pop ( bKey, None )
            if old is None: raise EditMergePruneError ( "Unknown bond to delete." )
    # . Finish up.
    return ( atoms, bondDictionary, oldToNew, oldToNewXYZ, index )

#===================================================================================================================================
# . Editing.
#===================================================================================================================================
# . Arguments:
#
#  - system             the system to be edited
#  - atomsToAdd         list of atom instances to add
#  - atomsToDelete      list of atom labels to delete
#  - bondsToDelete      list of ( aLabel1, aLabel2 ) pairs indicating the bonds to delete
#  - bondsToAdd         list of ( aLabel1, aLabel2, bondType ) tuples specifying the bonds to add to the edited system
#
def EditByAtom ( system               ,
                 atomsToAdd    = None ,
                 atomsToDelete = None ,
                 bondsToAdd    = None ,
                 bondsToDelete = None ):
    """Edit an input system."""
    # . Check input arguments.
    if atomsToAdd    is None: atomsToAdd    = []
    if atomsToDelete is None: atomsToDelete = []
    if bondsToDelete is None: bondsToDelete = []
    # . Process the system for old and new atoms and old bonds.
    ( atoms, bondDictionary, oldToNew, oldToNewXYZ, index ) = _ProcessSystem ( None          ,
                                                                               system        ,
                                                                               atomsToDelete ,
                                                                               bondsToDelete ,
                                                                               atomsToAdd    ,
                                                                               0             )
    # . Create bonds.
    bonds = [ Bond.WithNodes ( atom1, atom2, type = bondType ) for ( atom1, atom2, bondType ) in bondDictionary.values ( ) ]
    if bondsToAdd is not None:
        for ( iLabel, jLabel, bondType ) in bondsToAdd:
            if ( iLabel not in oldToNew ) or ( jLabel not in oldToNew ):
                raise EditMergePruneError ( "Adding bond for unknown atom." )
            bKey = _MakeBondKey ( iLabel, jLabel )
            if bKey in bondDictionary: raise EditMergePruneError ( "Adding duplicate bond." )
            bonds.append ( Bond.WithNodes ( oldToNew[iLabel], oldToNew[jLabel], type = bondType ) )
    # . Create coordinates.
    newXYZ = Coordinates3.WithExtent ( len ( atoms ) )
    newXYZ.Set ( _XYZZero )
    oldXYZ = system.coordinates3
    for ( old, new ) in oldToNewXYZ:
        for c in range ( 3 ):
            newXYZ[new,c] = oldXYZ[old,c]
    # . Create the edited and merged system.
    edited = System.FromAtoms ( atoms, bonds = bonds )
    edited.label        = "Edited system"
    edited.coordinates3 = newXYZ
    return edited

#===================================================================================================================================
# . Editing with merging.
#===================================================================================================================================
# . Arguments:
#
#  - systems            list of ( sKey, system ) pairs specifying the systems to be merged (not a dictionary, as order may be important)
#  - atomsToAdd         dictionary[sKey] = list of atom instances to add
#  - atomsToDelete      dictionary[sKey] = list of atom labels to delete
#  - bondsToDelete      dictionary[sKey] = list of ( aLabel1, aLabel2 ) pairs indicating the bonds to delete
#  - bondsToAdd         list of ( sKey1, aLabel1, sKey2, aLabel2, bondType ) tuples specifying the bonds to add to the merged system
#
def EditMergeByAtom ( systems              ,
                      atomsToAdd    = None ,
                      atomsToDelete = None ,
                      bondsToAdd    = None ,
                      bondsToDelete = None ):
    """Edit and merge systems."""
    # . Check input arguments.
    if atomsToAdd    is None: atomsToAdd    = {}
    if atomsToDelete is None: atomsToDelete = {}
    if bondsToDelete is None: bondsToDelete = {}
    sKeys = set ( [ sKey for ( sKey, system ) in systems ] )
    if len ( sKeys ) != len ( systems ): raise EditMergePruneError ( "System labels must be unique." )
    for argument in ( atomsToAdd, atomsToDelete, bondsToDelete ):
        keys = set ( argument.keys ( ) )
        if not keys.issubset ( sKeys ): raise EditMergePruneError ( "Unknown system key in input argument." )
    # . Process systems for old and new atoms and old bonds.
    atoms          = [] # . List of atoms in final system.
    bondDictionary = {} # . Dictionary of all bonds in final system.
    oldToNew       = {} # . oldToNew   [sKey] = mapping of old label to new atom.
    oldToNewXYZ    = {} # . oldToNewXYZ[sKey] = list of ( old index in system, new index in final system ) pairs.
    index          = 0
    for ( sKey, system ) in systems:
        ( atoms0, bondDictionary0, oldToNew0, oldToNewXYZ0, index0 ) = _ProcessSystem ( sKey                           ,
                                                                                        system                         ,
                                                                                        atomsToDelete.get ( sKey, [] ) ,
                                                                                        bondsToDelete.get ( sKey, [] ) ,
                                                                                        atomsToAdd.get    ( sKey, [] ) ,
                                                                                        index                          )
        atoms.extend          ( atoms0          )
        bondDictionary.update ( bondDictionary0 )
        oldToNew   [sKey] = oldToNew0
        oldToNewXYZ[sKey] = oldToNewXYZ0
        index             = index0
    # . Create bonds.
    bonds = [ Bond.WithNodes ( atom1, atom2, type = bondType ) for ( atom1, atom2, bondType ) in bondDictionary.values ( ) ]
    if bondsToAdd is not None:
        for ( ( sKey1, aLabel1 ), ( sKey2, aLabel2 ), bondType ) in bondsToAdd:
            if ( sKey1 not in oldToNew ) or ( aLabel1 not in oldToNew[sKey1] ) or \
               ( sKey2 not in oldToNew ) or ( aLabel2 not in oldToNew[sKey2] ):
                raise EditMergePruneError ( "Adding bond for unknown atom." )
            iLabel = _MakeAtomLabel ( sKey1, aLabel1 )
            jLabel = _MakeAtomLabel ( sKey2, aLabel2 )
            bKey   = _MakeBondKey   ( iLabel, jLabel )
            if bKey in bondDictionary: raise EditMergePruneError ( "Adding duplicate bond." )
            bonds.append ( Bond.WithNodes ( oldToNew[sKey1][aLabel1], oldToNew[sKey2][aLabel2], type = bondType ) )
    # . Create coordinates.
    newXYZ = Coordinates3.WithExtent ( len ( atoms ) )
    newXYZ.Set ( _XYZZero )
    for ( sKey, system ) in systems:
        oldXYZ = system.coordinates3
        for ( old, new ) in oldToNewXYZ[sKey]:
            for c in range ( 3 ):
                newXYZ[new,c] = oldXYZ[old,c]
    # . Create the edited and merged system.
    edited = System.FromAtoms ( atoms, bonds = bonds )
    edited.label        = "Edited and merged system"
    edited.coordinates3 = newXYZ
    return edited

#===================================================================================================================================
# . Merging.
#===================================================================================================================================
def MergeByAtom ( *arguments, **options ):
    """Merge the arguments."""
    merged = None
    if len ( arguments ) > 0:
        if isinstance ( arguments[0], ( list, tuple ) ): items = arguments[0]
        else:                                            items = arguments
        if hasattr ( items[0].__class__, "Merge" ): return items[0].__class__.Merge ( items, information = options )
        else: raise EditMergePruneError ( "Cannot merge objects of type {!r}.".format ( type ( items[0] ) ) )
    return merged

def MergeRepeatByAtom ( item, repeat, **options ):
    """Return a merged object consisting of |repeat| copies of |item|."""
    merged = None
    if repeat <= 0:
        raise EditMergePruneError ( "Invalid repeat argument." )
    elif repeat == 1:
        merged = item
    elif hasattr ( item, "MergeRepeat" ):
        merged = item.MergeRepeat ( repeat, information = options )
    elif hasattr ( item.__class__, "Merge" ):
        merged = item.__class__.Merge ( repeat * [ item ], information = options )
    else:
        raise EditMergePruneError ( "Cannot merge objects of type {!r}.".format ( type ( item ) ) )
    return merged

#===================================================================================================================================
# . Pruning.
#===================================================================================================================================
def PruneByAtom ( item, selection, **options ):
    """Prune |item| by |selection|."""
    if hasattr ( item, "Prune" ): return item.Prune ( selection, information = options )
    else: raise EditMergePruneError ( "Cannot prune object of type {!r}.".format ( type ( item ) ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
