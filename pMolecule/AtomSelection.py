"""Atom selection classes and functions."""

from pCore                 import Selection
from pScientific.Geometry3 import CrossPairList_FromSingleCoordinates3

#
# . Use as (for example):
#
#   from pMolecule import AtomSelection
#   newSelection = AtomSelection.ByComponent ( system, oldSelection )
#

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class AtomSelectionError ( Exception ): pass

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def ByBondedNeighbor ( system, selection, iterations = 1 ):
    """Expand the selection by bonded neighbors."""
    if ( len ( system.connectivity.bonds ) > 0 ) and ( iterations > 0 ):
        # . Get all connected atoms.
        atoms       = system.connectivity.nodes
        neighbors   = system.connectivity.adjacentNodes
        nodeToIndex = system.connectivity.nodeIndices
        indices     = set ( selection )
        for iteration in range ( iterations ):
            new = set ( )
            for index in indices:
                new.update ( [ nodeToIndex[n] for n in neighbors[atoms[index]] ] )
            indices.update ( new )
        # . Finish up.
        return Selection.FromIterable ( indices )
    else: return selection

def ByComponent ( system, selection ):
    """Expand the selection by component."""
    # . Get all the unique components.
    unique = set ( )
    for index in selection:
        unique.add ( system.atoms[index].parent )
    # . Get the indices of all atoms.
    indices = set ( )
    for component in unique:
        for atom in component.children: indices.add ( atom.index )
    # . Finish up.
    return Selection.FromIterable ( indices )

def ByEntity ( system, selection ):
    """Expand the selection by entity."""
    # . Get all the unique entities.
    unique = set ( )
    for index in selection:
        unique.add ( system.atoms[index].parent.parent )
    # . Get the indices of all atoms.
    indices = set ( )
    for entity in unique:
        for component in entity.children:
            for atom in component.children: indices.add ( atom.index )
    # . Finish up.
    return Selection.FromIterable ( indices )

def ByIsolate ( system, selection ):
    """Expand the selection by isolate."""
    if len ( system.connectivity.bonds ) > 0:
        isolates  = system.connectivity.isolateIndices
        inIsolate = { i : n for ( n, isolate ) in enumerate ( isolates ) for i in isolate }
        indices   = set ( )
        for s in selection:
            indices.update ( isolates[inIsolate[s]] )
        return Selection.FromIterable ( indices )
    else: return selection

def ByLinearPolymer ( system, selection ):
    """Expand the selection by linear polymer."""
    # . Get all the unique linear polymers.
    linearPolymerIndex = system.sequence.linearPolymerIndex
    unique             = set ( )
    for index in selection:
        polymer = linearPolymerIndex.get ( system.atoms[index].parent, None )
        if polymer is not None: unique.add ( polymer )
    # . Get the indices of all atoms.
    indices = set ( )
    for index in unique:
        polymer = system.sequence.linearPolymers[index]
        for component in polymer.ComponentIterator ( ):
            for atom in component.children: indices.add ( atom.index )
    # . Finish up.
    return Selection.FromIterable ( indices )

def ByRingSet ( system, selection ):
    """Expand the selection by ring set."""
    if len ( system.connectivity.bonds ) > 0:
        ringSets  = system.connectivity.ringSetIndices
        inRingSet = { i : n for ( n, ringSet ) in enumerate ( ringSets ) for i in ringSet }
        indices   = set ( )
        for s in selection:
            indices.update ( ringSets[inRingSet[s]] )
        return Selection.FromIterable ( indices )
    else: return selection

def FromAtomPattern ( system, atomPattern ):
    """Get a selection given an atom pattern."""
    # . Initialization.
    indices  = set ( )
    sequence = system.sequence
    # . Parse the pattern to get the label selection fields.
    ( eSFields, cSFields, aSFields ) = sequence.ParseAtomPattern ( atomPattern )
    # . Check for matches.
    for entity in sequence.children:
        if sequence.FieldsLabelMatch ( eSFields, entity.label ):
            for component in entity.children:
                if sequence.FieldsLabelMatch ( cSFields, component.label ):
                    for atom in component.children:
                        if sequence.FieldsLabelMatch ( aSFields, atom.label ): indices.add ( atom.index )
    # . Finish up.
    return Selection.FromIterable ( indices )

def Within ( system, selection, distance, excludeSelf = False ):
    """Selection of atoms within a given distance."""
    # . Find all atoms within distance of selection.
    pairlist = CrossPairList_FromSingleCoordinates3 ( system.coordinates3, selection1 = selection, safety = distance )
    indices  = set ( )
    for ( i, j ) in pairlist: indices.add ( j )
    # . Add self?
    if not excludeSelf:
        for j in selection: indices.add ( j )
    # . Finish up.
    return Selection.FromIterable ( indices )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
