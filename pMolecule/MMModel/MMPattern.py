"""MM patterns."""

from  .MMModelError        import MMModelError
from ..ConnectivityPattern import AtomPattern                  , \
                                  BondPattern                  , \
                                  ConnectivityPattern          , \
                                  ConnectivityPatternContainer

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMPattern ( ConnectivityPattern ):
    """An MM pattern used for assigning MM atom types."""

    _attributable = dict ( ConnectivityPattern._attributable )
    _attributable.update ( { "atomTypes" : None } )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MMPatternContainer ( ConnectivityPatternContainer ):
    """A container for MM patterns."""

    _attributable = dict ( ConnectivityPatternContainer._attributable )
    _patternClass = MMPattern
    _attributable.update ( { "mmAtomTypes" : None         ,
                             "termLabel"   : "MM Pattern" } )

    #yaml_tag = "!MMPatternContainer"

    def IndexAtomTypes ( self, mmAtomTypes ):
        """Index the atom types in the container's patterns."""
        if ( mmAtomTypes is not None ) and ( self.mmAtomTypes is None ):
            for item in self.items:
                # . Types.
                item.atomTypes = []
                for label in item.nodeResults["atomTypeLabel"]:
                    atomType = mmAtomTypes.GetItem ( label )
                    if ( atomType is None ): raise MMModelError ( "Unknown atom type in MM pattern: " + label + "." )
                    else:                    item.atomTypes.append ( atomType )
                # . Charges.
                charges = item.nodeResults.get ( "charge", None )
                if ( charges is None ) or ( len ( charges ) == 0 ):
                    item.nodeResults["charge"] = [ atomType.charge for atomType in item.atomTypes ]
            # . Finish up.
            self.mmAtomTypes = mmAtomTypes

    def TypeConnectivity ( self, connectivity, atomTypes, atomCharges, untypedAtoms ):
        """Type the atoms in a connectivity."""
        if self.mmAtomTypes is not None:
            for pattern in self.items:
                pattern.MakeConnections ( )
                matches = pattern.FindAllMatches ( connectivity, selection = untypedAtoms )
                # . Maybe here should check that there are not intersecting matches within the same call (selection takes care of the others).
                # . I.e. an atom has been matched more than once and with different atom types and/or charges - in this case which definition to take?
                # . Easiest if raise error. No problem if type/charge same.
                # . THIS NEEDS TO BE DONE.
                if ( matches is not None ) and ( len ( matches ) > 0 ):
                    charges      = pattern.nodeResults["charge"]
                    matchedAtoms = []
                    for match in matches:
                        for ( atomType, charge, isSupporting, iNode ) in zip ( pattern.atomTypes, charges, pattern.isSupporting, match ):
                            i = connectivity.nodeIndices[iNode]
                            if ( atomType is not None ) and ( atomTypes[i] is None ) and ( not isSupporting ):
                                atomCharges[i] = charge
                                atomTypes  [i] = atomType.label
                                matchedAtoms.append ( iNode )
                                hydrogenType   = atomType.hydrogenType
                                if hydrogenType is not None:
                                    hCharge = hydrogenType.charge
                                    for hNode in connectivity.adjacentNodes[iNode]:
                                        if hNode.atomicNumber == 1:
                                            h = connectivity.nodeIndices[hNode]
                                            atomCharges[h]  = hCharge
                                            atomCharges[i] -= hCharge
                                            atomTypes  [h]  = hydrogenType.label
                                            matchedAtoms.append ( hNode )
                    for node in matchedAtoms: untypedAtoms.remove ( node )
                    if len ( untypedAtoms ) <= 0: break

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
