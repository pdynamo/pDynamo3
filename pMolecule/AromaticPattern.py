"""Aromatic patterns."""

import os, os.path

from  pCore                         import YAMLMappingFile_ToObject     , \
                                           YAMLPickleFileExtension
from  pMolecule.ConnectivityPattern import ConnectivityPattern          , \
                                           ConnectivityPatternContainer

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class AromaticPattern ( ConnectivityPattern ):
    """An aromatic pattern used for assigning aromaticity."""
    pass

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class AromaticPatternContainer ( ConnectivityPatternContainer ):
    """A container for aromatic patterns."""

    _attributable = dict ( ConnectivityPatternContainer._attributable )
    _patternClass = AromaticPattern
    _attributable.update ( { "termLabel" : "Aromatic Pattern" } )

    #yaml_tag = "!AromaticPatternContainer"

    @classmethod
    def FromParameterDirectory ( selfClass, path = None ):
        """Constructor from parameter directory."""
        if path is None: path = os.path.join ( os.getenv ( "PDYNAMO3_PARAMETERS" ), "chemistry", "aromaticPatterns{:s}".format ( YAMLPickleFileExtension ) )
        return YAMLMappingFile_ToObject ( path, selfClass )

    def TypeConnectivity ( self, connectivity, untypedAtoms ):
        """Determine which atoms could be aromatic."""
        electrons = {}
        for pattern in self.items:
            pattern.MakeConnections ( )
            matches = pattern.FindAllMatches ( connectivity, selection = untypedAtoms )
            if ( matches is not None ) and ( len ( matches ) > 0 ):
                donatedElectrons = pattern.nodeResults["donatedElectrons"]
                matchedAtoms     = []
                for match in matches:
                    for ( donated, isSupporting, node ) in zip ( donatedElectrons, pattern.isSupporting, match ):
                        if not isSupporting:
                            electrons[node] = donated
                            matchedAtoms.append ( node )
                for node in matchedAtoms: untypedAtoms.discard ( node )
                if len ( untypedAtoms ) <= 0: break
        return electrons

#===================================================================================================================================
# . Standard aromatic patterns with and without implicit hydrogens.
#===================================================================================================================================
# . Paths.
path = os.path.join ( os.getenv ( "PDYNAMO3_PARAMETERS" ), "chemistry", "aromaticPatternsImplicitHydrogens{:s}".format ( YAMLPickleFileExtension ) )

# . Definitions.
AromaticPatterns                  = AromaticPatternContainer.FromParameterDirectory ( )
AromaticPatternsImplicitHydrogens = AromaticPatternContainer.FromParameterDirectory ( path = path )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :

    print ( "\nAromatic Patterns:" )
    for ( i, item ) in enumerate ( AromaticPatterns.items ):
        print ( "{:5d}   {:s}".format ( i+1, item.label ) )
    print ( "\nAromatic Patterns for Implicit Hydrogens:" )
    for ( i, item ) in enumerate ( AromaticPatternsImplicitHydrogens.items ):
        print ( "{:5d}   {:s}".format ( i+1, item.label ) )
