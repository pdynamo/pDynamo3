"""Connectivity pattern classes."""

import copy

from  pScientific.Graph import EdgePattern           , \
                               GraphPattern          , \
                               GraphPatternContainer , \
                               NodePattern
from .Bond              import BondType              , \
                               BondTypeFromLabel

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class AtomPattern ( NodePattern ):
    """A class to represent an atom pattern."""
    pass

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class BondPattern ( EdgePattern ):
    """A class to represent a bond pattern."""

    _attributable = dict ( EdgePattern._attributable )
    _attributable.update ( { "atomKey1"   : None ,
                             "atomKey2"   : None ,
                             "type"       : None ,
                             "typeObject" : None } )

    def Match ( self, bond ):
        """Match the pattern against a bond."""
        if self.typeObject is BondType.Undefined: isMatched = True
        else:  isMatched = ( self.typeObject is bond.type )
        return isMatched

    def SetOptions ( self, **options ):
        """Set options."""
        super ( BondPattern, self ).SetOptions ( **options )
        self.typeObject = BondTypeFromLabel.get ( self.type, None )
        if "atomKey1" in options: self.nodeKey1 = self.atomKey1
        if "atomKey2" in options: self.nodeKey2 = self.atomKey2

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class ConnectivityPattern ( GraphPattern ):
    """A class to represent a connectivity pattern."""

    _edgeClass = BondPattern
    _edgeTag   = "Bond"
    _nodeClass = AtomPattern
    _nodeTag   = "Atom"

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class ConnectivityPatternContainer ( GraphPatternContainer ):
    """A container for connectivity patterns."""

    _attributable = dict ( GraphPatternContainer._attributable )
    _patternClass = ConnectivityPattern
    _attributable.update ( { "termLabel" : "Connectivity Pattern" } )

    #yaml_tag = "!ConnectivityPatternContainer"

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
