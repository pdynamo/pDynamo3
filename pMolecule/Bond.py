"""Classes for handling bonds in a connectivity."""

import operator

from  enum              import Enum
from  pScientific.Graph import Edge
from .Atom              import Atom

# . Aromatic bond types are not defined. Instead bonds themselves are flagged if they form part of an aromatic system.

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
class BondType ( Enum ):
    """Bond definitions."""
    Undefined      = ( 0 , "Undefined" )
    Null           = ( 0 , "Null"      )
    Single         = ( 1 , "Single"    )
    Double         = ( 2 , "Double"    )
    Triple         = ( 3 , "Triple"    )
    Quadruple      = ( 4 , "Quadruple" )
#    Amide          = ( 1 , False , "Amide"          ) #????

    def __init__ ( self, bondOrder, label ):
        """Constructor."""
        self.bondOrder = bondOrder
        self.label     = label

# . Bond labels to definition type mapping.
BondTypeFromLabel = { member.label : member for member in BondType.__members__.values ( ) }

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class Bond ( Edge ):
    """A class to represent a bond."""

    _attributable = dict ( Edge._attributable )
    _attributable.update ( { "isAromatic" : False              ,
                             "type"       : BondType.Undefined } )

    def _CheckOptions ( self ):
        """Check options."""
        super ( Bond, self )._CheckOptions ( )
        if not ( isinstance ( self.node1 , Atom     ) and \
                 isinstance ( self.node2 , Atom     ) and \
                 isinstance ( self.type  , BondType ) ): raise TypeError ( "Invalid bond node or type." )

    @classmethod
    def FromIterable ( selfClass, value, nodes ):
        """A constructor for use with Connectivity's BondsFromIterable."""
        self = None
        # . A bond.
        if isinstance ( value, selfClass ) and ( value.node1 in nodes ) and ( value.node2 in nodes ):
            self = value
        # . A tuple of items ( i/node1, j/node2 [, type [, isAromatic ] ] ).
        else:
            bondType   = None
            isAromatic = False
            node1      = value[0]
            node2      = value[1]
            if isinstance ( node1, int ): node1 = nodes[node1]
            if isinstance ( node2, int ): node2 = nodes[node2]
            if len ( value ) > 2:
                bondType = value[2]
                if isinstance ( bondType, str ): bondType = BondTypeFromLabel[bondType]
                if len ( value ) > 3:
                    isAromatic = value[3]
            self = selfClass.WithNodes ( node1, node2, isAromatic = isAromatic, type = bondType )
        return self

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
