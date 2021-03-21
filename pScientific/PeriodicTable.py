"""Classes for handling elements."""

import glob, math, os, os.path

from collections.abc import Mapping
from pCore           import logFile                  , \
                            LogFileActive            , \
                            SummarizableObject       , \
                            YAMLPickle               , \
                            YAMLPickleFileExtension  , \
                            YAMLUnpickle             

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Unknown element data.
_UnknownElementNumber = -1
_UnknownElementSymbol = "*"

# . YAML tag.
_YAMLTag = "!Element"

# . It might be good to have defaults for some missing element attributes without having to check explicitly ...

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def _KeywordToAttributeLabel ( keyword ):
    """Convert a keyword to an attribute label."""
    tokens = keyword.split ( )
    first  = tokens.pop ( 0 ).lower ( )
    if len ( tokens ) == 0: return first
    else:                   return "".join ( [ first ] + tokens )

#===================================================================================================================================
# . Element class.
#===================================================================================================================================
class Element:
    """A chemical element."""

    def __lt__ ( self, other ):
        """Less than comparison."""
        return ( self.atomicNumber < other.atomicNumber )

    @classmethod
    def FromMapping ( selfClass, mapping, keysToAttributes ):
        """Constructor from a mapping."""
        self = selfClass ( )
        for ( key, value ) in mapping.items ( ):
            self.__dict__[keysToAttributes.get ( key, key )] = value
        return self

    def GetCoordinationAngle ( self, connections ):
        """Get a coordination angle given the number of connections."""
        angle  = None
        angles = self.__dict__.get ( "coordinationAngles", None )
        if angles is not None: angle = angles.get ( connections, None )
        return angle

    def GetSingleBondDistance ( self, atomicNumber ):
        """Get a single bond distance to another element."""
        r     = None
        bonds = self.__dict__.get ( "singleBondDistances", None )
        if bonds is not None:
            r = bonds.get ( atomicNumber, None )
        return r

    def ToMapping ( self, attributesToKeys ):
        """Return a mapping."""
        return { attributesToKeys.get ( name, name ) : value for ( name, value ) in self.__dict__.items ( ) }

#===================================================================================================================================
# . Element container class.
#===================================================================================================================================
# . Elements can be accessed by atomic number or symbol.
class ElementContainer ( Mapping, SummarizableObject ):
    """A container for chemical elements."""

    # . No other _summarizable.
    _attributable = dict ( SummarizableObject._attributable )
    _classLabel   = "Element Container"
    _attributable.update ( { "attributesToKeys" : None ,
                             "items"            : list ,
                             "keysToAttributes" : None ,
                             "numberIndex"      : dict ,
                             "symbolIndex"      : dict } )

    def __contains__ ( self, key ):
        """Contains using an integer atomic number or a string symbol."""
        return ( ( key in self.numberIndex ) or ( key in self.symbolIndex ) )

    def __getitem__ ( self, key ):
        """Get an item using an integer atomic number or a string symbol."""
        if   key in self.numberIndex: return self.numberIndex[key]
        elif key in self.symbolIndex: return self.symbolIndex[key]
        else: raise KeyError ( "Unknown atomic number or symbol \"{:s}\".".format ( repr ( key ) ) )

    def __iter__ ( self ):
        return iter ( self.items )

    def __len__ ( self ):
        return len ( self.items )

    def AtomicNumber ( self, inString ):
        """Find the atomic number of an element from a string."""
        number = _UnknownElementNumber
        if len ( inString ) > 0:
            # . The string contains a symbol.
            if inString[0:1].isalpha ( ):
                a = inString[0:1].upper ( )
                if ( len ( inString ) > 1 ) and inString[1:2].isalpha ( ): symbol = a + inString[1:2].lower ( )
                else:                                                      symbol = a
                if symbol in self.symbolIndex: number = self.symbolIndex[symbol].atomicNumber
            # . The string is an integer.
            else:
                ndigits = 0
                for c in inString:
                    if c.isdigit ( ): ndigits += 1
                    else:             break
                if ndigits > 0: number = int ( inString[0:ndigits] )
        return number

    def AtomicNumberFromMass ( self, mass ):
        """Find the atomic number given a mass."""
        # . This is not a fool-proof method. It returns the element with the
        # . minimum absolute difference between the input and reference masses.
        minimumDifference = math.fabs ( mass - self.items[0].mass )
        number            = self.items[0].atomicNumber
        for item in self.items[1:]:
            difference = math.fabs ( mass - item.mass )
            if difference < minimumDifference:
                minimumDifference = difference
                number            = item.atomicNumber
        return number

    def AtomicNumberFromSymbol ( self, inString ):
        """Find the atomic number of an element given a symbol."""
        number = _UnknownElementNumber
        if len ( inString ) > 0:
            if inString[0:1].isalpha ( ):
                a = inString[0:1].upper ( )
                if ( len ( inString ) > 1 ) and inString[1:2].isalpha ( ): symbol = a + inString[1:2].lower ( )
                else:                                                      symbol = a
                if  symbol in self.symbolIndex: number = self.symbolIndex[symbol].atomicNumber
        return number

    @staticmethod
    def CreateAttributesFromMappingKeys ( keys ):
        """Create attributes from mapping keys."""
        attributesToKeys = {}
        keysToAttributes = {}
        for key in keys:
            label = _KeywordToAttributeLabel ( key )
            attributesToKeys[label] = key
            keysToAttributes[key  ] = label
        return ( attributesToKeys, keysToAttributes )

    @classmethod
    def FromElements ( selfClass, items ):
        """Constructor from elements."""
        self       = selfClass ( )
        self.items = sorted ( items )
        for item in self.items:
            atomicNumber = item.atomicNumber
            symbol       = item.symbol
            if ( atomicNumber not in self.numberIndex ) and ( symbol not in self.symbolIndex ):
                self.numberIndex[atomicNumber] = item
                self.symbolIndex[symbol      ] = item
        return self

    @classmethod
    def FromParameterDirectory ( selfClass, path = None ):
        """Constructor from a directory of parameters."""
        if path is None: path = os.path.join ( os.getenv ( "PDYNAMO3_PARAMETERS" ), "elements" )
        filePaths = glob.glob ( os.path.join ( path, "*" + YAMLPickleFileExtension ) )
        keys      = set ( )
        mappings  = []
        for filePath in filePaths:
            mapping = YAMLUnpickle ( filePath )
            keys.update ( mapping.keys ( ) )
            mappings.append ( mapping )
        ( attributesToKeys, keysToAttributes ) = selfClass.CreateAttributesFromMappingKeys ( keys )
        items                 = [ Element.FromMapping ( mapping, keysToAttributes ) for mapping in mappings ]
        self                  = selfClass.FromElements ( items )
        self.attributesToKeys = attributesToKeys
        self.keysToAttributes = keysToAttributes
        return self

    def SummaryItems ( self ):
        """Summarizing."""
        return [ ( "Number of Elements", "{:d}".format ( len ( self.items ) ) ) ]

    def Symbol ( self, atomicNumber, index = None ):
        """Return a symbol given an atomic number with an optional index."""
        if atomicNumber in self.numberIndex: symbol = self.numberIndex[atomicNumber].symbol
        else:                                symbol = _UnknownElementSymbol
        if index is None: return symbol
        else:             return symbol + repr ( index )

    def Symbols ( self ):
        """Return the symbols of the elements in the container."""
        symbols = list ( self.symbolIndex.keys ( ) )
        symbols.sort ( )
        return symbols

    def ToParameterDirectory ( self, path ):
        """Create a parameter directory from the container."""
        if not os.path.exists ( path ): os.mkdir ( path )
        for ( key, value ) in self.symbolIndex.items ( ):
            mapping = value.ToMapping ( self.attributesToKeys )
            YAMLPickle ( os.path.join ( path, key + YAMLPickleFileExtension ), mapping, defaultFlowStyle = True, title = _YAMLTag ) 

    @property
    def atomicNumbers ( self ):
        """Return the atomic numbers of the elements in the container."""
        return sorted ( self.numberIndex.keys ( ) )

#===================================================================================================================================
# . The periodic table and related functions.
#===================================================================================================================================
# . Definition.
PeriodicTable       = ElementContainer.FromParameterDirectory ( )
PeriodicTable.label = "Periodic Table"

# . Element inquiry functions taking either an atomic number or symbol.
# . Qq none of these.
def IsMainGroup ( key ):
    """Is the element main-group?"""
    element = PeriodicTable.get ( key, None )
    return ( ( element is not None ) and ( element.group in ( 1, 2, 13, 14, 15, 16, 17, 18 ) ) )

def IsMetal ( key ):
    """Is the element a metal?"""
    element = PeriodicTable.get ( key, None )
    return ( ( element is not None ) and ( element.atomicNumber > 0 ) and ( not IsNonMetal ( key ) ) and ( not IsSemiMetal ( key ) ) )

def IsNonMetal ( key ):
    """Is the element a non-metal?"""
    element = PeriodicTable.get ( key, None )
    return ( ( element is not None ) and ( element.atomicNumber in ( 1, 2, 6, 7, 8, 9, 10, 15, 16, 17, 18, 34, 35, 36, 53, 54, 86 ) ) )

def IsSemiMetal ( key ):
    """Is the element a semi-metal?"""
    # . 118 uncertain for the moment.
    element = PeriodicTable.get ( key, None )
    return ( ( element is not None ) and ( element.atomicNumber in ( 5, 14, 32, 33, 51, 52, 85, 118 ) ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
