"""Classes and functions for reading SMILES strings.

Chirality is parsed but not handled.
"""

import string

from  pCore             import AttributableObject         , \
                               logFile
from  pMolecule         import AddExplicitHydrogens       , \
                               Atom                       , \
                               Bond                       , \
                               BondType                   , \
                               Connectivity               , \
                               ConvertInputConnectivity   , \
                               System
from .SMILESDefinitions import _BondTypeFromSymbol        , \
                               _ChargeTokens              , \
                               _ChiralityClassConnections , \
                               _ChiralityDefaultClasses   , \
                               _ChiralityMaximumNumbers   , \
                               _ChiralityFull             , \
                               _ChiralityReduced          , \
                               _ElementsAromatic          , \
                               _ElementsOrganic           , \
                               _ElementTokens             , \
                               _HCountTokens              , \
                               _ReducedElementTokens      , \
                               _ValenciesOrganic
                               
#===================================================================================================================================
# . Class for a parsable string.
#===================================================================================================================================
class _ParsableString:

    def __init__ ( self, inputString ):
        """Constructor."""
        # . Save the input string and set the current position.
        self.inputString = inputString
        self.position    = 0
        # . Convert the input string to a list of characters tagged with their positions in the input string.
        self.characters = [ ]
        for ( i, c ) in enumerate ( inputString ): self.characters.append ( ( c, i ) )
        # . Initialize the bracket data.
        # . Bracketpairs stores the length of each bracket pair.
        self.bracketCache = [ ]
        self.bracketPairs = { }

    def CacheBrackets ( self ):
        """Add a bracket pair to the bracket cache."""
        # . Decrement self.position.
        self.position = self.position - 1
        # . Get the position of the closing bracket corresponding to the opening bracket at the current position.
        closing = self.position + self.bracketPairs[self.characters[self.position][1]]
        # . Save the unparsed bracket string.
        self.bracketCache.extend ( self.characters[self.position:closing+1] )
        # . Delete details about the bracket pair.
        del ( self.bracketPairs[self.characters[self.position][1]] )
        # . Delete the unparsed characters.
        del ( self.characters[self.position:closing+1] )

    def Finalize ( self ):
        """Finish parsing."""
        # . Check that all the string has been parsed.
        if self.position < len ( self.characters ): raise SMILESReaderError ( "There are unparsed characters.", parseString = self )
        # . Check that the bracket cache is empty.
        if len ( self.bracketCache ) > 0: raise SMILESReaderError ( "There are unparsed brackets.", parseString = self )

    def GetCharacter ( self, characters, default = None ):
        """Get a single character from the current position."""
        if self.position < len ( self.characters ):
            c = self.characters[self.position][0]
            if c in characters:
                self.position = self.position + 1
                return c
        return default

    def GetCharacters ( self, characters, default = None ):
        """Get a string of characters from the current position."""
        characterList = [ ]
        while self.position < len ( self.characters ):
            c = self.characters[self.position][0]
            if c in characters:
                characterList.append ( c )
                self.position = self.position + 1
            else:
                break
        if len ( characterList ) > 0: return "".join ( characterList )
        else:                         return default

    def GetInteger ( self, default = None ):
        """Get an integer from the string."""
        integer = self.GetCharacters ( string.digits )
        if integer is None: return default
        else:               return int ( integer )

    def GetToken ( self, tokens, MaximumTokenLength, default = None, firstCharacters = None ):
        """Get a token from the string."""
        # . Check that there are more characters to parse.
        if self.position >= len ( self.characters ): return default
        # . Check the first characters (done for speed).
        if firstCharacters is not None:
            if self.characters[self.position][0] not in firstCharacters: return default
        # . Get MaximumTokenLength characters from the string (or as near as possible to it).
        characterList = [ ]
        position      = self.position
        for length in range ( MaximumTokenLength ):
            if position < len ( self.characters ): characterList.append ( self.characters[position][0] )
            else:                                  break
            position = position + 1
        tokenString = "".join ( characterList )
        tokenLength = len ( tokenString )
        # . Loop over token lengths in reverse order.
        for length in range ( tokenLength, 0, -1 ):
            if tokenString[0:length] in tokens:
                self.position = self.position + length
                return tokens[tokenString[0:length]]
        # . Nothing found.
        return default

    def GetTokenInteger ( self, tokens, MaximumTokenLength, defaultInteger = None, defaultToken = None, firstCharacters = None ):
        """Get a token from the string with a following integer."""
        # . Save the position.
        position = self.position
        # . Get the token.
        token = self.GetToken ( tokens, MaximumTokenLength, default = defaultToken, firstCharacters = firstCharacters )
        # . Find the following integer.
        if position != self.position: integer = self.GetInteger ( default = defaultInteger )
        else:                         integer = defaultInteger
        # . Return the results.
        return ( token, integer )

    def HasCharacters ( self ):
        """Check whether there are characters remaining in the string."""
        return self.position < len ( self.characters )

    def ParseBrackets ( self ):
        """Parse brackets to check syntax and set up appropriate data structures."""
        # . Initialization.
        copen = 0
        popen = 0
        sopen = 0
        plength = 0
        psize   = [ ]
        pstart  = [ ]
        # . Loop over the characters.
        self.position = 0
        while self.position < len ( self.characters ):
            c = self.characters[self.position][0]
            # . Check for bracket characters.
            # . Parentheses.
            if c == "(":
                if ( copen == 1 ) or ( sopen == 1 ): raise SMILESReaderError ( "Cannot nest ( within [] or {} pairs.", parseString = self )
                psize.append  ( plength       )
                pstart.append ( self.position )
                plength = 0
                popen   = popen + 1
            elif c == ")":
                if ( popen   <  1 ): raise SMILESReaderError ( "Unpaired closing ).", parseString = self )
                if ( plength <= 0 ): raise SMILESReaderError ( "Empty () pair.",      parseString = self )
                opening = pstart.pop ( )
                self.bracketPairs[opening] = self.position - opening
                plength = psize.pop ( )
                popen   = popen - 1
            # . Square brackets.
            elif c == "[":
                if ( copen == 1 ) or ( sopen == 1 ): raise SMILESReaderError ( "Cannot nest [ within [] or {} pairs.", parseString = self )
                sopen  = 1
                sstart = self.position
            elif c == "]":
                if ( sopen                      != 1 ): raise SMILESReaderError ( "Unpaired closing ].", parseString = self )
                if ( self.position - sstart - 1 <= 0 ): raise SMILESReaderError ( "Empty [] pair.",      parseString = self )
                sopen = 0
            # . Curly brackets.
            elif c == "}":
                if ( copen == 1 ) or ( sopen == 1 ): raise SMILESReaderError ( "Cannot nest { within [] or {} pairs.", parseString = self )
                copen  = 1
                cstart = self.position
            elif c == "{":
                if ( copen                      != 1 ): raise SMILESReaderError ( "Unpaired closing }.", parseString = self )
                if ( self.position - cstart - 1 <= 0 ): raise SMILESReaderError ( "Empty {} pair.",      parseString = self )
                copen = 0
            # . Increment plength.
            if ( c != "(" ) and ( c != ")" ): plength = plength + 1
            # . Increment the current position.
            self.position = self.position + 1
        # . Do some final checks.
        if ( copen != 0 ):
            self.position = cstart
            raise SMILESReaderError ( "Unpaired opening {.", parseString = self )
        if ( popen  > 0 ):
            self.position = pstart[-1]
            raise SMILESReaderError ( "Unpaired opening (.", parseString = self )
        if ( sopen != 0 ):
            self.position = sstart
            raise SMILESReaderError ( "Unpaired opening [.", parseString = self )
        # . Reinitialize the current position.
        self.position = 0

    def ParsedString ( self ):
        """Return the current parsed string as a string."""
        characterList = [ ]
        for character in self.characters: characterList.append ( character[0] )
        return "".join ( characterList )

    def RemoveWhiteSpace ( self ):
        """Remove white space."""
        for character in self.characters:
            if character[0] in string.whitespace: self.characters.remove ( character )

    def StatusString ( self, positions = None ):
        """Return the input string and its current position."""
        length = len ( self.inputString )
        if length > 0:
            # . Get the indices at which to place markers.
            if positions is None: indices = [ self.position ]
            else:                 indices = positions
            n        = len ( self.characters )
            cIndices = set ( )
            for index in indices:
                if   index >= n: cIndices.add ( n )
                elif index >= 0: cIndices.add ( self.characters[index][1] )
            # . Generate the marker string.
            cIndices = list ( cIndices )
            cIndices.sort ( )
            p = cIndices[-1]
            marker = ( p + 1 ) * [ " " ]
            for index in cIndices: marker[index] = "^"
            # . Create the string.
            s = " " + self.inputString + "\n" + " " + "".join ( marker )
        else:
            s = ""
        return s

    def UncacheBrackets ( self ):
        """Uncache any brackets."""
        # . Check for unparsed brackets.
        if len ( self.bracketCache ) > 0:
            # . Add the bracket cache at the current position.
            self.characters[self.position:self.position] = self.bracketCache
            # . Reinitialize the bracket cache.
            self.bracketCache = [ ]

#===================================================================================================================================
# . Error classes.
#===================================================================================================================================
class SMILESReaderError ( Exception ):

    def __init__ ( self, *arguments, **options ):
        """Constructor."""
        super ( SMILESReaderError, self ).__init__ ( *arguments )
        if ( len ( arguments ) == 0 ) or ( not isinstance ( arguments[0], str ) ): self.args = ( "SMILES reader error.", )
        parseString = options.get ( "parseString", None )
        positions   = options.get ( "positions"  , None )
        if parseString is not None:
            self.args = ( self.args[0] + "\n" + parseString.StatusString ( positions = positions ), )

#===================================================================================================================================
# . SMILES reader class.
#===================================================================================================================================
class SMILESReader ( AttributableObject ):
    """SMILESReader is the class for reading SMILES strings."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "atomPositions" : dict  ,
                             "connectivity"  : None  ,
                             "isParsed"      : False ,
                             "parsedSmiles"  : None  ,
                             "smilesString"  : None  } )

    def _CheckOptions ( self ):
        """Check options."""
        super ( SMILESReader, self )._CheckOptions ( )
        if not isinstance ( self.smilesString, str ): raise TypeError ( "Invalid input SMILES string." )

    def Parse ( self ):
        """Parse a SMILES string."""
        if not self.isParsed:
            try:
                if len ( self.smilesString ) > 0:

                    # . Start a connectivity.
                    self.connectivity = Connectivity ( )

                    # . Initialization.
                    bracketLastAtom = [ ]
                    crossLink       = { }
                    lastAtom        = None
                    newAtom         = None

                    # . Initialize the input string for parsing.
                    parseString = _ParsableString ( self.smilesString )
                    parseString.RemoveWhiteSpace ( )
                    parseString.ParseBrackets    ( )

                    # . Loop over all the characters in the string.
                    while parseString.HasCharacters ( ):

                        # . Check for parentheses.
                        c = parseString.GetCharacter ( [ "(", ")" ], default = " " )

                        # . An opening parenthesis.
                        if c == "(":
                            # . There is a root atom so parse this bracket.
                            if lastAtom is not None:
                                # . Save and reset lastAtom.
                                bracketLastAtom.append ( lastAtom )
                                lastAtom = None
                            # . There is no root atom.
                            else:
                                # . Cache this bracket pair for later.
                                parseString.CacheBrackets ( )

                        # . A closing parenthesis.
                        elif c == ")":
                            # . Reset lastAtom.
                            lastAtom = bracketLastAtom.pop ( )

                        # . An atom or crossLink specification.
                        else:

                            # . For the first atom in the system set a NullBond character.
                            if len ( self.connectivity ) == 0:
                                ( bondType, bondIsAromatic ) = ( BondType.Null, False )
                            # . Check for an explicit bond character.
                            else:
                                ( bondType, bondIsAromatic ) = parseString.GetToken ( _BondTypeFromSymbol, 1, default = ( BondType.Single, False ) )

                            # . Check for a crossLink or a full atom specification.
                            c = parseString.GetCharacter ( "%[" + string.digits, default = " " )

                            # . Parse the atom or crossLink string.
                            if ( c == "%" ) or ( c in string.digits ):

                                # . Check that lastAtom exists.
                                if lastAtom is None: raise SMILESReaderError ( "A cross-link specification must have a preceding atom.", parseString = parseString )

                                # . Parse the crossLink string.
                                if c == "%":
                                    link = parseString.GetInteger ( )
                                    if link is None: raise SMILESReaderError ( "The cross-link integer is missing.", parseString = parseString )
                                else:
                                    link = int ( c )

                                # . Check to see if this cross-link already exists.
                                if link in crossLink:

                                    # . Get the link data.
                                    ( connectedAtom, firstType, firstIsAromatic ) = crossLink[link]

                                    # . If this key is already associated with lastAtom there is an error.
                                    if connectedAtom is lastAtom:
                                        raise SMILESReaderError ( "An atom has two identical cross-link keys.", parseString = parseString )

                                    # . The key is not associated with lastAtom.
                                    else:

                                        # . Get the type of bond to add to the bond list.
                                        if bondType is BondType.Single:
                                            bondType = firstType
                                        # . The bond specifications are incompatible.
                                        elif ( bondType is firstType ) and ( firstType is not BondType.Single ):
                                            raise SMILESReaderError ( "A pair of cross-link keys has incompatible bond-type specifications.", parseString = parseString )

                                        # . Add this bond to the bond list (if it is not NullBond) and delete the associated dictionary entry.
                                        if bondType != BondType.Null:
                                            if lastAtom.isAromatic and connectedAtom.isAromatic and ( not bondIsAromatic ) and ( bondType is BondType.Single ):
                                                bondIsAromatic = True
                                                bondType       = BondType.Undefined
                                            self.connectivity.AddEdge ( Bond.WithNodes ( lastAtom, connectedAtom, isAromatic = bondIsAromatic, type = bondType ) )
                                        del crossLink[link]

                                # . The key does not exist so add it to the dictionary.
                                else:
                                    crossLink[link] = ( lastAtom, bondType, bondIsAromatic )

                            # . An atom specification.
                            else:

                                # . Full atom specification.
                                if c == "[":

                                    # . Set isReduced.
                                    isReduced = False

                                    # . Parse the individual atom fields.
                                    # . Isotope.
                                    isotope = parseString.GetInteger ( default = 0 )

                                    # . Atom symbol.
                                    smilesPosition               = parseString.position
                                    ( atomicNumber, isAromatic ) = parseString.GetToken ( _ElementTokens, 2, default = ( -1, False ) )
                                    if atomicNumber == -1: raise SMILESReaderError ( "Missing or unrecognized element symbol.", parseString = parseString )

                                    # . Chirality (full and reduced notations).
                                    ( chiralityClass, chiralityNumber ) = parseString.GetToken ( _ChiralityReduced, 4, default = ( None, 0 ), firstCharacters = [ "@" ] )
                                    if chiralityClass is not None:
                                        if chiralityNumber == 1:
                                            ( chiralityClass, chiralityNumber ) = parseString.GetTokenInteger ( _ChiralityFull, 2, defaultInteger = 0, defaultToken = None )
                                            if chiralityClass is None:
                                                chiralityClass  = "??"
                                                chiralityNumber = 1

                                    # . Hydrogen count.
                                    ( hCount, multiplier ) = parseString.GetTokenInteger ( _HCountTokens, 1, defaultInteger = 1, defaultToken = 0 )
                                    hCount *= multiplier

                                    # . Charge.
                                    charge = parseString.GetToken ( _ChargeTokens, 4, default = 0, firstCharacters = [ "+", "-" ] )
                                    if ( charge == 1 ) or ( charge == -1 ):
                                        multiplier = parseString.GetInteger ( default = 1 )
                                        charge = charge * multiplier

                                    # . Check the current character which must be a "]".
                                    if parseString.GetCharacter ( [ "]" ], default = " " ) != "]": raise SMILESReaderError ( "A closing ] is missing.", parseString = parseString )

                                # . Reduced atom specification.
                                else:

                                    # . Set isReduced and other defaults.
                                    charge    = 0
                                    hCount    = 0
                                    isotope   = 0
                                    isReduced = True

                                    # . Parse the individual atom fields.
                                    # . Atom symbol.
                                    smilesPosition               = parseString.position
                                    ( atomicNumber, isAromatic ) = parseString.GetToken ( _ReducedElementTokens, 2, default = ( -1, False ) )
                                    if atomicNumber == -1: raise SMILESReaderError ( "Missing or unrecognized organic-subset element symbol.", parseString = parseString )

                                    # . Chirality.
                                    ( chiralityClass, chiralityNumber ) = parseString.GetToken ( _ChiralityReduced, 4, default = ( None, 0 ), firstCharacters = [ "@" ] )

                                # . Create a new atom and set its properties.
                                newAtom = Atom.WithOptions ( atomicNumber = atomicNumber )
                                newAtom.chiralityClass        = chiralityClass
                                newAtom.chiralityNumber       = chiralityNumber
                                newAtom.formalCharge          = charge
                                newAtom.implicitHydrogens     = hCount
                                newAtom.isAromatic            = isAromatic
                                newAtom.isotope               = isotope
                                newAtom.isReduced             = isReduced
                                self.atomPositions[newAtom]   = smilesPosition

                                # . Put all the atom data in the atom list.
                                self.connectivity.AddNode ( newAtom )

                                # . Get lastAtom for the first atom of a branch.
                                if ( lastAtom is None ) and ( len ( bracketLastAtom ) > 0 ): lastAtom = bracketLastAtom[-1]

                                # . Save the current bond.
                                if ( bondType != BondType.Null ) and ( lastAtom is not None ):
                                    if lastAtom.isAromatic and newAtom.isAromatic and ( not bondIsAromatic ) and ( bondType is BondType.Single ):
                                        bondIsAromatic = True
                                        bondType       = BondType.Undefined
                                    self.connectivity.AddEdge ( Bond.WithNodes ( newAtom, lastAtom, isAromatic = bondIsAromatic, type = bondType ) )

                                # . Reset lastAtom.
                                lastAtom = newAtom

                                # . Uncache any unparsed bracket pairs.
                                parseString.UncacheBrackets ( )

                    # . Do some last checks.
                    # . Check parseString.
                    parseString.Finalize ( )

                    # . Check to see that the cross-link dictionary is empty.
                    if len ( crossLink ) > 0: raise SMILESReaderError ( "There are unparsed crossLink keys: " + " ".join ( [ "{:d}".format ( i ) for i in crossLink.keys ( ) ] ) )

                    # . Save the final parsed SMILES.
                    self.parsedSmiles = parseString.ParsedString ( )

                    # . Finish up.
                    ConvertInputConnectivity ( self.connectivity, _ValenciesOrganic )
                    self.VerifyChirality     ( parseString )

                # . Zero length string.
                else:
                    self.parsedSmiles = ""

            # . Errors.
            except Exception as error:
                raise SMILESReaderError ( error.args[0] )

            # . Flag the SMILES as having been parsed.
            self.isParsed = True

    @classmethod
    def StringToSystem ( selfClass, smiles, **options ):
        """Generate a system from a SMILES string."""
        reader = selfClass.WithOptions ( smilesString = smiles )
        reader.Parse ( )
        return reader.ToSystem ( **options )

    def ToSystem ( self, **options ):
        """Return a system."""
        if self.connectivity is not None:
            if options.get ( "addExplicitHydrogens", True ): AddExplicitHydrogens ( self.connectivity )
            system       = System.FromConnectivity ( self.connectivity )
            system.label = self.parsedSmiles
        else:
            system       = None
        return system

    def VerifyChirality ( self, parseString ):
        """Verify the chirality of the atoms in a connectivity but do nothing with it."""
        # . This method needs atom connections so do once these have been computed.
        for atom in self.connectivity.atoms:
            if atom.chiralityClass is not None:
                # . Reduced chirality class.
                if atom.chiralityClass == "??":
                    atom.chiralityClass = _ChiralityDefaultClasses.get ( atom.connections, None )
                    if atom.chiralityClass is None:
                        raise SMILESReaderError ( "An atom has a number of connections, {:d}, for which there is no default chirality class.".format ( atom.connections ), positions = [ self.atomPositions[atom] ] )
                # . Check the number of connections and the number of the class.
                if atom.connections != _ChiralityClassConnections[atom.chiralityClass]:
                    raise SMILESReaderError ( "An atom has a number of connections which is incompatible with the chirality class.", positions = [ self.atomPositions[atom] ] )
                if ( atom.chiralityNumber <= 0 ) or ( atom.chiralityNumber > _ChiralityMaximumNumbers[atom.chiralityClass] ):
                    raise SMILESReaderError ( "The number of a specified chirality class is invalid.", positions = [ self.atomPositions[atom] ] )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
