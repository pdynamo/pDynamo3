#===================================================================================================================================
# . Classes and functions to read Mopac input files.
#===================================================================================================================================

import math

from  pCore                 import logFile        , \
                                   LogFileActive  , \
                                   Selection      , \
                                   TextFileReader 
from  pMolecule             import System                                           
from  pScientific           import PeriodicTable                                    
from  pScientific.Geometry3 import Coordinates3   , \
                                   Vector3
from .ExportImport          import _Importer

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . Mopac element symbols.
_DUMMYATOM           = "XX"
_NONDUMMYATOMSYMBOLS = [ "TV", "CB", "+3", "++", "+", "Fr", "At", "-", "--", "-3" ]
_MOPACATOMSYMBOLS    = [ _DUMMYATOM ] + _NONDUMMYATOMSYMBOLS

# . Tolerance for determing whether points are co-linear.
#_LINEARTOLERANCE = 1.0e-1

# . Multiplicity keywords.
_Multiplicities = { "SINGLET" : 1 ,
                    "DOUBLET" : 2 ,
                    "TRIPLET" : 3 ,
                    "QUARTET" : 4 ,
                    "QUINTET" : 5 ,
                    "SEXTET"  : 6 ,
                    "SEPTET"  : 7 }

# . Reference data keywords.
_REFERENCEDATA    = { "CP" : "Constant Pressure Heat Capacity", "D" : "Dipole", "H" : "Enthalpy of Formation", "I" : "Ionization Potential", "IA" : "Ionization Potential", "IE" : "Ionization Potential", "IP" : "Ionization Potential", "S" : "Entropy" }
_REFERENCESOURCES = { "DR" : "Dipole Source", "GR" : "Free Energy Source", "HR" : "Enthalpy Source", "IR" : "Ionization Potential Source" }

# . Symmetry keywords.
_SYMMETRIES = set ( [ "SYM", "SYMM", "SYMMETRY" ] )

#===================================================================================================================================
# . Mopac input file reader class.
#===================================================================================================================================
class MopacInputFileReader ( TextFileReader ):
    """Read Mopac input files."""

    _classLabel = "Mopac Input File Reader"

    def FilterAtoms ( self, filterlevel ):
        """Filter out Mopac-specific atoms."""
        if self.isParsed and ( filterlevel > -3 ):
            atomicNumbers = []
            selected      = []
            for ( i, a ) in enumerate ( self.atomicNumbers ):
                if a >= filterlevel:
                    atomicNumbers.append ( a )
                    selected.append      ( i )
            if len ( atomicNumbers ) < len ( self.atomicNumbers ):
                coordinates3 = Coordinates3.WithExtent ( len ( atomicNumbers ) )
                coordinates3.Gather ( self.coordinates3, Selection.FromIterable ( selected ) )
                self.atomicNumbers = atomicNumbers
                self.coordinates3  = coordinates3

    @staticmethod
    def GetAtomicNumber ( token ):
        """Convert a token into an atomic number."""
        symbol = token.upper ( )
        if   symbol == _DUMMYATOM:           return -3
        elif symbol in _NONDUMMYATOMSYMBOLS: return -2
        else: return PeriodicTable.AtomicNumber ( token )

    def Parse ( self, log = logFile, filterlevel = -1, QGEOMETRY = True ):
        """Parse the data on the file."""
        if not self.isParsed:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
            # . Open the file.
            self.Open ( )
            try:
                # . Keyword line.
                self._ParseKeywords ( )
                # . Title line.
                self._ParseTitle    ( )
                # . Extra information.
                self._ParseExtras   ( )
                # . Geometry.
                if QGEOMETRY:
                    if self.QXYZ: self._ParseXYZ     ( )
                    else:         self._ParseZMatrix ( )
            except EOFError:
                pass
            # . Close the file.
            self.WarningStop ( )
            self.Close ( )
            # . Set the parsed flag and some other options.
            self.log     = None
            self.isParsed = True
            # . Filter out Mopac-specific atoms.
            self.FilterAtoms ( filterlevel )

    @classmethod
    def PathToCoordinates3 ( selfClass, path ):
        """Return the coordinates from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse ( )
        return inFile.ToCoordinates3 ( )

    @classmethod
    def PathToSystem ( selfClass, path ):
        """Return the system from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse ( )
        return inFile.ToSystem ( )

    def ToCoordinates3 ( self ):
        """Return a coordinates3 object."""
        if self.isParsed: return self.coordinates3
        else:            return None

    def ToElectronicState ( self, frameIndex = -1 ):
        """Return an electronic state object."""
        if self.isParsed:
            frame = self.frames[frameIndex]
            return { "Charge"       :                 getattr ( self, "charge"       , 0         )  ,
                     "Multiplicity" : _Multiplicities[getattr ( self, "multiplicity" , "SINGLET" )] }
        else:
            return None

    def ToSystem ( self ):
        """Return a system."""
        system = None
        if self.isParsed:
            system              = System.FromAtoms ( self.atomicNumbers )
            system.label        = self.title
            system.coordinates3 = self.ToCoordinates3 ( )
            system.Set ( "Electronic State", self.ToElectronicState ( ), attribute = "electronicState" )
        return system

    # . Private methods.
    def _GetKeywords ( self ):
        """Get a line and parse into keypairs and keywords."""
        tokens   = self.GetTokens ( )
        keypairs = {}
        keywords = set ( )
        for token in tokens:
            item = token.upper ( )
            pair = item.split ( "=", 1 )
            if len ( pair ) > 1:
                key   = pair[0]
                value = pair[1]
                if ( key in keypairs ) and ( keypairs[key] != value ):
                    self.Warning ( "Duplicate keyword (" + key + ") with different values (" + keypairs[key] + " and " + value + ").", False )
                else:
                    keypairs[key] = value
            else:
                keywords.add ( item )
        return ( keypairs, keywords )

    def _ParseExtras ( self ):
        """Parse the extras line."""
        # . Get keypairs and keywords.
        ( keypairs, keywords ) = self._GetKeywords ( )
        # . Reference data.
        self.referencedata = {}
        for ( key, value ) in keypairs.items ( ):
            if key in _REFERENCEDATA:
                try:
                    # . Take care of simple floating point values and combinations of the form "(x,y)".
                    value = value.replace ( "(", " " )
                    value = value.replace ( ",", " " )
                    value = value.replace ( ")", " " )
                    self.referencedata[_REFERENCEDATA[key]] = float ( value.split ( )[0] )
                except:
                    self.Warning ( "Error converting reference data value (" + key + " = " + value + ").", False )
            elif key  in _REFERENCESOURCES:
                self.referencedata[_REFERENCESOURCES[key]] = value

    def _ParseKeywords ( self ):
        """Parse the keyword line."""
        # . Initialization of items that are recognized.
        self.charge       = 0
        self.multiplicity = "SINGLET"
        self.QSYMMETRY    = False
        self.QXYZ         = False
        # . Get keypairs and keywords.
        ( keypairs, keywords ) = self._GetKeywords ( )
        # . Process the keywords.
        # . Charge.
        charge = keypairs.get ( "CHARGE", "0" )
        try:    self.charge = int ( charge )
        except: self.Warning ( "Invalid CHARGE value (" + charge + ").", False )
        # . Multiplicity.
        multiplicities = list ( keywords.intersection ( set ( _Multiplicities.keys ( ) ) ) )
        multiplicities.sort ( )
        if   len ( multiplicities ) == 0: self.multiplicity = "SINGLET"
        elif len ( multiplicities ) == 1: self.multiplicity = multiplicities[0]
        else: self.Warning ( "Duplicate MULTIPLICITY keywords (" + ", ".join ( multiplicities ) + ").", False )
        # . Symmetry.
        self.QSYMMETRY = ( len ( keywords.intersection ( _SYMMETRIES ) ) > 0 )
        # . XYZ input instead of
        self.QXYZ = ( "XYZ" in keywords )

    def _ParseTitle ( self ):
        """Parse the title line."""
        self.title = self.GetLine ( )

    def _ParseXYZ ( self ):
        """Parse a Cartesian coordinate geometry specification."""
        self.atomicNumbers = []
        xyz                = []
        while True:
            tokens  = self.GetTokens ( converters = [ MopacInputFileReader.GetAtomicNumber, float, None, float, None, float, None ] )
            ntokens = len ( tokens )
            if ntokens == 0:
                break
            elif ntokens >= 7:
                self.atomicNumbers.append ( tokens[0] )
                xyz.append ( [ tokens[1], tokens[3], tokens[5] ] )
            else:
                self.Warning ( "Insufficient number of tokens ({:d}) on Cartesian coordinate input line.".format ( ntokens ), True )
        # . Create the coordinates.
        self.coordinates3 = Coordinates3.WithExtent ( len ( self.atomicNumbers ) )
        self.coordinates3.Set ( 0.0 )
        for ( i, ( x, y, z ) ) in enumerate ( xyz ):
            self.coordinates3[i,0] = x
            self.coordinates3[i,1] = y
            self.coordinates3[i,2] = z

    def _ParseZMatrix ( self ):
        """Parse a Z-matrix geometry specification."""
        converters = [ MopacInputFileReader.GetAtomicNumber, float, None, float, None, float, None ]
        ncards     = 0
        oldfatal   = self.fatalErrors
        self.atomicNumbers = []
        zmatrixcards       = []
        while True:
            tokens  = self.GetTokens ( converters = converters )
            ntokens = len ( tokens )
            if ntokens == 0:
                break
            else:
                # . Append int to the converters for the first three cards.
                if ( ncards < 3 ): converters.append ( int )
                # . Basic checks on cards and tokens.
                if ( ( ncards == 0 ) and ( ntokens >=  7 ) ) or ( ( ncards == 1 ) and ( ntokens >=  8 ) ) or ( ( ncards == 2 ) and ( ntokens >=  9 ) ) or ( ( ncards >  2 ) and ( ntokens >= 10 ) ):
                    card = [ tokens[1], tokens[3], tokens[5] ]
                    if ( ncards > 0 ): card.append ( tokens[7] - 1 )
                    if ( ncards > 1 ): card.append ( tokens[8] - 1 )
                    if ( ncards > 2 ): card.append ( tokens[9] - 1 )
                    # . Check that the atom number values are valid.
                    QOK = True
                    for i in card[3:]:
                        QOK = QOK and ( i >= 0 ) and ( i < ncards )
                    if QOK:
                        ncards += 1
                        self.atomicNumbers.append ( tokens[0] )
                        zmatrixcards.append       ( card      )
                    else: self.Warning ( "Invalid atom number specification on Z-matrix input line.", True )
                else: self.Warning ( "Insufficient number of tokens ({:d}) on Z-matrix input line.".format ( ntokens ), True )
        # . Create the coordinates if there have been no fatal errors during the reading process.
        if self.fatalErrors == oldfatal:
            self.coordinates3 = Coordinates3.WithExtent ( len ( self.atomicNumbers ) )
            self.coordinates3.Set ( 0.0 )
            direction = Vector3.Null ( )
            for ( i, card ) in enumerate ( zmatrixcards ):
                if   i == 1:
                    direction.Set ( 0.0 )
                    direction[0] = 1.0
#                    direction.BasisVector ( 0 ) # . x-direction.
                    QOK = self.coordinates3.BuildPointFromDistance              ( i, card[3],                   card[0],          direction )
                elif i == 2:
                    direction.Set ( 0.0 )
                    direction[1] = 1.0
#                    direction.BasisVector ( 1 ) # . y-direction.
                    QOK = self.coordinates3.BuildPointFromDistanceAngle         ( i, card[3], card[4],          card[0], card[1], direction )
                elif i > 2:
                    QOK = self.coordinates3.BuildPointFromDistanceAngleDihedral ( i, card[3], card[4], card[5], card[0], card[1], card[2]   )
#                elif i > 1:
#                    # . Check for co-linear atoms.
#                    if math.fabs ( card[1] - 180.0 ) <= _LINEARTOLERANCE:
#                        j = card[3]
#                        k = card[4]
#                        for c in range ( 3 ):
#                            direction[c] = self.coordinates3[j,c] - self.coordinates3[k,c]
#                        direction.Normalize ( )
#                        QOK = self.coordinates3.BuildPointFromDistance ( i, j, card[0], direction )
#                    else:
#                        if i == 2:
#                            direction.BasisVector ( 1 ) # . y-direction.
#                            QOK = self.coordinates3.BuildPointFromDistanceAngle         ( i, card[3], card[4],          card[0], card[1], direction )
#                        else:
#                    # . Check for linear atoms.
#                            QOK = self.coordinates3.BuildPointFromDistanceAngleDihedral ( i, card[3], card[4], card[5], card[0], card[1], card[2]   )
                if not QOK:
                    self.Warning ( "Unable to build Cartesian coordinates from Z-matrix card number {:d}.".format ( i+1 ), True )
                    break

#===================================================================================================================================
# . Importer definitions.
#===================================================================================================================================
_Importer.AddHandler ( { Coordinates3 : MopacInputFileReader.PathToCoordinates3 ,
                         System       : MopacInputFileReader.PathToSystem       } ,
                       [ "mopin", "MOPIN" ], "MOPAC Input", defaultFunction = MopacInputFileReader.PathToSystem )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
