"""Classes and functions for reading Jaguar input files."""

from  pCore                 import logFile        , \
                                   LogFileActive  , \
                                   TextFileReader
from  pMolecule             import System
from  pScientific           import PeriodicTable
from  pScientific.Arrays    import Array          , \
                                   StorageType
from  pScientific.Geometry3 import Coordinates3
from .ExportImport          import _Importer

#===================================================================================================================================
# . Class for storing orbital data.
#===================================================================================================================================
class JaguarOrbitalData:
    """A class for storing orbital data."""

    def __init__ ( self ):
        """Constructor."""
        pass

#===================================================================================================================================
# . Jaguar input file reader class.
#===================================================================================================================================
class JaguarInputFileReader ( TextFileReader ):
    """JaguarInputFileReader is the class for Jaguar input files that are to be read."""

    _classLabel = "Jaguar Input File Reader"

    @staticmethod
    def ToFloat ( token ):
        """Convert a token to a float."""
        if isinstance ( token, str ) and token.endswith ( "#" ): token = token[0:-1]
        return float ( token )

    def GetLine ( self, signalWarnings = True ):
        """Get a non-empty line with no warnings if an EOF is reached."""
        try:
            QFOUND = False
            while not QFOUND:
                line = next ( self.file ).strip ( )
                self.linesParsed += 1
                QFOUND = ( len ( line ) > 0 )
            return line
        except:
            if signalWarnings: self.Warning ( "Unexpected end-of-file.", True )
            raise EOFError

    def Parse ( self, log = logFile ):
        """Parse the data on the file."""
        if not self.isParsed:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
            # . Open the file.
            self.Open ( )
            try:
                # . Loop over the lines.
                while True:
                    line = self.GetLine ( signalWarnings = False )
                    # . Entry name.
                    if line.startswith ( "entry name:" ):
                        entryName = line[11:].strip ( )
                        if len ( entryName ) > 0: self.entryName = entryName
                    # . Other sections.
                    elif line.startswith ( "&gen"   ) or line.startswith ( "$gen"   ): self.ParseGenSection     ( line )
                    elif line.startswith ( "&guess" ) or line.startswith ( "$guess" ): self.ParseGuessSection   ( line )
                    elif line.startswith ( "&hess"  ) or line.startswith ( "$hess"  ): self.ParseHessianSection ( line )
                    elif line.startswith ( "&zmat"  ) or line.startswith ( "$zmat"  ): self.ParseZmatSection    ( line )
            except:
                pass
            # . Close the file.
            self.WarningStop ( )
            self.Close ( )
            # . Set the parsed flag and some other options.
            self.log     = None
            self.isParsed = True

    def ParseGenSection ( self, line ):
        """Parse general commands."""
        # . Get the terminator as the opening character of the current line.
        terminator = line[0:1]
        # . Loop until a terminator is reached
        QCONTINUE = True
        QFIRST    = True
        while QCONTINUE:
            if QFIRST:
                items  = line.split ( )
                QFIRST = False
            else:
                items  = self.GetTokens ( )
            for item in items:
                if   item == terminator:           QCONTINUE         = False
                elif item.startswith ( "molchg" ): self.charge       = self.ParseKeywordOption ( item, int )
                elif item.startswith ( "multip" ): self.multiplicity = self.ParseKeywordOption ( item, int )

    def ParseGuessSection ( self, line ):
        """Parse the guess section."""
        # . Get the terminator as the opening character of the current line.
        terminator = line[0:1]
        # . Get the guess basis.
        words = line.split ( None, 1 )
        if len ( words ) > 1:
            tokens = words[1].split ( "=", 1 )
            if ( len ( tokens ) > 1 ) and ( tokens[0] == "basgss" ): guessbasis = tokens[1].strip ( )
        # . Initialization.
        coefficients = None
        nbasis       = None
        orbitals     = None
        orbitalsets  = {}
        key          = ""
        warnings0   = self.warnings
        # . Loop until an end of section is reached
        while True:
            items = self.GetTokens ( )
            # . The end found.
            if ( len ( items ) == 1 ) and ( items[0] == terminator ):
                break
            # . An orbital set header.
            elif ( len ( items ) > 2 ) and ( items[-2] == "Molecular" ) and ( items[-1] == "Orbitals" ):
                key      = " ".join ( items[0:-2] )
                orbitals = None
            # . An orbital line.
            elif ( "Orbital" in items ):
                # . Orbital data.
                if "Energy"     in items: energy    = float ( items[items.index ( "Energy"     ) + 1] )
                else:                     energy    = 0.0
                if "Occupation" in items: occupancy = float ( items[items.index ( "Occupation" ) + 1] )
                else:                     occupancy = 0.0
                if ( nbasis is None ) and ( coefficients is not None ): nbasis = len ( coefficients )
                coefficients = []
                # . Orbital set data.
                if orbitals is None: orbitals = []
                orbitals.append ( ( energy, occupancy, coefficients ) )
                # . Orbital collection data.
                if key is not None:
                    if key in orbitalsets: self.Warning ( "Orbital sets with duplicate names.", False )
                    else:
                        orbitalsets[key] = orbitals
                        key              = None
            # . A coefficient line.
            elif ( len ( items ) > 0 ) and ( coefficients is not None ):
                for item in items: coefficients.append ( float ( item ) )
                if ( nbasis is not None ) and ( len ( coefficients ) > nbasis ): self.Warning ( "There are orbitals with differing numbers of coefficients.", False )
            # . Other lines.
            else:
                self.Warning ( "Unable to interpret guess section line.", False )
        # . Save the data if there have been no warnings.
        if ( warnings0 == self.warnings ) and ( nbasis is not None ):
            # . Initialization.
            self.orbitalsets = {}
            # . Process each of the orbital sets in turn.
            for ( key, orbitals ) in orbitalsets.items ( ):
                norbitals   = len ( orbitals )
                energies    = Array.WithExtent  ( norbitals )
                occupancies = Array.WithExtent  ( norbitals )
                vectors     = Array.WithExtents ( nbasis, norbitals )
                for ( i, ( energy, occupancy, coefficients ) ) in enumerate ( orbitals ):
                    energies[i]    = energy
                    occupancies[i] = occupancy
                    for ( j, v ) in enumerate ( coefficients ): vectors[j,i] = v
                self.orbitalsets[key.lower ( )] = ( norbitals, nbasis, energies, occupancies, vectors )

    def ParseHessianSection ( self, line ):
        """Pass the Hessian section."""
        # . Get the terminator as the opening character of the current line.
        terminator = line[0:1]
        # . Initialization.
        n                 = len ( self.atomicNumbers )
        n3                = 3 * n
        hessian           = Array.WithExtent ( n3, storageType = StorageType.Symmetric )
        numberOfWarnings0 = self.warnings
        numberRead        = 0
        # . Loop until an end of section is reached
        while True:
            items = self.GetTokens ( )
            # . A single item.
            if len ( items ) == 1:
                # . The end.
                if items[0] == terminator:
                    break
                # . Next column block.
                else:
                    columnIncrement = int ( items[0] ) - 1
            # . A data line.
            elif len ( items ) > 1:
                row = int ( items.pop ( 0 ) ) - 1
                for ( c, item ) in enumerate ( items ):
                    hessian[row,columnIncrement+c] = float ( item )
                numberRead += len ( items )
        # . Check for correct number of items.
        if numberRead != ( ( n3 * ( n3 + 1 ) ) // 2 ):
            self.Warning ( "Incorrect number of Hessian elements read.", False )
        # . Save the data if there have been no warnings.
        if numberOfWarnings0 == self.warnings:
            self.hessian = hessian

    def ParseKeywordOption ( self, option, converter ):
        """Parse a keyword option of the form xxx=yyy."""
        QOK    = False
        value  = None
        tokens = option.split ( "=", 1 )
        if len ( tokens ) == 2:
            QOK = True
            if converter is None:
                value = tokens[1]
            else:
                try:    value = converter ( tokens[1] )
                except: QOK   = False
        if not QOK: self.Warning ( "Unable to parse keyword option: " + option + ".", False )
        return value

    def ParseZmatSection ( self, line ):
        """Parse the geometry (XYZ only for the moment)."""
        # . Get the terminator as the opening character of the current line.
        terminator = line[0:1]
        # . Initialization.
        atomicNumbers = []
        xyz           = []
        warnings0    = self.warnings
        # . Loop until an end of section is reached.
        while True:
            items = self.GetTokens ( converters = [ None, JaguarInputFileReader.ToFloat, JaguarInputFileReader.ToFloat, JaguarInputFileReader.ToFloat ] )
            # . The end found.
            if ( len ( items ) == 1 ) and ( items[0] == terminator ):
                break
            # . An XYZ line.
            elif len ( items ) == 4:
                atomicNumber = PeriodicTable.AtomicNumber ( items[0] )
                if atomicNumber > 0:
                    atomicNumbers.append ( atomicNumber )
                    xyz.append ( ( items[1], items[2], items[3] ) )
            # . Other lines.
            else:
                self.Warning ( "Unable to interpret geometry line.", True )
        # . Save the data if there have been no warnings.
        if warnings0 == self.warnings:
            self.atomicNumbers = atomicNumbers
            self.coordinates3  = Coordinates3.WithExtent ( len ( xyz ) )
            for ( i, ( x, y, z ) ) in enumerate ( xyz ):
                self.coordinates3[i,0] = x
                self.coordinates3[i,1] = y
                self.coordinates3[i,2] = z

    @classmethod
    def PathToCoordinates3 ( selfClass, path, log = logFile ):
        """Return the coordinates from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse ( log = log )
        return inFile.ToCoordinates3 ( )

    @classmethod
    def PathToSystem ( selfClass, path, log = logFile ):
        """Return the system from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse ( log = log )
        return inFile.ToSystem ( )

    def ToCoordinates3 ( self ):
        """Return a coordinates3 object."""
        if self.isParsed and hasattr ( self, "coordinates3" ): return self.coordinates3
        else:                                                 return None

    def ToElectronicState ( self ):
        """Return an electronic state object."""
        if self.isParsed:
            return { "Charge"       : getattr ( self, "charge"       , 0 ) ,
                     "Multiplicity" : getattr ( self, "multiplicity" , 1 ) }
        else:
            return None

    def ToSystem ( self ):
        """Return a system."""
        system = None
        if self.isParsed and hasattr ( self, "atomicNumbers" ):
            system              = System.FromAtoms ( self.atomicNumbers )
            system.coordinates3 = self.ToCoordinates3 ( )
            system.Set ( "Electronic State", self.ToElectronicState ( ), attribute = "electronicState" )
            if hasattr ( self, "entryName" ): system.label = self.entryName
        return system

#===================================================================================================================================
# . Importer definitions.
#===================================================================================================================================
_Importer.AddHandler ( { Coordinates3 : JaguarInputFileReader.PathToCoordinates3 ,
                         System       : JaguarInputFileReader.PathToSystem       } ,
                       [ "jagin", "jin", "JAGIN", "JIN" ], "Jaguar Input", defaultFunction = JaguarInputFileReader.PathToSystem )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
