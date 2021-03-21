"""Read data from a CHARMM CRD file."""

# . It may be good to have an extract method here for coordinates (i.e. given sequence).

from  pCore                 import Clone         , \
                                   logFile       , \
                                   LogFileActive , \
                                   TextFileReader
from  pMolecule             import Sequence      , \
                                   System
from  pScientific.Geometry3 import Coordinates3
from .ExportImport          import _Importer

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class CHARMMCRDFileReader ( TextFileReader ):
    """Class for reading CHARMM CRD files."""

    _classLabel = "CHARMM CRD File Reader"

    def GetAtomLineFormat ( self, extendedFormat ):
        """Get the format of the atom lines.

        This is:
        ATOMNO RESNO   RES  TYPE  X     Y     Z    SEGID RESID Weighting
        I5    I5  1X A4 1X A4 F10.5  F10.5  F10.5  1X A4 1X A4 F10.5
        I10   I10 2X A8 2X A8 F20.10 F20.10 F20.10 2X A8 2X A8 F20.10
        """
        a = 4 ; f = 10 ; i = 5 ; x = 1
        if extendedFormat:
            a *= 2 ; f *= 2 ; i *= 2 ; x *= 2
        p = i
        format =      [ ( p, p+i, int  , 0   ) ] ; p += (i+x) # . Res. number.
        format.append ( ( p, p+a, None , ""  ) ) ; p += (a+x) # . Res. name.
        format.append ( ( p, p+a, None , ""  ) ) ; p +=  a    # . Atom name.
        format.append ( ( p, p+f, float, 0.0 ) ) ; p +=  f    # . X.
        format.append ( ( p, p+f, float, 0.0 ) ) ; p +=  f    # . Y.
        format.append ( ( p, p+f, float, 0.0 ) ) ; p += (f+x) # . Z.
        format.append ( ( p, p+a, None , ""  ) ) ; p += (a+x) # . Seg. ID.
        format.append ( ( p, p+a, None , ""  ) ) ; p +=  a    # . Res. ID.
        format.append ( ( p, p+f, float, 0.0 ) )              # . Weighting.
        return format

    def GetLine ( self, signalWarnings = True ):
        """Get a non-empty line removed of comments."""
        try:
            QFOUND = False
            while not QFOUND:
                line  = next ( self.file )
                index = line.find ( "!" )
                if index >= 0: line = line[:index]
                line  = line.strip ( )
                self.linesParsed += 1
                QFOUND = ( len ( line ) > 0 )
            return line
        except:
            if signalWarnings: self.Warning ( "Unexpected end-of-file.", True )
            raise EOFError

    def Parse ( self, log = logFile ):
        """Parsing."""
        if not self.isParsed:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
            # . Initialization.
            extendedFormat = False
            natoms         =  0
            self.atoms     = []
            self.title     = []
            # . Open the file.
            self.Open ( )
            # . Parse all the lines.
            try:
                # . Title and number of atoms lines.
                while True:
                    line = self.GetLine ( ).strip ( )
                    if line.startswith ( "*" ):
                        line = line[1:].strip ( )
                        if len ( line ) > 0: self.title.append ( line )
                    else:
                        tokens = line.split ( )
                        try:
                            natoms = int ( tokens[0] )
                            if len ( tokens ) > 1:
                                if tokens[1].upper ( ) == "EXT": extendedFormat = True
                                else: raise
                        except:
                            self.Warning ( "Invalid atom counter line.", True )
                        break
                # . Get the atom line format.
                atomlineformat = self.GetAtomLineFormat ( extendedFormat )
                # . Atom lines - store atom name, residue name, resID, segID, x, y, z, w.
                for i in range ( natoms ):
                    tokens = self.GetFixedFormatTokens ( *atomlineformat )
                    self.atoms.append ( ( tokens[2], tokens[1], tokens[7], tokens[6], tokens[3], tokens[4], tokens[5], tokens[8] ) )
            except EOFError:
                pass
            # . Close the file.
            self.WarningStop ( )
            self.Close ( )
            # . Set the parsed flag and some other options.
            self.log     = None
            self.isParsed = True

    @classmethod
    def PathToCoordinates3 ( selfClass, path, log = logFile ):
        """Return the coordinates from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse   ( log = log )
        inFile.Summary ( log = log )
        return inFile.ToCoordinates3 ( )

    @classmethod
    def PathToSequence ( selfClass, path, log = logFile ):
        """Return the sequence from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse   ( log = log )
        inFile.Summary ( log = log )
        return inFile.ToSequence ( )

    @classmethod
    def PathToSystem ( selfClass, path, log = logFile ):
        """Return the system from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse   ( log = log )
        inFile.Summary ( log = log )
        return inFile.ToSystem ( )

    def SummaryItems ( self ):
        """Summary items."""
        items = []
        if self.isParsed: items.append ( ( "Atom Lines", "{:d}".format ( len ( self.atoms ) ) ) )
        return items

    def ToCoordinates3 ( self ):
        """Return the coordinates."""
        if self.isParsed:
            coordinates3 = Coordinates3.WithExtent ( len ( self.atoms ) )
            for ( i, datum ) in enumerate ( self.atoms ):
                coordinates3[i,0] = datum[4]
                coordinates3[i,1] = datum[5]
                coordinates3[i,2] = datum[6]
            return coordinates3
        else:
            return None

    def ToSequence ( self ):
        """Return a sequence."""
        if self.isParsed:
            # . Initialization.
            majorSeparator = Sequence._attributable["labelSeparator"]
            minorSeparator = Sequence._attributable["fieldSeparator"]
            # . Create atom paths - segID (datum[3]), resname (datum[1]), resID (datum[2]), atom name (datum[0]).
            atomPaths = []
            for datum in self.atoms:
                items = []
                for i in range ( 4 ):
                    item = datum[i]
                    if len ( item ) <= 0: item = " "
                    items.append ( item )
                atomPaths.append ( items[3] + majorSeparator + items[1] + minorSeparator + items[2] + majorSeparator + items[0] )
            # . Get the sequence.
            return Sequence.FromAtomPaths ( atomPaths )

    def ToSystem ( self ):
        """Return a system from the CHARMM CRD data."""
        if self.isParsed:
            sequence            = self.ToSequence ( )
            system              = System.FromSequence ( sequence )
            system.coordinates3 = self.ToCoordinates3 ( )
            system.label        = " ".join ( self.title ).strip ( )
            return system
        else:
            return None

#===================================================================================================================================
# . Importer definitions.
#===================================================================================================================================
_Importer.AddHandler ( { Coordinates3 : CHARMMCRDFileReader.PathToCoordinates3 ,
                         System       : CHARMMCRDFileReader.PathToSystem       } ,
                       [ "chm", "CHM" ], "Charmm Coordinates", defaultFunction = CHARMMCRDFileReader.PathToSystem )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass

