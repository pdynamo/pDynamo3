"""Read data from an fDynamo CRD file."""

from  pCore                 import Clone, logFile, LogFileActive, TextFileReader
from  pMolecule             import Sequence, System
from  pScientific.Geometry3 import Coordinates3
from .ExportImport          import _Importer

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class fDynamoCRDFileReader ( TextFileReader ):
    """Class for reading fDynamo CRD files."""

    _classLabel = "fDynamo CRD File Reader"

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
            # . Read data.
            self.atoms      = []
            self.residues   = []
            self.subsystems = []
            # . Open the file.
            self.Open ( )
            # . Parse all the lines.
            try:
                ( natoms, nresidues, nsubsystems ) = self.__GetCounters ( )
                # . Subsystems.
                for isub in range ( nsubsystems ):
                    nres = self.__SubsystemHeader ( )
                    for ires in range ( nres ):
                        natm = self.__ResidueHeader ( )
                        for iatm in range ( natm ):
                            self.__AtomRecord ( )
                # . Do a final logic check.
                if ( len ( self.atoms ) != natoms ) or ( len ( self.residues ) != nresidues ) or ( len ( self.subsystems ) != nsubsystems ):
                    self.Warning ( "Counter mismatch after parsing: {:d}/{:d}, {:d}/{:d}, {:d}/{:d}.".format ( len ( self.atoms      ) , natoms      , \
                                                                                                               len ( self.residues   ) , nresidues   , \
                                                                                                               len ( self.subsystems ) , nsubsystems ), True )
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
        if self.isParsed:
            items.extend ( [ ( "Atoms",      "{:d}".format ( len ( self.atoms      ) ) ) ,
                             ( "Residues",   "{:d}".format ( len ( self.residues   ) ) ) ,
                             ( "Subsystems", "{:d}".format ( len ( self.subsystems ) ) ) ] )
        return items

    def ToCoordinates3 ( self ):
        """Return the coordinates."""
        if self.isParsed:
            coordinates3 = Coordinates3.WithExtent ( len ( self.atoms ) )
            for ( i, ( aname, atomicNumber, x, y, z ) ) in enumerate ( self.atoms ):
                coordinates3[i,0] = x
                coordinates3[i,1] = y
                coordinates3[i,2] = z
            return coordinates3
        else:
            return None

    def ToSequence ( self ):
        """Return a sequence."""
        sequence = None
        if self.isParsed:
            majorSeparator = Sequence._attributable["labelSeparator"]
            minorSeparator = Sequence._attributable["fieldSeparator"]
            atomPaths = []
            for ( sname, sresidues ) in self.subsystems:
                for ires in range ( sresidues ):
                    ( rname, ratoms ) = residues.pop ( 0 )
                    pathHead = sname + majorSeparator + rname + minorSeparator + "{:d}".format ( ires + 1 ) + majorSeparator
                    for iatom in range ( ratoms ):
                        atom    = atoms.pop ( 0 )
                        atomPaths.append ( pathHead + atom )
            sequence = Sequence.FromAtomPaths ( atomPaths )
        return sequence

    def ToSystem ( self ):
        """Return a system from the CRD data but without sequence information."""
        if self.isParsed:
            atoms = []
            for ( aname, atomicNumber, x, y, z ) in self.atoms:
                atoms.append ( atomicNumber )
            system = System.FromAtoms ( atoms )
            system.coordinates3 = self.ToCoordinates3 ( )
            return system
        else:
            return None

    # . Private methods.
    def __AtomRecord ( self ):
        """Extract atom data from the next record."""
        tokens = self.GetTokens ( converters = [ int, None, int, float, float, float ] )
        if ( len ( tokens ) >= 6 ):
            self.atoms.append ( ( tokens[1], tokens[2], tokens[3], tokens[4], tokens[5] ) )
        else:
            self.Warning ( "Invalid atom data line." )
            self.atoms.append ( ( "", -1, 0.0, 0.0, 0.0 ) )

    def __GetCounters ( self ):
        """Get the counters from the first line of the file."""
        tokens = self.GetTokens ( converters = [ int, int, int ] )
        if ( len ( tokens ) >= 3 ) and ( tokens[0] >= 0 ) and ( tokens[1] >= 0 ) and ( tokens[2] >= 0 ):
            return ( tokens[0], tokens[1], tokens[2] ) # . natoms, nresidues, nsubsystems.
        else:
            self.Warning ( "Invalid natoms, nresidues or nsubsystems counter.", True )
            return ( 0, 0, 0 )

    def __ResidueHeader ( self ):
        """Parse a residue header."""
        # . Initialization.
        name = ""
        natm = ""
        # . Residue line.
        tokens = self.GetTokens ( converters = [ None, int, None ] )
        if ( len ( tokens ) >= 3 ) and ( tokens[0].upper ( ) == "RESIDUE" ): name = tokens[2]
        else: self.Warning ( "Invalid residue header line.", True )
        # . Atom counter line.
        tokens = self.GetTokens ( converters = [ int ] )
        if ( len ( tokens ) >= 1 ) and ( tokens[0] >= 0 ): natm = tokens[0]
        else: self.Warning ( "Invalid residue atom counter.", True )
        self.residues.append ( ( name, natm ) )
        return natm

    def __SubsystemHeader ( self ):
        """Parse a subsystem header."""
        # . Initialization.
        name = ""
        nres = ""
        # . Subsystem line.
        tokens = self.GetTokens ( converters = [ None, int, None ] )
        if ( len ( tokens ) >= 3 ) and ( tokens[0].upper ( ) == "SUBSYSTEM" ): name = tokens[2]
        else: self.Warning ( "Invalid subsystem header line.", True )
        # . Residue counter line.
        tokens = self.GetTokens ( converters = [ int ] )
        if ( len ( tokens ) >= 1 ) and ( tokens[0] >= 0 ): nres = tokens[0]
        else: self.Warning ( "Invalid subsystem residue counter.", True )
        self.subsystems.append ( ( name, nres ) )
        return nres

#===================================================================================================================================
# . Importer definitions.
#===================================================================================================================================
_Importer.AddHandler ( { Coordinates3 : fDynamoCRDFileReader.PathToCoordinates3 ,
                         System       : fDynamoCRDFileReader.PathToSystem       } ,
                       [ "fcrd", "FCRD" ], "Fortran Dynamo Coordinates", defaultFunction = fDynamoCRDFileReader.PathToSystem )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass

