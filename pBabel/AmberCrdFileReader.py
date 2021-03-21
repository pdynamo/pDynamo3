"""Classes and functions for reading Amber crd files."""

from  pCore                 import logFile            , \
                                   LogFileActive      , \
                                   TextFileReader
from  pScientific.Geometry3 import Coordinates3
from  pScientific.Symmetry  import SymmetryParameters
from .ExportImport          import _Importer

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . The width of the natoms field (second line).
_NATOMSWIDTH = 6

# . The number and the widths of XYZ elements on a line.
_XYZNUMBER =  6
_XYZWIDTH  = 12

# . The symmetry width.
_SYMMETRYNUMBER =  3
_SYMMETRYWIDTH  = 12

#===================================================================================================================================
# . AmberCrd file reader class.
#===================================================================================================================================
class AmberCrdFileReader ( TextFileReader ):
    """AmberCrdFileReader is the class for Amber crd files that are to be read."""

    _classLabel = "Amber Crd File Reader"

    def Parse ( self, log = logFile ):
        """Parse the data on the file."""
        if not self.isParsed:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
            # . Open the file.
            self.Open ( )
            try:
                # . Keyword line.
                self.title = self.GetLine ( )
                # . Number of atoms.
                items    = self.GetTokens ( converters = ( int, ) )
                natoms   = items[0]
                # . The coordinate data.
                items    = self.GetFixedFormatArray ( 3 * natoms, _XYZNUMBER, _XYZWIDTH, converter = float, default = 0.0 )
                self.xyz = Coordinates3.WithExtent ( natoms )
                for n in range ( natoms ):
                    for i in range ( 3 ): self.xyz[n,i] = items[3*n+i]
                # . Symmetry data - optional.
                items = self.GetFixedFormatArray ( 3, _SYMMETRYNUMBER, _SYMMETRYWIDTH, converter = float, default = 0.0, signalWarnings = False )
                self.symmetryParameters = SymmetryParameters ( )
                self.symmetryParameters.SetCrystalParameters ( items[0], items[1], items[2], 90.0, 90.0, 90.0 )
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
    def PathToSymmetryParameters ( selfClass, path, log = logFile ):
        """Return the symmetry paramters from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse   ( log = log )
        inFile.Summary ( log = log )
        return inFile.ToSymmetryParameters ( )

    def SummaryItems ( self ):
        """Summary items."""
        items = []
        if self.isParsed: items.append ( ( "Atoms", "{:d}".format ( self.xyz.rows ) ) )
        return items

    def ToCoordinates3 ( self ):
        """Return the coordinates."""
        if self.isParsed: return self.xyz
        else:            return None

    def ToSymmetryParameters ( self ):
        """Return the symmetry parameters."""
        if self.isParsed: return self.symmetryParameters
        else:            return None

#===================================================================================================================================
# . Importer definitions.
#===================================================================================================================================
_Importer.AddHandler ( { Coordinates3 : AmberCrdFileReader.PathToCoordinates3 } , [ "crd", "CRD" ], "Amber Coordinates" )

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass
