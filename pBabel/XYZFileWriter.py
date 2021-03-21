"""Classes and functions for writing XYZ files."""

import os

from  pCore                 import TextFileWriter      , \
                                   TextFileWriterError
from  pMolecule             import System
from  pScientific           import PeriodicTable
from  pScientific.Geometry3 import Coordinates3
from .ExportImport          import _Exporter

# . Maximum label length (needed to avoid errors in some programs - e.g. VMD).
_MaximumLabelLength = 80

#===================================================================================================================================
# . XYZ file writer class.
#===================================================================================================================================
class XYZFileWriter ( TextFileWriter ):
    """XYZFileWriter is the class for XYZ files that are to be written."""

    _classLabel = "XYZ File Writer"

    @classmethod
    def PathFromSystem ( selfClass, path, system, label = None, xyz = None ):
        """Create a file given a path and system."""
        outFile = selfClass.FromPath ( path )
        outFile.WriteSingleSystem ( system, label = label, xyz = xyz )

    def WriteFrame ( self, system, label = None, xyz = None ):
        """Write a single frame."""
        # . Check the data.
        if xyz == None: xyz = system.coordinates3
        if ( xyz is None ) or ( not isinstance ( xyz, Coordinates3 ) ) or ( xyz.rows != len ( system.atoms ) ):
            print ( len ( system.atoms ), xyz.rows )
            raise TextFileWriterError ( "Invalid or missing data to write to XYZ file." )
        # . Check the label.
        if label is None: label = system.label
        # . Write the frame.
        self.file.write ( "{:6d}\n".format ( len ( system.atoms ) ) )
        if label is None:
            self.file.write ( "\n" )
        else:
            if len ( label ) > _MaximumLabelLength: label = label[0:_MaximumLabelLength]
            self.file.write ( label + "\n" )
        for ( i, atom ) in enumerate ( system.atoms ):
            symbol = PeriodicTable.Symbol ( atom.atomicNumber )
            self.file.write ( "{:<5s}{:25.15f}{:25.15f}{:25.15f}\n".format ( symbol, xyz[i,0], xyz[i,1], xyz[i,2] ) )

    def WriteSingleSystem ( self, system, label = None, xyz = None ):
        """Write a complete file."""
        self.Open  ( )
        self.WriteFrame ( system, label = label, xyz = xyz )
        self.Close ( )

#===================================================================================================================================
# . Exporter definitions.
#===================================================================================================================================
_Exporter.AddHandler ( { System : XYZFileWriter.PathFromSystem } , [ "xyz" ], "XYZ" )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
