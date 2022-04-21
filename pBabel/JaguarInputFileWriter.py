"""Classes and functions for writing Jaguar input files."""

from  pCore                 import TextFileWriter      , \
                                   TextFileWriterError
from  pMolecule             import System
from  pScientific           import PeriodicTable
from  pScientific.Geometry3 import Coordinates3
from .ExportImport          import _Exporter

#===================================================================================================================================
# . Jaguar input file writer class.
#===================================================================================================================================
class JaguarInputFileWriter ( TextFileWriter ):
    """JaguarInputFileWriter is the class for Jaguar input files that are to be written."""

    _classLabel = "Jaguar Input File Writer"

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
            raise TextFileWriterError ( "Invalid or missing data to write to Jaguar input file." )
        # . Check the label.
        if label is None: label = system.label
        if label is None: label = "Jaguar input file."
        # . Check the electronic state.
        if system.qcModel is None:
            charge = 0
            multip = 1
        else:
            charge = system.qcModel.charge
            multip = system.qcModel.multiplicity
        # . Header.
        self.file.write ( label + "\n\n" )
        # . General section.
        self.file.write ( "&gen\n" )
        self.file.write ( "molchg={:d}\n".format ( charge ) )
        self.file.write ( "multip={:d}\n".format ( multip ) )
        self.file.write ( "&\n" )
        # . Coordinates.
        self.file.write ( "&zmat\n" )
        for ( i, atom ) in enumerate ( system.atoms ):
            symbol = PeriodicTable.Symbol ( atom.atomicNumber, index = i+1 )
            self.file.write ( "{:<5s}{:25.15f}{:25.15f}{:25.15f}\n".format ( symbol, xyz[i,0], xyz[i,1], xyz[i,2] ) )
        self.file.write ( "&\n" )

    def WriteSingleSystem ( self, system, label = None, xyz = None ):
        """Write a complete file."""
        self.Open  ( )
        self.WriteFrame ( system, label = label, xyz = xyz )
        self.Close ( )

#===================================================================================================================================
# . Exporter definitions.
#===================================================================================================================================
_Exporter.AddHandler ( { System : JaguarInputFileWriter.PathFromSystem } , [ "jagin", "jin", "JAGIN", "JIN" ], "Jaguar Input" )

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass
