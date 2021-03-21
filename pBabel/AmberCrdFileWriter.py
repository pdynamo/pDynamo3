"""Classes and functions for writing Amber crd files."""

from  pCore                 import TextFileWriter
from  pMolecule             import System
from  pScientific.Geometry3 import Coordinates3
from .ExportImport          import _Exporter

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . Number of atoms line (second line).
_AtomLineFormat = "{:6d}\n"

# . Coordinate lines.
_AtomCoordinateFormat = "{:12.7f}{:12.7f}{:12.7f}"
_NumberAtomsPerLine   = 2

# . Symmetry line.
_SymmetryLineFormat   = "{:12.7f}{:12.7f}{:12.7f}"

#===================================================================================================================================
# . AmberCrd file writer class.
#===================================================================================================================================
class AmberCrdFileWriter ( TextFileWriter ):
    """AmberCrdFileWriter is the class for Amber crd files that are to be written."""

    _classLabel = "Amber Crd File Writer"

    @classmethod
    def PathFromSystem ( selfClass, path, system, includeSymmetryParameters = True, label = None, xyz = None ):
        """Create a file given a path and system."""
        outFile = selfClass.FromPath ( path )
        outFile.WriteSystem ( system, includeSymmetryParameters = includeSymmetryParameters, label = label, xyz = xyz )

    def WriteFrame ( self, system, includeSymmetryParameters = None, label = None, xyz = None ):
        """Write a single frame."""
        # . Check the data.
        if xyz == None: xyz = system.coordinates3
        if ( xyz is None ) or ( not isinstance ( xyz, Coordinates3 ) ) or ( xyz.rows != len ( system.atoms ) ): raise TextFileWriterError ( "Invalid or missing data to write to Amber crd file." )
        # . Check the label.
        if label is None: label = system.label
        # . Title.
        self.file.write ( label + "\n" )
        # . Number of atoms.
        self.file.write ( _AtomLineFormat.format ( len ( system.atoms ) ) )
        # . Coordinates.
        for i in range ( len ( system.atoms ) ):
            self.write ( _AtomCoordinateFormat.format ( xyz[i,0], xyz[i,1], xyz[i,2] ) )
            if ( ( i + 1 ) % _NumberAtomsPerLine ) == 0: self.write ( "\n" )
        if ( ( len ( system.atoms ) + 1 ) % _NumberAtomsPerLine ) == 0: self.write ( "\n" )
        # . Symmetry.
        if includeSymmetryParameters:
            symmetryParameters = system.symmetryParameters
            if symmetryParameters is not None:
                self.write ( _SymmetryLineFormat.format ( symmetryParameters.a, symmetryParameters.b, symmetryParameters.c ) )

    def WriteSystem ( self, system, includeSymmetryParameters = True, label = None, xyz = None ):
        """Write a system."""
        self.Open  ( )
        self.WriteFrame ( system, includeSymmetryParameters = includeSymmetryParameters, label = label, xyz = xyz )
        self.Close ( )

#===================================================================================================================================
# . Exporter definitions.
#===================================================================================================================================
_Exporter.AddHandler ( { System : AmberCrdFileWriter.PathFromSystem } , [ "crd", "CRD" ], "Amber Coordinates" )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
