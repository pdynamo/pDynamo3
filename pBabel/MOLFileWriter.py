"""Classes and functions for writing MOL files."""

from  pCore                 import TextFileWriter       , \
                                   TextFileWriterError
from  pMolecule             import BondType             , \
                                   CheckForKekuleOutput , \
                                   System
from  pScientific           import PeriodicTable
from  pScientific.Geometry3 import Coordinates3
from .ExportImport          import _Exporter

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . The maximum line length.
_MaximumLineLength = 80

# . Bond type definitions.
_MOLBondTypes = { ( BondType.Single , False ) : 1 ,
                  ( BondType.Double , False ) : 2 ,
                  ( BondType.Triple , False ) : 3 ,
                  ( BondType.Single , True  ) : 4 ,
                  ( BondType.Double , True  ) : 4 }

#===================================================================================================================================
# . MOL file writer class.
#===================================================================================================================================
class MOLFileWriter ( TextFileWriter ):
    """MOLFileWriter is the class for MOL files that are to be written."""

    _classLabel = "MOL File Writer"

    @classmethod
    def PathFromSystem ( selfClass, path, system, label = None, useKekule = False, xyz = None ):
        """Create a file given a path and system."""
        outFile = selfClass.FromPath ( path )
        outFile.WriteSingleSystem ( system, label = label, useKekule = useKekule, xyz = xyz )

    def WriteFrame ( self, system, label = None, useKekule = False, xyz = None ):
        """Write a single frame."""
        # . Check the data.
        useKekule = CheckForKekuleOutput ( system.connectivity, useKekule )
        if xyz == None: xyz = system.coordinates3
        if ( xyz is None ) or ( not isinstance ( xyz, Coordinates3 ) ) or ( xyz.rows != len ( system.atoms ) ): raise TextFileWriterError ( "Invalid or missing data to write to MOL file." )
        # . Check the label.
        if label is None: label = system.label
        # . Write the header block.
        if label is None: self.file.write ( "\n\n\n" )
        else:             self.file.write ( label.strip ( )[0:min ( len ( label ), _MaximumLineLength )] + "\n\n\n" )
        # . Counts line.
        self.file.write ( "{:3d}{:3d}".format ( len ( system.atoms ), len ( system.connectivity.bonds ) ) + 8 * "  0" + "999 V2000\n" )
        # . Atom block.
        for ( i, atom ) in enumerate ( system.atoms ):
            self.file.write ( "{:10.4f}{:10.4f}{:10.4f} {:<2s}".format ( xyz[i,0], xyz[i,1], xyz[i,2], PeriodicTable.Symbol ( atom.atomicNumber ) ) + 5 * "  0" + "\n" )
        # . Bond block.
        mapping = system.connectivity.nodeIndices
        for ( i, bond ) in enumerate ( system.connectivity.bonds ):
            if useKekule: bondKey = ( bond.type, False           )
            else:         bondKey = ( bond.type, bond.isAromatic )
            self.file.write ( "{:3d}{:3d}{:3d}".format ( mapping[bond.node1] + 1, mapping[bond.node2] + 1, _MOLBondTypes[bondKey] ) + 3 * "  0" + "\n" )
        # . Properties block.
        for ( i, atom ) in enumerate ( system.atoms ):
            formalCharge = getattr ( atom, "formalCharge", 0 )
            if formalCharge != 0: self.file.write ( "M  CHG  1{:4d}{:4d}\n".format ( i+1, formalCharge ) )
        self.file.write ( "M  END\n" )

    def WriteSingleSystem ( self, system, label = None, useKekule = False, xyz = None ):
        """Write a complete file."""
        self.Open  ( )
        self.WriteFrame ( system, label = label, useKekule = useKekule, xyz = xyz )
        self.Close ( )

#===================================================================================================================================
# . Exporter definitions.
#===================================================================================================================================
_Exporter.AddHandler ( { System : MOLFileWriter.PathFromSystem } , [ "mol", "MOL" ], "MDL MOL" )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
