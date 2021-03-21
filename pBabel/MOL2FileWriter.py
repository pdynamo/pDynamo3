"""Classes and functions for writing MOL2 files."""

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
# . Atom type.
_AtomType = "Any"

# . The charge type.
_ChargeType = "NO_CHARGES"

# . The maximum line length.
_MaximumLineLength = 80

# . Bond type definitions.
_MOL2BondTypes = { ( BondType.Single    , False ) : "1"  ,
                   ( BondType.Double    , False ) : "2"  ,
                   ( BondType.Triple    , False ) : "3"  ,
                   ( BondType.Single    , True  ) : "ar" ,
                   ( BondType.Double    , True  ) : "ar" ,
                   ( BondType.Undefined , False ) : "un" }

# . The molecule type.
_MoleculeType = "SMALL"

# . Section header.
_RTIString  = "@<TRIPOS>"

#===================================================================================================================================
# . MOL2 file writer class.
#===================================================================================================================================
class MOL2FileWriter ( TextFileWriter ):
    """MOL2FileWriter is the class for MOL2 files that are to be written."""

    _classLabel = "MOL2 File Writer"

    @classmethod
    def PathFromSystem ( selfClass, path, system, label = None, useKekule = False, xyz = None ):
        """Create a file given a path and system."""
        outFile = selfClass.FromPath ( path )
        outFile.WriteSingleSystem ( system, label = label, useKekule = useKekule, xyz = xyz )

    def WriteFrame ( self, system, label = None, useAtomicNumberTypes = True, useKekule = False, xyz = None ):
        """Write a single frame."""
        # . Get data.
        natoms = len ( system.atoms )
        nbonds = len ( system.connectivity.bonds )
        # . Check for Kekule output.
        useKekule = CheckForKekuleOutput ( system.connectivity, useKekule )
        # . Check the data.
        if xyz == None: xyz = system.coordinates3
        if ( xyz is None ) or ( not isinstance ( xyz, Coordinates3 ) ) or ( xyz.rows != natoms ): raise TextFileWriterError ( "Invalid or missing data to write to MOL2 file." )
        # . Check the label.
        if label is None: label = system.label
        # . Molecule block.
        self.file.write ( _RTIString + "MOLECULE\n" )
        if label is None: self.file.write ( "\n" )
        else:             self.file.write ( label.strip ( )[0:min ( len ( label ), _MaximumLineLength )] + "\n" )
        self.file.write ( "{:d} {:d}\n".format ( natoms, nbonds ) )
        self.file.write ( _MoleculeType + "\n" )
        self.file.write ( _ChargeType   + "\n" )
        self.file.write ( "\n" )
        # . Atom block.
        self.file.write ( _RTIString + "ATOM\n" )
        for ( i, atom ) in enumerate ( system.atoms ):
            if useAtomicNumberTypes: sType = PeriodicTable.Symbol ( atom.atomicNumber )
            else:                    sType = _AtomType
            self.file.write ( "{:6d} {:<6s} {:10.4f} {:10.4f} {:10.4f}   {:s}\n".format ( i+1, atom.label, xyz[i,0], xyz[i,1], xyz[i,2], sType ) )
        self.file.write ( "\n" )
        # . Bond block.
        mapping = system.connectivity.nodeIndices
        self.file.write ( _RTIString + "BOND\n" )
        for ( i, bond ) in enumerate ( system.connectivity.bonds ):
            if useKekule: bondKey = ( bond.type, False           )
            else:         bondKey = ( bond.type, bond.isAromatic )
            self.file.write ( "{:6d} {:6d} {:6d} {:<s}\n".format ( i+1, mapping[bond.node1] + 1, mapping[bond.node2] + 1, _MOL2BondTypes[bondKey] ) )
        self.file.write ( "\n" )
        # . Symmetry block.
        if system.symmetryParameters is not None:
            p = system.symmetryParameters
            self.file.write ( _RTIString + "CRYSIN\n" )
            self.file.write ( "{:.3f} {:.3f} {:.3f} {:.1f} {:.1f} {:.1f}  1  1\n".format ( p.a, p.b, p.c, p.alpha, p.beta, p.gamma ) )

    def WriteSingleSystem ( self, system, label = None, useKekule = False, xyz = None ):
        """Write a complete file."""
        self.Open  ( )
        self.WriteFrame ( system, label = label, useKekule = useKekule, xyz = xyz )
        self.Close ( )

#===================================================================================================================================
# . Exporter definitions.
#===================================================================================================================================
_Exporter.AddHandler ( { System : MOL2FileWriter.PathFromSystem } , [ "mol2", "MOL2" ], "Tripos MOL2" )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
