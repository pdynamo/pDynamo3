"""Classes and functions for writing DYL files."""

# . DYL files are pDynamo YAML files for single systems.
# . Here they are written explicitly but they can also be dumped using YAML.

from  pCore                 import TextFileWriter, TextFileWriterError
from  pMolecule             import BondType, System
from  pScientific.Geometry3 import Coordinates3
from .ExportImport          import _Exporter

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Header and footer.
_Header = "---\n# !DYL File\n"
_Footer = "...\n"

# . Atom and bond fields with and without indices.
_AtomFieldI  = [ "Index" ]
_AtomFields  = [ "Label", "Atomic Number", "X", "Y", "Z" ]
_BondFieldsI = [ "Index 1", "Index 2", "Bond Type" ]
_BondFieldsL = [ "Label 1", "Label 2", "Bond Type" ] 

# . Other options.
_JoinFormat = " , "
_XYZFormat  = "{:10.5f}"

#===================================================================================================================================
# . DYL file writer class.
#===================================================================================================================================
class DYLFileWriter ( TextFileWriter ):
    """DYLFileWriter is the class for DYL files that are to be written."""

    _classLabel = "DYL File Writer"

    @classmethod
    def PathFromSystem ( selfClass, path, system, **options ):
        """Create a file given a path and system."""
        outFile = selfClass.FromPath ( path )
        outFile.WriteSingleSystem ( system, **options )

    # . Options with defaults include:
    #
    #   - documentation   None
    #   - label           None
    #   - useIndices      False
    #   - xyz             None
    #   - xyzFormat       _XYZFormat
    #
    def WriteFrame ( self, system, **options ):
        """Write a single frame."""
        # . Options.
        useIndices = options.get ( "useIndices", False               )
        xyz        = options.get ( "xyz"       , system.coordinates3 )
        xyzFormat  = options.get ( "xyzFormat" , _XYZFormat          )
        # . Index and label arrays.
        length = max ( [ len ( atom.path ) for atom in system.atoms ] )
        labels = [ atom.path.ljust ( length ) for atom in system.atoms ]
        if useIndices:
            indices = [ "{:d}".format ( i+1 ) for i in range ( len ( system.atoms ) ) ] # . Start from 1.
            length  = len ( indices[-1] )
            indices = [ "{:s}".format ( index.rjust ( length ) ) for index in indices ]
        # . Header.
        self.file.write ( _Header )
        # . Label and documentation.
        label = options.get ( "label", None )
        if ( label is     None ) or  ( len ( label ) <= 0 ): label = system.label
        if ( label is not None ) and ( len ( label ) >  0 ):
            self.file.write ( "Label: {:s}\n".format ( label ) )
        documentation = options.get ( "documentation", None )
        if documentation is not None:
            self.file.write ( "Documentation: >-\n{:s}\n".format ( documentation ) ) # . Needs checking (e.g. for indentation).
        # . Atoms.
        atomFields  = _AtomFields
        atomItems   = []
        for ( a, atom ) in enumerate ( system.atoms ):
            atomItems.append ( [ labels[a]                            ,
                                 "{:3d}".format ( atom.atomicNumber ) ,
                                 xyzFormat.format ( xyz[a,0] )        ,
                                 xyzFormat.format ( xyz[a,1] )        ,
                                 xyzFormat.format ( xyz[a,2] )      ] )
        if useIndices:
            atomFields = _AtomFieldI + atomFields
            for ( items, index ) in zip ( atomItems, indices ):
                items.insert ( 0, index )
        self.file.write ( "Atom Fields:\n" )
        for field in atomFields:
            self.file.write ( "  - {:s}\n".format ( field ) )
        self.file.write ( "Atoms:\n" )
        for items in atomItems:
            self.file.write ( "  - [ {:s} ]\n".format ( _JoinFormat.join ( items ) ) )
        # . Bonds.
        bonds   = system.connectivity.bonds
        mapping = system.connectivity.nodeIndices
        if ( bonds is not None ) and ( len ( bonds ) > 0 ):
            length    = max ( [ len ( bond.type.label ) for bond in bonds ] )
            types     = [ bond.type.label.ljust ( length ) for bond in bonds ]
            bondItems = []
            if useIndices:
                bondFields = _BondFieldsI
                for ( b, bond ) in enumerate ( bonds ):
                    bondItems.append ( [ indices[mapping[bond.node1]], indices[mapping[bond.node2]], types[b] ] )
            else:
                bondFields = _BondFieldsL
                for ( b, bond ) in enumerate ( bonds ):
                    bondItems.append ( [ labels[mapping[bond.node1]], labels[mapping[bond.node2]], types[b] ] )
            self.file.write ( "Bond Fields:\n" )
            for field in bondFields:
                self.file.write ( "  - {:s}\n".format ( field ) )
            self.file.write ( "Bonds:\n" )
            for items in bondItems:
                self.file.write ( "  - [ {:s} ]\n".format ( _JoinFormat.join ( items ) ) )
        # . Other properties.
        # . Formal charges.
        formalCharges = [ ( i, atom.formalCharge ) for ( i, atom ) in enumerate ( system.atoms ) if ( atom.formalCharge != 0 ) ]
        if len ( formalCharges ) > 0:
            self.file.write ( "Formal Charges:\n" )
            if useIndices: items = indices
            else:          items = labels
            for ( i, q ) in formalCharges:
                self.file.write ( "  {:s}: {:4d}\n".format ( items[i], q ) )
        # . Footer.
        self.file.write ( _Footer )

    def WriteSingleSystem ( self, system, **options ):
        """Write a complete file."""
        self.Open  ( )
        self.WriteFrame ( system, **options )
        self.Close ( )

#===================================================================================================================================
# . Exporter definitions.
#===================================================================================================================================
_Exporter.AddHandler ( { System : DYLFileWriter.PathFromSystem } , [ "dyl", "DYL" ], "pDynamo YAML" )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
