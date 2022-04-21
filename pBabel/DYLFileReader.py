#===================================================================================================================================
# . Classes and functions to read DYL files.
#===================================================================================================================================

from  pCore                 import logFile        , \
                                   LogFileActive  , \
                                   TextFileReader , \
                                   YAMLUnpickle
from  pMolecule             import Atom           , \
                                   Bond           , \
                                   BondType       , \
                                   System
from  pScientific.Geometry3 import Coordinates3
from .ExportImport          import _Importer

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================

#===================================================================================================================================
# . DYL file reader class.
#===================================================================================================================================
class DYLFileReader ( TextFileReader ):
    """DYLFileReader is the class for DYL files that are to be read."""

    _classLabel = "DYL File Reader"

    def Parse ( self, log = logFile ):
        """Parsing."""
        if not self.isParsed:
            # . Read the file.
            if LogFileActive ( log ): self.log = log
            try   : self.data = YAMLUnpickle ( self.path )
            except: self.Warning ( "Error interpreting YAML in DYL file." )
            # . Finish up.
            self.WarningStop ( )
            self.log      = None
            self.isParsed = True

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

    def ToAtoms ( self ):
        """Get a list of atoms."""
        atoms = []
        if self.isParsed:
            # . Atoms.
            atomFields   = self.data["Atom Fields"]
            atomItems    = self.data["Atoms"      ]
            atomicNumber = atomFields.index ( "Atomic Number" )
            label        = atomFields.index ( "Label"         )
            for items in atomItems:
                atoms.append ( Atom.WithOptions ( atomicNumber = items[atomicNumber], label = items[label] ) )
            useIndices = ( "Index" in atomFields )
            if not useIndices:
                labelToIndex = { atom.label : index for ( index, atom ) in enumerate ( atoms ) }
            # . Properties.
            formalCharges = self.data.get ( "Formal Charges", None )
            if ( formalCharges is not None ) and ( len ( formalCharges ) > 0 ):
                if useIndices:
                    for ( key, value ) in formalCharges.items ( ): atoms[             key ].formalCharge = value
                else:
                    for ( key, value ) in formalCharges.items ( ): atoms[labelToIndex[key]].formalCharge = value
        return atoms

    def ToBonds ( self, atoms ):
        """Get a list of atoms."""
        bonds = None
        if self.isParsed:
            bondFields   = self.data.get ( "Bond Fields", None )
            bondItems    = self.data.get ( "Bonds"      , None )
            if ( bondFields is not None ) and ( bondItems is not None ):
                bondType = bondFields.index ( "Bond Type" )
                if ( "Index 1" in bondFields ) and ( "Index 2" in bondFields ):
                    index1 = bondFields.index ( "Index 1" )
                    index2 = bondFields.index ( "Index 2" )
                    bonds  = [ Bond.WithNodes ( atoms[items[index1]], atoms[items[index2]], type = BondType[items[bondType].replace ( " ", "" )] ) for items in bondItems ]
                else:
                    labelToIndex = { atom.label : index for ( index, atom ) in enumerate ( atoms ) }
                    label1       = bondFields.index ( "Label 1" )
                    label2       = bondFields.index ( "Label 2" )
                    bonds        = [ Bond.WithNodes ( atoms[labelToIndex[items[label1]]], atoms[labelToIndex[items[label2]]], type = BondType[items[bondType].replace ( " ", "" )] ) for items in bondItems ]
        return bonds

    def ToCoordinates3 ( self ):
        """Get a set of coordinates."""
        if self.isParsed:
            atomFields = self.data["Atom Fields"]
            atomItems  = self.data["Atoms"      ]
            n   = len ( atomItems )
            x   = atomFields.index ( "X" )
            y   = atomFields.index ( "Y" )
            z   = atomFields.index ( "Z" )
            xyz = Coordinates3.WithExtent ( n )
            for ( i, items ) in enumerate ( atomItems ):
                xyz[i,0] = items[x]
                xyz[i,1] = items[y]
                xyz[i,2] = items[z]
            return xyz
        else:
            return None

    def ToSystem ( self ):
        """Return a system."""
        if self.isParsed:
            atoms = self.ToAtoms (       )
            bonds = self.ToBonds ( atoms )
            system              = System.FromAtoms ( atoms, bonds = bonds )
            system.label        = self.data.get ( "Label", None )
            system.coordinates3 = self.ToCoordinates3 ( )
            return system
        else:
            return None

#===================================================================================================================================
# . Importer definitions.
#===================================================================================================================================
_Importer.AddHandler ( { Coordinates3 : DYLFileReader.PathToCoordinates3 ,
                         System       : DYLFileReader.PathToSystem       } ,
                       [ "dyl", "DYL" ], "pDynamo YAML", defaultFunction = DYLFileReader.PathToSystem )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
