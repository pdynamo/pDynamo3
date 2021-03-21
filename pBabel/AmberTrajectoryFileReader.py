"""Classes and functions for reading Amber Trajectory files."""

from  pCore           import TextFileReader
from .ExportImport    import _Importer
from .TrajectoryMixin import TrajectoryMixin

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . The number and the widths of symmetry elements on a line.
_SymmetryNumber = 3
_SymmetryWidth  = 8

# . The number and the widths of XYZ elements on a line.
_XYZNumber = 10
_XYZWidth  =  8

#===================================================================================================================================
# . AmberTrajectory file reader class.
#===================================================================================================================================
class AmberTrajectoryFileReader ( TextFileReader ):
    """Class for reading Amber Trajectory files."""

    _attributable = dict ( TextFileReader._attributable )
    _classLabel   = "Amber Trajectory File Reader"
    _attributable.update ( { "atomIndices"    : None  ,
                             "freeAtomsOnly"  : True  ,
                             "includeBox"     : False ,
                             "isTrajectory"   : True  , # . A fudge until the mixin is used directly.
                             "numberOfFrames" : -1    , # . Number of frames unknown until the complete trajectory has been read.
                             "numberOfReads"  :  0    ,
                             "owner"          : None  ,
                             "title"          : None  } )

    def __len__ ( self ):
        return self.numberOfFrames

    def _CheckOptions ( self ):
        """Check options."""
        super ( AmberTrajectoryFileReader, self )._CheckOptions ( )
        # . Include a box?
        if hasattr ( self.owner, "symmetry" ) and ( self.owner.symmetry is not None ):
            crystalSystem   = getattr ( self.owner.symmetry, "crystalSystem", None )
            self.includeBox = ( crystalSystem is not None ) and crystalSystem.IsOrthogonal ( )
        # . Atom indices.
        if self.freeAtomsOnly and ( self.owner.freeAtoms is not None ):
            self.atomIndices = self.owner.freeAtoms
        else:
            self.atomIndices = range ( len ( self.owner.atoms ) )

    @classmethod
    def FromPathAndOwner ( selfClass, path, owner, freeAtomsOnly = False ):
        """Constructor given path, owner and other options."""
        self = selfClass.WithOptions ( freeAtomsOnly = freeAtomsOnly ,
                                       owner         = owner         ,
                                       path          = path          )
        self.Open ( )
        return self

    def ReadFooter ( self ):
        """Read the trajectory footer."""
        pass

    def ReadHeader ( self ):
        """Read the trajectory header."""
        self.title = self.GetLine ( )

    def RestoreOwnerData ( self ):
        """Restore data from a frame to the owner."""
        try:
            natoms = len ( self.atomIndices )
            xyz    = self.owner.coordinates3
            items  = self.GetFixedFormatArray ( 3 * natoms, _XYZNumber, _XYZWidth, converter = float, default = 0.0 )
            for ( i, s ) in enumerate ( self.atomIndices ):
                for j in range ( 3 ): xyz[s,j] = items[3*i+j]
            if self.includeBox:
                items = self.GetFixedFormatArray ( 3, _SymmetryNumber, _SymmetryWidth, converter = float, default = 0.0 )
                self.owner.symmetryParameters.SetCrystalParameters ( items[0], items[1], items[2], 90.0, 90.0, 90.0 )
            self.numberOfReads += 1
            return True
        except EOFError:
            self.numberOfFrames = self.numberOfReads
            return False

#===================================================================================================================================
# . Importer definitions.
#===================================================================================================================================
_Importer.AddHandler ( { TrajectoryMixin : AmberTrajectoryFileReader.FromPathAndOwner } , [ "mdcrd", "MDCRD" ], "Amber Trajectory" )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
