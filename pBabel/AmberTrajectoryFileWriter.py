"""Classes and functions for writing Amber Trajectory files.

Taken from the Amber web page:

AMBER trajectory (coordinate or velocity) file specification
This file is optionally written during dynamics in SANDER or GIBBS.

FORMAT(20A4) ITITL
  ITITL  : the title of the current run, from the AMBER
           parameter/topology file

The following snapshot is written every NTWX steps in the trajectory
(specified in the control input file):

FORMAT(10F8.3)  (X(i), Y(i), Z(i), i=1,NATOM)
  X,Y,Z  : coordinates or velocities (velocity units: Angstroms per 1/20.455 ps)

If constant pressure (in 4.1, also constant volume)
For each snapshot:

FORMAT(10F8.3)  BOX(1), BOX(2), BOX(3)
  BOX    : size of periodic box
"""

from  pCore           import TextFileWriter
from .ExportImport    import _Exporter
from .TrajectoryMixin import TrajectoryMixin

#===================================================================================================================================
# . AmberTrajectory file writer class.
#===================================================================================================================================
class AmberTrajectoryFileWriter ( TextFileWriter ):
    """Class for writing Amber Trajectory files."""

    _attributable = dict ( TextFileWriter._attributable )
    _classLabel   = "Amber Trajectory File Writer"
    _attributable.update ( { "atomIndices"    : None  ,
                             "freeAtomsOnly"  : True  ,
                             "includeBox"     : False ,
                             "isTrajectory"   : True  , # . A fudge until the mixin is used directly.
                             "numberOfFrames" : 0     , # . Number of frames currently on trajectory = number of writes.
                             "numberOfWrites" : 0     ,
                             "owner"          : None  ,
                             "title"          : None  } )

    def __len__ ( self ):
        return self.numberOfFrames

    def _CheckOptions ( self ):
        """Check options."""
        super ( AmberTrajectoryFileWriter, self )._CheckOptions ( )
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
    def FromPathAndOwner ( selfClass, path, owner, freeAtomsOnly = False, title = None ):
        """Constructor given path, owner and other options."""
        self = selfClass.WithOptions ( freeAtomsOnly = freeAtomsOnly ,
                                       owner         = owner         ,
                                       path          = path          ,
                                       title         = title         )
        self.Open ( )
        return self

    def WriteFooter ( self ):
        """Write the trajectory footer."""
        pass

    def WriteHeader ( self ):
        """Write the trajectory header."""
        if ( self.title is None ):
            self.file.write ( "Amber Trajectory File\n" )
        else:
            self.file.write ( self.title[0:min ( 80, len ( self.title ) )] + "\n" )

    def WriteOwnerData ( self ):
        """Write data from the owner to a frame."""
        xyz = self.owner.coordinates3
        n   = 0
        for s in self.atomIndices:
            for j in range ( 3 ):
                self.file.write ( "{:8.3f}".format ( xyz[s,j] ) )
                n += 1
                if ( n >= 10 ):
                    n = 0
                    self.file.write ( "\n" )
        if n != 0: self.file.write ( "\n" )
        if self.includeBox:
            sp = self.owner.symmetryParameters
            self.file.write ( "{:8.3f}{:8.3f}{:8.3f}\n".format ( sp.a, sp.b, sp.c ) )
        self.numberOfFrames += 1
        self.numberOfWrites += 1

#===================================================================================================================================
# . Exporter definitions.
#===================================================================================================================================
_Exporter.AddHandler ( { TrajectoryMixin : AmberTrajectoryFileWriter.FromPathAndOwner } , [ "mdcrd", "MDCRD" ], "Amber Trajectory" )

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass
