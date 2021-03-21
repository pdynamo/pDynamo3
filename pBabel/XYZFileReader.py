"""Classes for reading XYZ files and trajectories."""

from  pCore                 import logFile         , \
                                   LogFileActive   , \
                                   TextFileReader
from  pMolecule             import System
from  pScientific           import PeriodicTable
from  pScientific.Geometry3 import Coordinates3
from .ExportImport          import _Importer
from .TrajectoryMixin       import TrajectoryMixin

#===================================================================================================================================
# . XYZ file reader class.
#===================================================================================================================================
class XYZFileReader ( TextFileReader ):
    """XYZFileReader is the class for XYZ files that are to be read."""

    _classLabel = "XYZ File Reader"

    def Parse ( self, log = logFile ):
        """Parse the data on the file."""
        if not self.isParsed:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
            # . Open the file.
            self.Open ( )
            try:
                # . Number of atoms.
                items         = self.GetTokens ( converters = [ int ] )
                numberOfAtoms = items[0]
                # . Title line.
                self.title    = self.GetLine ( )
                # . XYZ lines.
                self.atomicNumbers = []
                self.coordinates3  = Coordinates3.WithExtent ( numberOfAtoms )
                for i in range ( numberOfAtoms ):
                    items = self.GetTokens ( converters = [ PeriodicTable.AtomicNumber, float, float, float ] )
                    self.atomicNumbers.append ( items[0] )
                    self.coordinates3[i,0] = items[1]
                    self.coordinates3[i,1] = items[2]
                    self.coordinates3[i,2] = items[3]
            except EOFError:
                pass
            # . Close the file.
            self.WarningStop ( )
            self.Close       ( )
            # . Set the parsed flag and some other options.
            self.isParsed = True
            self.log      = None

    @classmethod
    def PathToCoordinates3 ( selfClass, path ):
        """Return the coordinates from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse ( )
        return inFile.ToCoordinates3 ( )

    @classmethod
    def PathToSystem ( selfClass, path ):
        """Return the system from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse ( )
        return inFile.ToSystem ( )

    def ToCoordinates3 ( self ):
        """Return a coordinates3 object."""
        if self.isParsed and hasattr ( self, "coordinates3" ): return self.coordinates3
        else:                                                  return None

    def ToSystem ( self ):
        """Return a system."""
        system = None
        if self.isParsed and hasattr ( self, "atomicNumbers" ):
            system              = System.FromAtoms ( self.atomicNumbers )
            system.label        = self.title
            system.coordinates3 = self.ToCoordinates3 ( )
        return system

#===================================================================================================================================
# . XYZTrajectory file reader class.
#===================================================================================================================================
class XYZTrajectoryFileReader ( TextFileReader ):
    """Class for reading XYZ trajectory files."""

    _attributable = dict ( TextFileReader._attributable )
    _classLabel   = "XYZ Trajectory File Reader"
    _attributable.update ( { "isComplete"     : True ,
                             "isTrajectory"   : True ,
                             "numberOfReads"  :  0   ,
                             "numberOfFrames" : -1   , # . Number of frames unknown until reading finished.
                             "owner"          : None } )

    def __len__ ( self ): return self.numberOfFrames

    def _NextFrame ( self ):
        """Get the next frame on the file."""
        # . Number of atoms.
        items           = self.GetTokens ( converters = [ int ] )
        numberOfAtoms   = items[0]
        self.isComplete = False
        # . Title line.
        self.title = self.GetLine ( )
        # . XYZ lines.
        atoms = self.owner.atoms
        xyz   = self.owner.coordinates3
        isOK  = ( numberOfAtoms == len ( atoms ) )
        if isOK:
            for i in range ( numberOfAtoms ):
                items = self.GetTokens ( converters = [ PeriodicTable.AtomicNumber, float, float, float ] )
                if items[0] != atoms[i].atomicNumber:
                    isOK = False
                    break
                xyz[i,0] = items[1]
                xyz[i,1] = items[2]
                xyz[i,2] = items[3]
            self.isComplete = True
        if not isOK: raise IOError ( "Frame data does not match owner definition." )

    @classmethod
    def FromPathAndOwner ( selfClass, path, owner, **options ):
        """Constructor given path, owner and other options."""
        options          = dict ( options )
        options["path" ] = path
        options["owner"] = owner
        self = selfClass.WithOptions ( **options )
        self.Open ( )
        return self

    def ReadFooter ( self ):
        """Read the trajectory footer."""
        pass

    def ReadHeader ( self ):
        """Read the trajectory header."""
        # . Activate the trajectory if necessary.
        if not self.isActive:
            self.numberOfReads = 0
            self.Open ( )

    def RestoreOwnerData ( self, index = -1 ):
        """Restore data from a frame to the owner."""
        try:
            self._NextFrame ( )
            self.numberOfReads += 1
            return True
        except EOFError:
            if self.isComplete:
                self.numberOfFrames = self.numberOfReads
                return False
            else:
                raise IOError ( "Incomplete frame." )

#===================================================================================================================================
# . Importer definitions.
#===================================================================================================================================
_Importer.AddHandler ( { Coordinates3    : XYZFileReader.PathToCoordinates3         ,
                         System          : XYZFileReader.PathToSystem               ,
                         TrajectoryMixin : XYZTrajectoryFileReader.FromPathAndOwner } ,
                       [ "xyz" ], "XYZ", defaultFunction = XYZFileReader.PathToSystem )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
