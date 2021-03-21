"""Classes and functions for manipulating system geometry trajectories."""

import glob, os, os.path

from  pCore                 import AttributableObject  , \
                                   Clone               , \
                                   Pickle              , \
                                   PickleFileExtension , \
                                   Unpickle
from  pScientific.Geometry3 import Coordinates3
from  pScientific.Symmetry  import SymmetryParameters
from .ExportImport          import _Exporter           , \
                                   _Importer
from .TrajectoryMixin       import TrajectoryMixin

# . Really need an index array in the header and some way of renumbering the frames after random access if frames added or removed.

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Frame naming.
_FramePrefix  = "frame"
_FramePostfix = PickleFileExtension

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class SystemGeometryTrajectory ( AttributableObject ):
    """Class for system geometry trajectories."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "isReadable"            : False ,
                             "isTrajectory"          :  True ,
                             "isWritable"            : False ,
                             "hasSymmetryParameters" : False ,
                             "numberOfFrames"        :     0 ,
                             "numberOfReads"         :     0 ,
                             "numberOfWrites"        :     0 ,
                             "owner"                 :  None ,
                             "path"                  :  None ,
                             "position"              :    -1 } )

    def __getitem__ ( self, index ):
        """Get an item."""
        if ( index < 0 ) or ( index >= len ( self ) ): raise IndexError
        else:                                          return self.ReadFrame ( frame = index )

    def __len__ ( self ):
        return self.numberOfFrames

    def Close ( self ):
        """Close the trajectory."""
        pass

    def Open ( self ):
        """Open the trajectory."""
        # . Check to see if the path exists.
        pathExists = os.access ( self.path, os.F_OK )
        if pathExists:
            if not os.path.isdir ( self.path ): raise IOError ( "Trajectory exists that is not a directory." )
        else:
            if not self.isWritable:
                if self.isReadable: raise IOError ( "Read-only mode specified for a trajectory that does not exist." )
                else:               raise IOError ( "Neither read nor write modes specified for a trajectory." )
            os.mkdir ( self.path )
        # . Check for readability and writeability.
        isReadable = os.access ( self.path, os.R_OK )
        isWritable = os.access ( self.path, os.W_OK )
        # . Check for consistency with the modes.
        if self.isReadable and ( not isReadable ): raise IOError ( "Read mode specified for a trajectory that is unreadable."   )
        if self.isWritable and ( not isWritable ): raise IOError ( "Write mode specified for a trajectory that is unwriteable." )
        # . Check the contents of an existing trajectory.
        if pathExists:
            frames = glob.glob ( os.path.join ( self.path, _FramePrefix + "*" + _FramePostfix ) )
            # . Set the number of frames for a readable trajectory.
            if self.isReadable:
                self.numberOfFrames = len ( frames )
            # . Clear the files if the trajectory is only writable.
            else:
                for frame in frames: os.remove ( frame )
        # . Check to see if the owner has symmetry parameters.
        self.hasSymmetryParameters = hasattr ( self.owner, "symmetry" ) and ( self.owner.symmetry is not None )

    @classmethod
    def ReaderFromPathAndOwner ( selfClass, path, owner ):
        """Constructor given path, owner and other options."""
        self = selfClass.WithOptions ( isReadable = True  ,
                                       owner      = owner ,
                                       path       = path  )
        self.Open ( )
        return self

    def ReadFooter ( self ):
        """Read the trajectory footer."""
        return {}

    def ReadFrame ( self, frame = None, identifier = None, instance = None ):
        """Read a frame."""
        if self.isReadable:
            data = Unpickle ( os.path.join ( self.path, _FramePrefix + "{:d}".format ( frame ) + _FramePostfix ) )
            self.numberOfReads += 1
            if ( instance is None ) or isinstance ( data, instance ):
                return data
            else:
                for item in data:
                    if isinstance ( item, instance ): return item
        else: raise IOError ( "Reading from trajectory that is not readable." )

    def ReadHeader ( self ):
        """Read the trajectory header."""
        return {}

    def RestoreOwnerData ( self, index = -1 ):
        """Restore data from a frame to the owner."""
        # . Initialization.
        data = None
        # . Restore data from a specific frame.
        if index >= 0:
            data = self.ReadFrame ( frame = index )
            self.position = index + 1
        # . Restore data from the next frame in the sequence.
        else:
            self.position += 1
            if self.position >= self.numberOfFrames:
                self.position = -1
            else:
                data = self.ReadFrame ( frame = self.position )
        # . Restore any data.
        if data is None:
            return False
        else:
            if isinstance ( data, Coordinates3 ):
                self.owner.coordinates3       = data
            elif isinstance ( data, tuple ):
                self.owner.coordinates3       = data[0]
                self.owner.symmetryParameters = data[1]
            return True

    def WriteFooter ( self ):
        """Write the trajectory footer."""
        pass

    def WriteFrame ( self, data, frame = None, identifier = None ):
        """Write a frame."""
        if self.isWritable:
            if frame is None: f = self.numberOfFrames
            else:             f = frame
            Pickle ( os.path.join ( self.path, _FramePrefix + "{:d}".format ( f ) + _FramePostfix ), data )
            if frame is None: self.numberOfFrames += 1
            self.numberOfWrites += 1
        else: raise IOError ( "Writing to trajectory that is not writeable." )

    def WriteHeader ( self ):
        """Write the trajectory header."""
        pass

    def WriteOwnerData ( self, index = -1 ):
        """Write data from the owner to a frame."""
        if index >= 0 : frame = index
        else:           frame = None
        if self.hasSymmetryParameters: data = ( self.owner.coordinates3, self.owner.symmetryParameters )
        else:                          data =   self.owner.coordinates3
        self.WriteFrame ( data, frame = frame )

    @classmethod
    def WriterFromPathAndOwner ( selfClass, path, owner, append = False ):
        """Constructor given path, owner and other options."""
        self = selfClass.WithOptions ( isReadable = append ,
                                       isWritable = True   ,
                                       owner      = owner  ,
                                       path       = path   )
        self.Open ( )
        return self

#===================================================================================================================================
# . Exporter and importer definitions.
#===================================================================================================================================
_Exporter.AddHandler ( { TrajectoryMixin : SystemGeometryTrajectory.WriterFromPathAndOwner } , [ "ptGeo" ], "System Geometry Trajectory" )
_Importer.AddHandler ( { TrajectoryMixin : SystemGeometryTrajectory.ReaderFromPathAndOwner } , [ "ptGeo" ], "System Geometry Trajectory" )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
