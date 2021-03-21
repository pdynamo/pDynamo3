"""Trajectory object mixin or template."""

# . TrajectoryMixin is currently not used directly but shows the structure that all trajectories (should) follow.
# . Note that trajectories that read and write need only implement read and write attributes and methods.

# . Some of these attributes and methods are taken care of by other classes (e.g. TextFile).
# . Some may also be unnecessary but must nevertheless be defined for general use.

# . Trajectories may be read only, write only or both read and write.
# . The frame index for read and write operations.

# . It would be nice to use this directly but the other classes (including Cython ones) need to be sorted out.

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class TrajectoryMixin:
    """Mixin for trajectory classes."""

    # . Attributes and/or properties - they can be implemented either way.
    _attributable = { "isTrajectory" : True , # . A fudge until the mixin is used directly.
                      "owner"        : None , # . The owner of the trajectory (where data is read to or written from).
                      "path"         : None } # . The trajectory path.

    @property
    def numberOfFrames ( self ):
        """The number of frames on the trajectory."""
        return -1 # . If not known (e.g. when reading and no header information) return -1.

    @property
    def numberOfReads ( self ):
        """The number of frames read."""
        return 0

    @property
    def numberOfWrites ( self ):
        """The number of frames written."""
        return 0

    # . General methods.
    def __len__ ( self ): return self.numberOfFrames

    @classmethod
    def FromPathAndOwner ( selfClass, path, owner, **options ):
        """Constructor given path and owner and other options."""
        return None

    def Close ( self ):
        """Close the trajectory."""
        pass

    def Open ( self ):
        """Open the trajectory."""
        pass

    # . Reading.
    def ReadFooter ( self ):
        """Read a trajectory footer."""
        # . Return dictionary of data.
        return {}

    def ReadHeader ( self ):
        """Read a trajectory header."""
        # . Return dictionary of data.
        return {} 

    def RestoreOwnerData ( self, index = -1 ): # . Index is the frame index.
        """Restore owner data."""
        # . Return False if no data left on the trajectory, true otherwise.
        return False 

    # . Writing.
    def WriteFooter ( self ):
        """Write a trajectory footer."""
        pass

    def WriteHeader ( self ):
        """Write a trajectory header."""
        pass

    def WriteOwnerData ( self, index = -1 ):
        """Write owner data."""
        pass

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
