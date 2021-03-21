"""Classes and functions for manipulating system internal coordinate trajectories."""

# . Currently:
#
#   - Only distances are supported, although it would be straightforward to add other ICs.
#
#   - All ICs need to be stated specifically. It would be nice to have also classes of IC too
#     (e.g. all water oxygens within X Angstroms of a given atom) without specific indexation.
#     Would need pairlist and updating mechanism.
#

import glob, os, os.path, random

from  pCore                  import AttributableObject  , \
                                    Pickle              , \
                                    PickleFileExtension , \
                                    Unpickle
from .ExportImport           import _Exporter           , \
                                    _Importer
from .TrajectoryMixin        import TrajectoryMixin

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Default block size.
_DefaultBlockSize = 1024

# . Block naming.
_BlockPrefix  = "block"
_BlockPostfix = PickleFileExtension

# . Header and footer naming.
_FooterName = "footer"
_HeaderName = "header"

#===================================================================================================================================
# . Reader class.
#===================================================================================================================================
class SystemICTrajectoryReader ( AttributableObject ):
    """Class for reading system IC trajectories."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "blocks"          :     0 ,
                             "blockSize"       :     0 ,
                             "currentBlock"    :    -1 ,
                             "current"         :     0 ,
                             "data"            :  None ,
                             "distanceCutOff"  :   0.0 ,
                             "distanceIndices" :  None ,
                             "isTrajectory"    :  True ,
                             "numberOfFrames"  :    -1 ,
                             "numberOfReads"   :     0 ,
                             "owner"           :  None ,
                             "path"            :  None } )

    def Close ( self ):
        """Close the trajectory."""
        pass

    @classmethod
    def FromPathAndOwner ( selfClass, path, owner ):
        """Constructor given path, owner and other options."""
        self = selfClass.WithOptions ( owner = owner, path = path )
        self.Open ( )
        return self

    def Open ( self ):
        """Open the trajectory."""
        # . Check that the trajectory exists, is a directory and is readable.
        if os.access ( self.path, os.F_OK ) and os.access ( self.path, os.R_OK ) and os.path.isdir ( self.path ):
            # . Check for a valid header.
            if not os.path.exists ( os.path.join ( self.path, _HeaderName + _BlockPostfix ) ): raise IOError ( "Unable to find trajectory header." )
            # . Find the number of blocks.
            self.blocks = len ( glob.glob ( os.path.join ( self.path, _BlockPrefix + "*" + _BlockPostfix ) ) )
        # . Invalid trajectory.
        else: raise IOError ( "Invalid or non-existent trajectory." )

    def ReadBlock ( self ):
        """Read a block of data."""
        if self.currentBlock < self.blocks:
            self.data = Unpickle ( os.path.join ( self.path, _BlockPrefix + "{:d}".format ( self.currentBlock ) + _BlockPostfix ) )
            self.blockSize     = len ( self.data )
            self.current       = 0
            self.currentBlock += 1
        else: raise IndexError ( "Invalid block index." )

    def ReadFooter ( self ):
        """Read the footer."""
        pass

    def ReadHeader ( self ):
        """Read the trajectory header."""
        header = Unpickle ( os.path.join ( self.path, _HeaderName + _BlockPostfix ) )
        for ( key, value ) in header.items ( ): setattr ( self, key, value )
        return header

    def RestoreOwnerData ( self ):
        """Restore data from a frame to the owner."""
        if self.current >= self.blockSize:
            if self.currentBlock >= self.blocks:
                self.numberOfFrames = self.numberOfReads
                return False
            else:
                self.ReadBlock ( )
        self.owner.scratch.icTerms = self.data[current] # . Put data into scratch.
        self.current       += 1
        self.numberOfReads += 1
        return True

    def ReturnAllFrameData ( self ):
        """Return all frame data as a list."""
        data = []
        self.currentBlock = 0
        for i in range ( self.blocks ):
            self.ReadBlock ( )
            data.extend ( self.data )
        return data

#===================================================================================================================================
# . Writer class.
#===================================================================================================================================
class SystemICTrajectoryWriter ( AttributableObject ):
    """Class for writing system IC trajectories."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "blocks"          :                 0 ,
                             "blockSize"       : _DefaultBlockSize ,
                             "current"         :                 0 ,
                             "data"            :              None ,
                             "distanceCutOff"  :               0.0 ,
                             "distanceIndices" :              None ,
                             "isAppendable"    :             False ,
                             "isTrajectory"    :              True ,
                             "numberOfFrames"  :                 0 ,
                             "numberOfWrites"  :                 0 ,
                             "owner"           :              None ,
                             "path"            :              None } )

    def _CheckOptions ( self ):
        """Check options."""
        super ( SystemICTrajectoryWriter, self )._CheckOptions ( )
        isOK = ( self.distanceCutOff  is not None ) and (       self.distanceCutOff    > 0.0 ) and \
               ( self.distanceIndices is not None ) and ( len ( self.distanceIndices ) > 0   )
        if not isOK: raise Exception ( "Invalid options to IC trajectory." )

    def Close ( self ):
        """Close the trajectory."""
        self.WriteBlock  ( )
        self.WriteFooter ( )

    @classmethod
    def FromPathAndOwner ( selfClass, path, owner, append = False, distanceCutOff = None, distanceIndices = None ):
        """Constructor given path, owner and other options."""
        self = selfClass.WithOptions ( distanceCutOff  = distanceCutOff  ,
                                       distanceIndices = distanceIndices ,
                                       isAppendable    = append          ,
                                       owner           = owner           ,
                                       path            = path            )
        self.Open ( )
        return self

    def Open ( self ):
        """Open the trajectory."""
        pathExists = os.access ( self.path, os.F_OK )
        if pathExists:
            if not os.path.isdir ( self.path ): raise IOError ( "Trajectory exists that is not a directory." )
        else:
            os.mkdir ( self.path )
        if not os.access ( self.path, os.W_OK ): raise IOError ( "Trajectory is not writeable." )
        if pathExists:
            if self.isAppendable:
                self.blocks = len ( glob.glob ( os.path.join ( self.path, _BlockPrefix + "*" + _BlockPostfix ) ) )
            else:
                for target in glob.glob ( os.path.join ( self.path, "*" ) ): os.remove ( target )

    def WriteBlock ( self ):
        """Write a block of data."""
        if self.current > 0:
            Pickle ( os.path.join ( self.path, _BlockPrefix + "{:d}".format ( self.blocks ) + _BlockPostfix ), self.data )
            self.blocks += 1
            self.current = 0
            self.data    = []

    def WriteFooter ( self ):
        """Write a footer."""
        pass

    def WriteHeader ( self ):
        """Write the trajectory header."""
        header = { "distanceCutoff"  : self.distanceCutOff  ,
                   "distanceIndices" : self.distanceIndices }
        Pickle ( os.path.join ( self.path, _HeaderName + _BlockPostfix ), header )
        self.data = []

    def WriteOwnerData ( self ):
        """Write data from the owner to a frame."""
        if self.current >= self.blockSize: self.WriteBlock ( )
        crd3 = self.owner.coordinates3
        ics  = []
        for ( i, j ) in self.distanceIndices:
            d = crd3.Distance ( i, j )
            if d <= self.distanceCutOff:
                ics.append ( ( i, j, d ) )
        self.data.append ( ics )
        self.current        += 1
        self.numberOfFrames += 1
        self.numberOfWrites += 1

#===================================================================================================================================
# . Exporter and importer definitions.
#===================================================================================================================================
_Exporter.AddHandler ( { TrajectoryMixin : SystemICTrajectoryWriter.FromPathAndOwner } , [ "ptIC" ], "System Internal Coordinate Trajectory" )
_Importer.AddHandler ( { TrajectoryMixin : SystemICTrajectoryReader.FromPathAndOwner } , [ "ptIC" ], "System Internal Coordinate Trajectory" )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
