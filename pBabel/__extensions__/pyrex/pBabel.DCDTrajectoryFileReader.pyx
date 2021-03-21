"""Classes and functions for reading DCD trajectory files."""

from  pCore           import logFile, LogFileActive
from .ExportImport    import _Importer
from .TrajectoryMixin import TrajectoryMixin

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class DCDTrajectoryFileReader:
    """DCD trajectory file reader."""

    def __dealloc__ ( self ):
        """Finalization."""
        self.Close ( )

    def __init__ ( self ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate   ( )

    def __len__ ( self ): return self.numberOfFrames

    def _Allocate ( self ):
        """Allocation."""
        self.cObject = DCDHandle_Allocate ( )
        if self.cObject == NULL: DCDStatus_Check ( CDCDStatus_OutOfMemory )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject      = NULL
        self.isOpen       = False
        self.isTrajectory = True
        self.owner        = None
        self.path         = None

    def AssignOwnerData ( self ):
        """Assign owner data to the trajectory."""
        cdef Coordinates3       data3
        cdef SymmetryParameters symmetryParameters
        # . Get objects.
        data3              = self.owner.coordinates3
        symmetryParameters = self.owner.symmetryParameters
        # . Assignment.
        if data3              is not None: DCDStatus_Check ( DCDHandle_SetData3              ( self.cObject, data3.cObject              ) )
        if symmetryParameters is not None: DCDStatus_Check ( DCDHandle_SetSymmetryParameters ( self.cObject, symmetryParameters.cObject ) )

    def Close ( self ):
        """Close the file."""
        if self.isOpen:
            DCDRead_Close ( &self.cObject )
            self.isOpen = False

    @classmethod
    def FromPathAndOwner ( selfClass, path, owner ):
        """Constructor given path, owner and other options."""
        self       = selfClass ( )
        self.path  = path
        self.owner = owner
        self.AssignOwnerData ( )
        self.Open ( )
        return self

    def Open ( self ):
        """Open the file."""
        cdef char *path
        if not self.isOpen:
            byteString  = self.path.encode ( "UTF-8" )
            path        = byteString
            DCDStatus_Check ( DCDRead_Open ( self.cObject, path ) )
            self.isOpen = True

# . Useful if owner data sources change?
#    def ResetOwnerData ( self ):
#        """Reset owner data."""
#        pass

    def ReadFooter ( self ):
        """Read the trajectory footer."""
        pass

    def ReadHeader ( self ):
        """Read the trajectory header."""
        DCDStatus_Check ( DCDRead_Header               ( self.cObject ) )
        DCDStatus_Check ( DCDHandle_CheckNumberOfAtoms ( self.cObject, len ( self.owner.atoms ) ) )
        DCDStatus_Check ( DCDHandle_AllocateQW         ( self.cObject ) )

    def RestoreOwnerData ( self ):
        """Restore data from a frame to the owner."""
        try:
            if self.currentFrame >= self.numberOfFrames: raise EOFError
            DCDStatus_Check ( DCDRead_Frame ( self.cObject ) )
            return True
        except EOFError:
            DCDStatus_Check ( DCDRead_GotoFrame ( self.cObject, 0 ) )
            return False

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ) and ( self.cObject != NULL ) and self.isOpen:
            log.SummaryOfItems ( self.SummaryItems ( ), title = "DCD Trajectory File Reader" )

    def SummaryItems ( self ):
        """Summary items."""
        return [ ( "Atoms"        , "{:d}".format ( self.cObject.numberOfAtoms        ) ) ,
                 ( "Frames"       , "{:d}".format ( self.cObject.numberOfFrames       ) ) ,
                 ( "Atom Indices" , "{:d}".format ( self.cObject.numberOfAtomIndices  ) ) ,
                 ( "Has Symmetry" , "{!r}".format ( self.cObject.hasUnitCell == CTrue ) ) ]

    @property
    def currentFrame   ( self ): return DCDHandle_CurrentFrame   ( self.cObject )
    @property
    def numberOfFrames ( self ): return DCDHandle_NumberOfFrames ( self.cObject )
    @property
    def numberOfReads  ( self ): self.currentFrame

#===================================================================================================================================
# . Importer definitions.
#===================================================================================================================================
_Importer.AddHandler ( { TrajectoryMixin : DCDTrajectoryFileReader.FromPathAndOwner } , [ "dcd", "DCD" ], "DCD Trajectory" )
