"""Classes and functions for writing DCD trajectory files."""

from  pCore           import logFile, LogFileActive
from .ExportImport    import _Exporter
from .TrajectoryMixin import TrajectoryMixin

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class DCDTrajectoryFileWriter:
    """DCD trajectory file writer."""

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
        self.title        = None

    def AssignOwnerData ( self ):
        """Assign owner data to the trajectory."""
        cdef Coordinates3       data3
        cdef Selection          atomIndices
        cdef SymmetryParameters symmetryParameters
        # . Get objects.
        atomIndices        = self.owner.freeAtoms
        data3              = self.owner.coordinates3
        symmetryParameters = self.owner.symmetryParameters
        # . Assignment.
        if atomIndices        is not None: DCDStatus_Check ( DCDHandle_SetAtomIndices        ( self.cObject, atomIndices.cObject        ) )
        if data3              is not None: DCDStatus_Check ( DCDHandle_SetData3              ( self.cObject, data3.cObject              ) )
        if symmetryParameters is not None: DCDStatus_Check ( DCDHandle_SetSymmetryParameters ( self.cObject, symmetryParameters.cObject ) )

    def Close ( self ):
        """Close the file."""
        if self.isOpen:
            DCDWrite_Close ( &self.cObject )
            self.isOpen = False

    @classmethod
    def FromPathAndOwner ( selfClass, path, owner, title = "DCD Trajectory created by pDynamo" ):
        """Constructor given path, owner and other options."""
        self       = selfClass ( )
        self.path  = path
        self.owner = owner
        self.title = title
        self.AssignOwnerData ( )
        self.Open ( )
        return self

    def Open ( self ):
        """Open the file."""
        cdef char *path
        if not self.isOpen:
            byteString  = self.path.encode ( "UTF-8" )
            path        = byteString
            DCDStatus_Check ( DCDWrite_Open ( self.cObject, path ) )
            self.isOpen = True

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ) and ( self.cObject != NULL ) and self.isOpen:
            log.SummaryOfItems ( self.SummaryItems ( ), title = "DCD Trajectory File Writer" )

    def SummaryItems ( self ):
        """Summary items."""
        return [ ( "Atoms"        , "{:d}".format ( self.cObject.numberOfAtoms        ) ) ,
                 ( "Frames"       , "{:d}".format ( self.cObject.numberOfFrames       ) ) ,
                 ( "Atom Indices" , "{:d}".format ( self.cObject.numberOfAtomIndices  ) ) ,
                 ( "Has Symmetry" , "{!r}".format ( self.cObject.hasUnitCell == CTrue ) ) ]

    def WriteFooter ( self ):
        """Write the trajectory footer."""
        pass

    def WriteHeader ( self ):
        """Write the trajectory header."""
        cdef char *title
        byteString = self.title.encode ( "UTF-8" )
        title      = byteString
        DCDWrite_Header ( self.cObject, title )

    def WriteOwnerData ( self ):
        """Write data from the owner to a frame."""
        self.AssignOwnerData ( )
        DCDStatus_Check ( DCDWrite_Frame ( self.cObject ) )

    @property
    def currentFrame   ( self ): return DCDHandle_CurrentFrame   ( self.cObject )
    @property
    def numberOfFrames ( self ): return DCDHandle_NumberOfFrames ( self.cObject )
    @property
    def numberOfWrites ( self ): self.currentFrame

#===================================================================================================================================
# . Exporter definitions.
#===================================================================================================================================
_Exporter.AddHandler ( { TrajectoryMixin : DCDTrajectoryFileWriter.FromPathAndOwner } , [ "dcd", "DCD" ], "DCD Trajectory" )
