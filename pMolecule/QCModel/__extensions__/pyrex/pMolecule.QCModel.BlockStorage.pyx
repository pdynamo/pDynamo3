"""Handle block storage."""

from  pScientific import Magnitude_Adjust

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class BlockStorage:

    # . Public methods.
    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner: BlockStorage_Deallocate ( &self.cObject )

    def __init__ ( self ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate   ( )

    def __len__ ( self ): return self.count

    def _Allocate ( self ):
        """Allocation."""
        self.cObject = BlockStorage_Allocate ( NULL )
        self.isOwner = True

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    def Empty ( self ):
        """Empty the storage."""
        BlockStorage_Empty ( self.cObject )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    @property
    def byteSize ( self ):
        return Magnitude_Adjust ( BlockStorage_ByteSize ( self.cObject ) )

    @property
    def count ( self ):
        return BlockStorage_Count ( self.cObject )
