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

    # . For debugging only!
    def __getstate__ ( self ):
        """Return the state."""
        cdef CBlock   *cBlock ;
        cdef CInteger  c, i, i16, i32, iR, n16, n32, nR
        # . Basic data.
        n16 = self.cObject.nIndices16
        n32 = self.cObject.nIndices32
        nR  = self.cObject.nReal
        state = { "Block Size" : self.cObject.blockSize ,
                  "Count"      : self.cObject.count     ,
                  "Indices 16" : n16                    ,
                  "Indices 32" : n32                    ,
                  "Real"       : nR                     }
        if self.cObject.checkUnderFlow == CTrue:
            state["Underflow"] = self.cObject.underFlow
        # . Blocks.
        records = []
        List_Iterate_Initialize ( self.cObject.blocks )
        while True:
            cBlock = BlockStorage_Iterate ( self.cObject )
            if cBlock == NULL: break
            i16 = 0
            i32 = 0
            iR  = 0
            for c from 0 <= c < cBlock.count:
                record = []
                for i from 0 <= i < n16:
                    record.append ( cBlock.indices16[i16+i] )
                for i from 0 <= i < n32:
                    record.append ( cBlock.indices32[i32+i] )
                for i from 0 <= i < nR:
                    record.append ( cBlock.data[iR+i] )
                i16 += n16
                i32 += n32
                iR  += nR
                records.append ( record )
        # . Finish up.
        state["Records"] = records
        return state

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
