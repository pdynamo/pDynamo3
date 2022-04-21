from pCore.CPrimitiveTypes cimport CBoolean    , \
                                   CCardinal16 , \
                                   CCardinal32 , \
                                   CFalse      , \
                                   CInteger    , \
                                   CReal       , \
                                   CTrue
from pCore.List            cimport CList       , \
                                   List_Iterate_Initialize
from pCore.Status          cimport CStatus     , \
                                   CStatus_OK

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "BlockStorage.h":

    ctypedef struct CBlock "Block":
        CInteger     count
        CCardinal16 *indices16
        CCardinal32 *indices32
        CReal       *data

    ctypedef struct CBlockStorage "BlockStorage":
        CBoolean  checkUnderFlow
        CInteger  blockSize
        CInteger  count
        CInteger  nIndices16
        CInteger  nIndices32
        CInteger  nReal
        CReal     underFlow
        CList    *blocks

    cdef CBlockStorage *BlockStorage_Allocate   ( CStatus        *status )
    cdef CReal          BlockStorage_ByteSize   ( CBlockStorage  *self   )
    cdef CInteger       BlockStorage_Count      ( CBlockStorage  *self   )
    cdef void           BlockStorage_Deallocate ( CBlockStorage **self   )
    cdef void           BlockStorage_Empty      ( CBlockStorage  *self   )
    cdef CBlock        *BlockStorage_Iterate    ( CBlockStorage  *self   )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class BlockStorage:

    cdef CBlockStorage *cObject
    cdef public object  isOwner


