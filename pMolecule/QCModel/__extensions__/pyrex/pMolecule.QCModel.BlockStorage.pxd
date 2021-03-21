from pCore.CPrimitiveTypes cimport CBoolean   , \
                                   CFalse     , \
                                   CInteger   , \
                                   CReal      , \
                                   CTrue
from pCore.Status          cimport CStatus    , \
                                   CStatus_OK

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "BlockStorage.h":

    ctypedef struct CBlockStorage "BlockStorage":
        pass

    cdef CBlockStorage *BlockStorage_Allocate   ( CStatus        *status )
    cdef CReal          BlockStorage_ByteSize   ( CBlockStorage  *self   )
    cdef CInteger       BlockStorage_Count      ( CBlockStorage  *self   )
    cdef void           BlockStorage_Deallocate ( CBlockStorage **self   )
    cdef void           BlockStorage_Empty      ( CBlockStorage  *self   )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class BlockStorage:

    cdef CBlockStorage *cObject
    cdef public object  isOwner
