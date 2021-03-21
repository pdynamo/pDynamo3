from pCore.CPrimitiveTypes          cimport CBoolean      , \
                                            CFalse        , \
                                            CInteger      , \
                                            CReal         , \
                                            CTrue            
from pCore.Status                   cimport CStatus       , \
                                            CStatus_OK       
from pMolecule.QCModel.BlockStorage cimport BlockStorage  , \
                                            CBlockStorage

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "BlockStorageContainer.h":

    ctypedef struct CBlockStorageContainer "BlockStorageContainer":
        CBoolean        isOwner 
        CInteger        capacity
        CBlockStorage **entries 

    cdef CBlockStorageContainer *BlockStorageContainer_Allocate   ( CInteger                 capacity ,
                                                                    CStatus                 *status   )
    cdef void                    BlockStorageContainer_Deallocate ( CBlockStorageContainer **self     )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class BlockStorageContainer:

    cdef CBlockStorageContainer *cObject
    cdef public object           isOwner
