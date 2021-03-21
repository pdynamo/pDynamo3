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
cdef extern from "Slice.h":

    ctypedef struct CMultiSlice "MultiSlice":
        pass

    cdef CMultiSlice *MultiSlice_Allocate        ( CInteger      capacity  ,
                                                   CStatus      *status    )
    cdef void         MultiSlice_Deallocate      ( CMultiSlice **self      )
    cdef CInteger     MultiSlice_GetCapacity     ( CMultiSlice  *self      )
    cdef CInteger     MultiSlice_GetExtent       ( CMultiSlice  *self      ,
                                                   CInteger      dimension ,
                                                   CStatus      *status    )
    cdef CInteger     MultiSlice_GetRank         ( CMultiSlice  *self      )
    cdef CInteger     MultiSlice_GetSize         ( CMultiSlice  *self      )
    cdef void         MultiSlice_SetRank         ( CMultiSlice  *self      )
    cdef void         MultiSlice_SetSliceIndices ( CMultiSlice  *self      ,
                                                   CInteger      dimension ,
                                                   CBoolean      isScalar  ,
                                                   CInteger      start     ,
                                                   CInteger      stop      ,
                                                   CInteger      stride    ,
                                                   CInteger      n         )

    cdef CInteger     SliceIndices_CheckScalar   ( CInteger      index     ,
                                                   CInteger      extent    ,
                                                   CStatus      *status    )
    cdef CInteger     SliceIndices_CheckSlice    ( CInteger     *pStart    ,
                                                   CInteger     *pStop     ,
                                                   CInteger     *pStride   ,
                                                   CInteger      extent    ,
                                                   CInteger     *qStart    ,
                                                   CInteger     *qStop     ,
                                                   CInteger     *qStride   ,
                                                   CStatus      *status    )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class MultiSlice:

    cdef CMultiSlice   *cObject
    cdef public object  isOwner

cdef CInteger CProcessMultiSlice ( shape, indices, CMultiSlice *multiSlice )
