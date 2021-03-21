from pCore.CPrimitiveTypes cimport CBoolean, CFalse, CInteger, CReal, CTrue
from pCore.Selection       cimport Selection
from pCore.Status          cimport CStatus, CStatus_OK

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "RealBlock.h":

    # . The array type.
    ctypedef struct CRealBlock "RealBlock":
        CInteger  capacity
        CReal    *items

    # . Functions.
    cdef CRealBlock *RealBlock_Allocate    ( CInteger     capacity ,
                                             CStatus     *status   )
    cdef CRealBlock *RealBlock_Clone       ( CRealBlock  *self     ,
                                             CStatus     *status   )
    cdef void        RealBlock_CopyTo      ( CRealBlock  *self     ,
                                             CRealBlock  *other    ,
                                             CStatus     *status   )
    cdef void        RealBlock_Deallocate  ( CRealBlock **self     )
    cdef void        RealBlock_Dereference ( CRealBlock  *self     )
    cdef CReal       RealBlock_GetItem     ( CRealBlock  *self     ,
                                             CInteger     index    ,
                                             CStatus     *status   )
    cdef void        RealBlock_Reference   ( CRealBlock  *self     )
    cdef void        RealBlock_Set         ( CRealBlock  *self     ,
                                             CReal        value    )
    cdef void        RealBlock_SetItem     ( CRealBlock  *self     ,
                                             CInteger     index    ,
                                             CReal        value    ,
                                             CStatus     *status   )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class RealBlock:

    cdef CRealBlock *cObject
