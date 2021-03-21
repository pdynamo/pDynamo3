from pCore.CPrimitiveTypes cimport CBoolean, CFalse, CInteger, CReal, CTrue
from pCore.Status          cimport CStatus, CStatus_OK

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "BooleanBlock.h":

    # . The array type.
    ctypedef struct CBooleanBlock "BooleanBlock":
        CInteger   capacity
        CBoolean  *items

    # . Functions.
    cdef CBooleanBlock *BooleanBlock_Allocate    ( CInteger        capacity ,
                                                   CStatus        *status   )
    cdef CBooleanBlock *BooleanBlock_Clone       ( CBooleanBlock  *self     ,
                                                   CStatus        *status   )
    cdef void           BooleanBlock_CopyTo      ( CBooleanBlock  *self     ,
                                                   CBooleanBlock  *other    ,
                                                   CStatus        *status   )
    cdef void           BooleanBlock_Deallocate  ( CBooleanBlock **self     )
    cdef void           BooleanBlock_Dereference ( CBooleanBlock  *self     )
    cdef CBoolean       BooleanBlock_GetItem     ( CBooleanBlock  *self     ,
                                                   CInteger        index    ,
                                                   CStatus        *status   )
    cdef void           BooleanBlock_Negate      ( CBooleanBlock  *self     )
    cdef void           BooleanBlock_Reference   ( CBooleanBlock  *self     )
    cdef void           BooleanBlock_Set         ( CBooleanBlock  *self     ,
                                                   CBoolean        value    )
    cdef void           BooleanBlock_SetItem     ( CBooleanBlock  *self     ,
                                                   CInteger        index    ,
                                                   CBoolean        value    ,
                                                   CStatus        *status   )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class BooleanBlock:

    cdef CBooleanBlock *cObject
