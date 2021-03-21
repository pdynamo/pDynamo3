from pCore.CPrimitiveTypes cimport CBoolean, CFalse, CInteger, CReal, CTrue
from pCore.Selection       cimport Selection
from pCore.Status          cimport CStatus, CStatus_OK

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "IntegerBlock.h":

    # . The array type.
    ctypedef struct CIntegerBlock "IntegerBlock":
        CInteger  capacity
        CInteger *items

    # . Functions.
    cdef CIntegerBlock *IntegerBlock_Allocate    ( CInteger        capacity ,
                                                   CStatus        *status   )
    cdef CIntegerBlock *IntegerBlock_Clone       ( CIntegerBlock  *self     ,
                                                   CStatus        *status   )
    cdef void           IntegerBlock_CopyTo      ( CIntegerBlock  *self     ,
                                                   CIntegerBlock  *other    ,
                                                   CStatus        *status   )
    cdef void           IntegerBlock_Deallocate  ( CIntegerBlock **self     )
    cdef void           IntegerBlock_Dereference ( CIntegerBlock  *self     )
    cdef CInteger       IntegerBlock_GetItem     ( CIntegerBlock  *self     ,
                                                   CInteger        index    ,
                                                   CStatus        *status   )
    cdef void           IntegerBlock_Reference   ( CIntegerBlock  *self     )
    cdef void           IntegerBlock_Set         ( CIntegerBlock  *self     ,
                                                   CInteger        value    )
    cdef void           IntegerBlock_SetItem     ( CIntegerBlock  *self     ,
                                                   CInteger        index    ,
                                                   CInteger        value    ,
                                                   CStatus        *status   )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class IntegerBlock:

    cdef CIntegerBlock *cObject
