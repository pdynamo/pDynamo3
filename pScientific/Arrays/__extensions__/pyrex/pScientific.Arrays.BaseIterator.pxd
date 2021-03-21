from pCore.CPrimitiveTypes  cimport CBoolean           , \
                                    CFalse             , \
                                    CInteger           , \
                                    CReal              , \
                                    CTrue
from pCore.IntegerUtilities cimport Integer_Allocate   , \
                                    Integer_Deallocate
from pCore.Status           cimport CStatus            , \
                                    CStatus_OK

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "Iterator.h":

    # . The iterator type.
    ctypedef struct CIterator "Iterator":
        pass

    # . Functions.
    cdef CIterator *Iterator_Allocate     ( CStatus    *status ) 
    cdef CIterator *Iterator_Clone        ( CIterator  *self   , 
                                            CStatus    *status ) 
    cdef CInteger   Iterator_CurrentIndex ( CIterator  *self   )
    cdef CInteger   Iterator_DataOffSet   ( CIterator  *self   )
    cdef void       Iterator_Deallocate   ( CIterator **self   )
    cdef CInteger  *Iterator_Dump         ( CIterator  *self   ,
                                            CInteger   *n      ,
                                            CStatus    *status )
    cdef CInteger   Iterator_GetSize      ( CIterator  *self   )
    cdef CIterator *Iterator_Load         ( CInteger    n      ,
                                            CInteger   *state  ,
                                            CStatus    *status )
    cdef CInteger   Iterator_NextIndex    ( CIterator  *self   )
    cdef void       Iterator_Reset        ( CIterator  *self   )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class BaseIterator:

    cdef CIterator     *cIterator
    cdef public object  block
    cdef public object  isOwner

    cdef _SetDataPointer ( self )
