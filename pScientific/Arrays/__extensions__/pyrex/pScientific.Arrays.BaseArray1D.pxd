from pCore.CPrimitiveTypes           cimport CBoolean           , \
                                             CFalse             , \
                                             CInteger           , \
                                             CReal              , \
                                             CTrue
from pCore.IntegerUtilities          cimport Integer_Allocate   , \
                                             Integer_Deallocate
from pCore.Selection                 cimport Selection
from pCore.Status                    cimport CStatus            , \
                                             CStatus_OK
from pScientific.Arrays.BaseIterator cimport BaseIterator       , \
                                             CIterator

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "View1D.h":

    # . The N-D view type.
    ctypedef struct CView1D "View1D":
        pass

    # . Functions.
    cdef void       View1D_CopyTo       ( CView1D     *self   ,
                                          CView1D     *other  )
    cdef void       View1D_Deallocate   ( CView1D    **self   )
    cdef CInteger   View1D_GetExtent    ( CView1D     *self   )
    cdef CInteger   View1D_GetOffset    ( CView1D     *self   )
    cdef CInteger   View1D_GetStride    ( CView1D     *self   )
    cdef void       View1D_Initialize   ( CView1D     *self   ,
                                          CInteger     extent ,
                                          CStatus     *status )
    cdef CIterator *View1D_MakeIterator ( CView1D     *self   ,
                                          CStatus     *status )
    cdef void       View1D_SetState     ( CView1D     *self   ,
                                          CInteger     extent ,
                                          CInteger     offset ,
                                          CInteger     size   ,
                                          CInteger     stride )
    cdef void       View1D_View         ( CView1D     *self   ,
                                          CInteger     start  ,
                                          CInteger     extent ,
                                          CInteger     stride ,
                                          CView1D     *view   ,
                                          CStatus     *status )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class BaseArray1D:

    cdef CView1D       *cView
    cdef public object  block
    cdef        object _iterator

    cdef void _AllocateCObject   ( self )
    cdef void _DeallocateCObject ( self )
