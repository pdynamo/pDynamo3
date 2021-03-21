from pCore.CPrimitiveTypes           cimport CBoolean              , \
                                             CFalse                , \
                                             CInteger              , \
                                             CReal                 , \
                                             CTrue
from pCore.IntegerUtilities          cimport Integer_Allocate      , \
                                             Integer_Deallocate
from pCore.Selection                 cimport CSelection            , \
                                             Selection
from pCore.Status                    cimport CStatus               , \
                                             CStatus_OK
from pScientific.Arrays.BaseArray1D  cimport BaseArray1D           , \
                                             CView1D
from pScientific.Arrays.BaseIterator cimport BaseIterator          , \
                                             CIterator
from pScientific.Arrays.Slicing      cimport CMultiSlice           , \
                                             CProcessMultiSlice    , \
                                             MultiSlice_Allocate   , \
                                             MultiSlice_Deallocate

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "View2D.h":

    # . The N-D view type.
    ctypedef struct CView2D "View2D":
        pass

    # . Functions.
    cdef void       View2D_CopyTo           ( CView2D     *self       ,
                                              CView2D     *other      )
    cdef void       View2D_Deallocate       ( CView2D    **self       )
    cdef CInteger   View2D_GetColumns       ( CView2D     *self       )
    cdef CInteger   View2D_GetOffset        ( CView2D     *self       )
    cdef CInteger   View2D_GetRows          ( CView2D     *self       )
    cdef CInteger   View2D_GetSize          ( CView2D     *self       )
    cdef CInteger   View2D_GetStride        ( CView2D     *self       ,
                                              CInteger     dimension  ,
                                              CStatus     *status     )
    cdef void       View2D_Initialize       ( CView2D     *self       ,
                                              CInteger     rows       ,
                                              CInteger     columns    ,
                                              CStatus     *status     )
    cdef CIterator *View2D_MakeIterator     ( CView2D     *self       ,
                                              CStatus     *status     )
    cdef CIterator *View2D_MakeRowIterator  ( CView2D     *self       ,
                                              CSelection  *rows       ,
                                              CStatus     *status     )
    cdef void       View2D_SetState         ( CView2D     *self       ,
                                              CInteger     extent0    ,
                                              CInteger     extent1    ,
                                              CInteger     offset     ,
                                              CInteger     size       ,
                                              CInteger     stride0    ,
                                              CInteger     stride1    )
    cdef void       View2D_View1DMultiSlice ( CView2D     *self       ,
                                              CMultiSlice *multiSlice ,
                                              CView1D     *view       ,
                                              CStatus     *status     )
    cdef void       View2D_ViewMultiSlice   ( CView2D     *self       ,
                                              CMultiSlice *multiSlice ,
                                              CView2D     *view       ,
                                              CStatus     *status     )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class BaseArray2D:

    cdef CMultiSlice   *cMultiSlice
    cdef CView2D       *cView
    cdef public object  block
    cdef        object _iterator

    cdef void _AllocateCObject   ( self )
    cdef void _DeallocateCObject ( self )
