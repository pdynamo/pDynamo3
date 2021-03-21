from pCore.CPrimitiveTypes           cimport CBoolean              , \
                                             CFalse                , \
                                             CInteger              , \
                                             CReal                 , \
                                             CTrue
from pCore.IntegerUtilities          cimport Integer_Allocate      , \
                                             Integer_Deallocate
from pCore.Status                    cimport CStatus               , \
                                             CStatus_OK
from pScientific.Arrays.BaseArray1D  cimport BaseArray1D           , \
                                             CView1D
from pScientific.Arrays.BaseArray2D  cimport BaseArray2D           , \
                                             CView2D
from pScientific.Arrays.BaseIterator cimport BaseIterator          , \
                                             CIterator
from pScientific.Arrays.Slicing      cimport CMultiSlice           , \
                                             CProcessMultiSlice    , \
                                             MultiSlice_Allocate   , \
                                             MultiSlice_Deallocate

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "ViewND.h":

    # . The N-D view type.
    ctypedef struct CViewND "ViewND":
        pass

    # . Functions.
    cdef CViewND   *ViewND_AllocateWithShape  ( CInteger     rank       ,
                                                CInteger    *shape      ,
                                                CStatus     *status     )
    cdef CViewND   *ViewND_Clone              ( CViewND     *self       ,
                                                CStatus     *status     )
    cdef void       ViewND_Deallocate         ( CViewND    **self       )
    cdef CViewND   *ViewND_FromState          ( CInteger     rank       ,
                                                CInteger     offset     ,
                                                CInteger     size       ,
                                                CStatus     *status     )
    cdef CInteger   ViewND_GetExtent          ( CViewND     *self       ,
                                                CInteger     dimension  ,
                                                CStatus     *status     )
    cdef CInteger   ViewND_GetIndexMultiSlice ( CViewND     *self       ,
                                                CMultiSlice *multiSlice ,
                                                CStatus     *status     )
    cdef CInteger   ViewND_GetOffset          ( CViewND     *self       )
    cdef CInteger   ViewND_GetRank            ( CViewND     *self       )
    cdef CInteger   ViewND_GetSize            ( CViewND     *self       )
    cdef CInteger   ViewND_GetStride          ( CViewND     *self       ,
                                                CInteger     dimension  ,
                                                CStatus     *status     )
    cdef CIterator *ViewND_MakeIterator       ( CViewND     *self       ,
                                                CStatus     *status     )
    cdef void       ViewND_SetExtentStride    ( CViewND     *self       ,
                                                CInteger     dimension  ,
                                                CInteger     extent     ,
                                                CInteger     stride     ,
                                                CStatus     *status     )
    cdef void       ViewND_View1DMultiSlice   ( CViewND     *self       ,
                                                CMultiSlice *multiSlice ,
                                                CView1D     *view       ,
                                                CStatus     *status     )
    cdef void       ViewND_View2DMultiSlice   ( CViewND     *self       ,
                                                CMultiSlice *multiSlice ,
                                                CView2D     *view       ,
                                                CStatus     *status     )
    cdef CViewND   *ViewND_ViewMultiSlice     ( CViewND     *self       ,
                                                CMultiSlice *multiSlice ,
                                                CStatus     *status     )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class BaseArrayND:

    cdef CMultiSlice   *cMultiSlice
    cdef CViewND       *cView
    cdef public object  block
    cdef        object _iterator

    cdef void _DeallocateCObject ( self )
    cdef void _MakeCObject       ( self )
