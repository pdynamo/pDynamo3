from pCore.CPrimitiveTypes          cimport CBoolean            , \
                                            CFalse              , \
                                            CInteger            , \
                                            CReal               , \
                                            CTrue
from pCore.RealBlock                cimport CRealBlock          , \
                                            RealBlock
from pCore.Status                   cimport CStatus             , \
                                            CStatus_OK
from pScientific.Arrays.BaseArrayND cimport BaseArrayND         , \
                                            CViewND
from pScientific.Arrays.RealArray1D cimport RealArray1D
from pScientific.Arrays.RealArray2D cimport RealArray2D
from pScientific.Arrays.Slicing     cimport CMultiSlice         , \
                                            MultiSlice_Allocate

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "RealArrayND.h":

    # . The array type.
    ctypedef struct CRealArrayND "RealArrayND":
        CViewND    *view
        CRealBlock *block
        CReal      *data

    # . Functions.
    cdef void          RealArrayND_CopyTo            ( CRealArrayND  *self          ,
                                                       CRealArrayND  *other         ,
                                                       CStatus       *status        )
    cdef void          RealArrayND_Deallocate        ( CRealArrayND **self          )
    cdef CRealArrayND *RealArrayND_FromViewBlock     ( CViewND       *view          ,
                                                       CRealBlock    *block         ,
                                                       CBoolean       withReference ,
                                                       CStatus       *status        )
    cdef CReal         RealArrayND_GetItemMultiSlice ( CRealArrayND  *self          ,
                                                       CMultiSlice   *multiSlice    ,
                                                       CStatus       *status        )
    cdef void          RealArrayND_Set               ( CRealArrayND  *self          ,
                                                       CReal          value         ,
                                                       CStatus       *status        )
    cdef void          RealArrayND_SetItemMultiSlice ( CRealArrayND  *self          ,
                                                       CMultiSlice   *multiSlice    ,
                                                       CReal          value         ,
                                                       CStatus       *status        )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class RealArrayND ( BaseArrayND ):

    cdef CRealArrayND *cObject
