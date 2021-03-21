from pCore.CPrimitiveTypes             cimport CBoolean            , \
                                               CFalse              , \
                                               CInteger            , \
                                               CReal               , \
                                               CTrue
from pCore.IntegerBlock                cimport CIntegerBlock       , \
                                               IntegerBlock
from pCore.Status                      cimport CStatus             , \
                                               CStatus_OK
from pScientific.Arrays.BaseArrayND    cimport BaseArrayND         , \
                                               CViewND
from pScientific.Arrays.IntegerArray1D cimport IntegerArray1D
from pScientific.Arrays.IntegerArray2D cimport IntegerArray2D
from pScientific.Arrays.Slicing        cimport CMultiSlice         , \
                                               MultiSlice_Allocate

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "IntegerArrayND.h":

    # . The array type.
    ctypedef struct CIntegerArrayND "IntegerArrayND":
        CViewND       *view
        CIntegerBlock *block
        CInteger      *data

    # . Functions.
    cdef void             IntegerArrayND_CopyTo            ( CIntegerArrayND  *self          ,
                                                             CIntegerArrayND  *other         ,
                                                             CStatus          *status        )
    cdef void             IntegerArrayND_Deallocate        ( CIntegerArrayND **self          )
    cdef CIntegerArrayND *IntegerArrayND_FromViewBlock     ( CViewND          *view          ,
                                                             CIntegerBlock    *block         ,
                                                             CBoolean          withReference ,
                                                             CStatus          *status        )
    cdef CInteger         IntegerArrayND_GetItemMultiSlice ( CIntegerArrayND  *self          ,
                                                             CMultiSlice      *multiSlice    ,
                                                             CStatus          *status        )
    cdef void             IntegerArrayND_Set               ( CIntegerArrayND  *self          ,
                                                             CInteger          value         ,
                                                             CStatus          *status        )
    cdef void             IntegerArrayND_SetItemMultiSlice ( CIntegerArrayND  *self          ,
                                                             CMultiSlice      *multiSlice    ,
                                                             CInteger          value         ,
                                                             CStatus          *status        )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class IntegerArrayND ( BaseArrayND ):

    cdef CIntegerArrayND *cObject
