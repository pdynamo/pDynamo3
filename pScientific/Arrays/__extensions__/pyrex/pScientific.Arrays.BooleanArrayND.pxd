from pCore.CPrimitiveTypes             cimport CBoolean            , \
                                               CFalse              , \
                                               CInteger            , \
                                               CReal               , \
                                               CTrue
from pCore.BooleanBlock                cimport CBooleanBlock , \
                                               BooleanBlock
from pCore.Status                      cimport CStatus             , \
                                               CStatus_OK
from pScientific.Arrays.BaseArrayND    cimport BaseArrayND         , \
                                               CViewND
from pScientific.Arrays.BooleanArray1D cimport BooleanArray1D
from pScientific.Arrays.BooleanArray2D cimport BooleanArray2D
from pScientific.Arrays.Slicing        cimport CMultiSlice         , \
                                               MultiSlice_Allocate

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "BooleanArrayND.h":

    # . The array type.
    ctypedef struct CBooleanArrayND "BooleanArrayND":
        CViewND       *view
        CBooleanBlock *block
        CBoolean      *data

    # . Functions.
    cdef void             BooleanArrayND_CopyTo            ( CBooleanArrayND  *self          ,
                                                             CBooleanArrayND  *other         ,
                                                             CStatus          *status        )
    cdef void             BooleanArrayND_Deallocate        ( CBooleanArrayND **self          )
    cdef CBooleanArrayND *BooleanArrayND_FromViewBlock     ( CViewND          *view          ,
                                                             CBooleanBlock    *block         ,
                                                             CBoolean          withReference ,
                                                             CStatus          *status        )
    cdef CBoolean         BooleanArrayND_GetItemMultiSlice ( CBooleanArrayND  *self          ,
                                                             CMultiSlice      *multiSlice    ,
                                                             CStatus          *status        )
    cdef void             BooleanArrayND_Set               ( CBooleanArrayND  *self          ,
                                                             CBoolean          value         ,
                                                             CStatus          *status        )
    cdef void             BooleanArrayND_SetItemMultiSlice ( CBooleanArrayND  *self          ,
                                                             CMultiSlice      *multiSlice    ,
                                                             CBoolean          value         ,
                                                             CStatus          *status        )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class BooleanArrayND ( BaseArrayND ):

    cdef CBooleanArrayND *cObject
