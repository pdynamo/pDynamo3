from pCore.CPrimitiveTypes             cimport CBoolean       , \
                                               CFalse         , \
                                               CInteger       , \
                                               CReal          , \
                                               CTrue
from pCore.BooleanBlock                cimport CBooleanBlock  , \
                                               BooleanBlock
from pCore.Status                      cimport CStatus        , \
                                               CStatus_OK
from pScientific.Arrays.BaseArray2D    cimport BaseArray2D    , \
                                               CView2D
from pScientific.Arrays.BooleanArray1D cimport BooleanArray1D
from pScientific.Arrays.Slicing        cimport CMultiSlice

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "BooleanArray2D.h":

    # . The array type.
    ctypedef struct CBooleanArray2D "BooleanArray2D":
        pass

    # . Functions.
    cdef CBooleanArray2D *BooleanArray2D_Allocate            ( CStatus          *status        )
    cdef CBooleanArray2D *BooleanArray2D_AllocateWithExtents ( CInteger          rows          ,
                                                               CInteger          columns       ,
                                                               CStatus          *status        )
    cdef void             BooleanArray2D_AssignBlock         ( CBooleanArray2D  *self          ,
                                                               CBooleanBlock    *block         ,
                                                               CBoolean          withReference ,
                                                               CStatus          *status        )
    cdef void             BooleanArray2D_CopyTo              ( CBooleanArray2D  *self          ,
                                                               CBooleanArray2D  *other         ,
                                                               CStatus          *status        )
    cdef void             BooleanArray2D_Deallocate          ( CBooleanArray2D **self          )
    cdef CBoolean         BooleanArray2D_GetItem             ( CBooleanArray2D  *self          ,
                                                               CInteger          i             ,
                                                               CInteger          j             ,
                                                               CStatus          *status        )
    cdef CBoolean         BooleanArray2D_GetItemMultiSlice   ( CBooleanArray2D  *self          ,
                                                               CMultiSlice      *multiSlice    ,
                                                               CStatus          *status        )
    cdef CBoolean        *BooleanArray2D_PointerToData       ( CBooleanArray2D  *self          )
    cdef void             BooleanArray2D_Set                 ( CBooleanArray2D  *self          ,
                                                               CBoolean          value         )
    cdef void             BooleanArray2D_SetItem             ( CBooleanArray2D  *self          ,
                                                               CInteger          i             ,
                                                               CInteger          j             ,
                                                               CBoolean          value         ,
                                                               CStatus          *status        )
    cdef void             BooleanArray2D_SetItemMultiSlice   ( CBooleanArray2D  *self          ,
                                                               CMultiSlice      *multiSlice    ,
                                                               CBoolean          value         ,
                                                               CStatus          *status        )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class BooleanArray2D ( BaseArray2D ):

    cdef CBooleanArray2D *cObject
