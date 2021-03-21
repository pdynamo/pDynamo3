from pCore.CPrimitiveTypes             cimport CBoolean       , \
                                               CFalse         , \
                                               CInteger       , \
                                               CReal          , \
                                               CTrue
from pCore.IntegerBlock                cimport CIntegerBlock  , \
                                               IntegerBlock
from pCore.Status                      cimport CStatus        , \
                                               CStatus_OK
from pScientific.Arrays.BaseArray2D    cimport BaseArray2D    , \
                                               CView2D
from pScientific.Arrays.IntegerArray1D cimport IntegerArray1D
from pScientific.Arrays.Slicing        cimport CMultiSlice

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "IntegerArray2D.h":

    # . The array type.
    ctypedef struct CIntegerArray2D "IntegerArray2D":
        pass

    # . Functions.
    cdef CIntegerArray2D *IntegerArray2D_Allocate            ( CStatus          *status        )
    cdef CIntegerArray2D *IntegerArray2D_AllocateWithExtents ( CInteger          rows          ,
                                                               CInteger          columns       ,
                                                               CStatus          *status        )
    cdef void             IntegerArray2D_AssignBlock         ( CIntegerArray2D  *self          ,
                                                               CIntegerBlock    *block         ,
                                                               CBoolean          withReference ,
                                                               CStatus          *status        )
    cdef void             IntegerArray2D_CopyTo              ( CIntegerArray2D  *self          ,
                                                               CIntegerArray2D  *other         ,
                                                               CStatus          *status        )
    cdef void             IntegerArray2D_Deallocate          ( CIntegerArray2D **self          )
    cdef CInteger         IntegerArray2D_GetItem             ( CIntegerArray2D  *self          ,
                                                               CInteger          i             ,
                                                               CInteger          j             ,
                                                               CStatus          *status        )
    cdef CInteger         IntegerArray2D_GetItemMultiSlice   ( CIntegerArray2D  *self          ,
                                                               CMultiSlice      *multiSlice    ,
                                                               CStatus          *status        )
    cdef CInteger        *IntegerArray2D_PointerToData       ( CIntegerArray2D  *self          )
    cdef void             IntegerArray2D_Set                 ( CIntegerArray2D  *self          ,
                                                               CInteger          value         )
    cdef void             IntegerArray2D_SetItem             ( CIntegerArray2D  *self          ,
                                                               CInteger          i             ,
                                                               CInteger          j             ,
                                                               CInteger          value         ,
                                                               CStatus          *status        )
    cdef void             IntegerArray2D_SetItemMultiSlice   ( CIntegerArray2D  *self          ,
                                                               CMultiSlice      *multiSlice    ,
                                                               CInteger          value         ,
                                                               CStatus          *status        )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class IntegerArray2D ( BaseArray2D ):

    cdef CIntegerArray2D *cObject
