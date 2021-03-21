from pCore.CPrimitiveTypes          cimport CBoolean      , \
                                            CFalse        , \
                                            CInteger      , \
                                            CReal         , \
                                            CTrue
from pCore.BooleanBlock             cimport CBooleanBlock , \
                                            BooleanBlock
from pCore.Status                   cimport CStatus       , \
                                            CStatus_OK
from pScientific.Arrays.BaseArray1D cimport BaseArray1D   , \
                                            CView1D

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "BooleanArray1D.h":

    # . The array type.
    ctypedef struct CBooleanArray1D "BooleanArray1D":
        pass

    # . Functions.
    cdef CBooleanArray1D *BooleanArray1D_Allocate           ( CStatus          *status        )
    cdef CBooleanArray1D *BooleanArray1D_AllocateWithExtent ( CInteger          extent        ,
                                                              CStatus          *status        )  
    cdef void             BooleanArray1D_AssignBlock        ( CBooleanArray1D  *self          ,
                                                              CBooleanBlock    *block         ,
                                                              CBoolean          withReference ,
                                                              CStatus          *status        )
    cdef void             BooleanArray1D_CopyTo             ( CBooleanArray1D  *self          ,
                                                              CBooleanArray1D  *other         ,
                                                              CStatus          *status        )
    cdef void             BooleanArray1D_Deallocate         ( CBooleanArray1D **self          )
    cdef CBoolean         BooleanArray1D_GetItem            ( CBooleanArray1D  *self          ,
                                                              CInteger          index         ,
                                                              CStatus          *status        )
    cdef CBoolean        *BooleanArray1D_PointerToData      ( CBooleanArray1D  *self          )
    cdef void             BooleanArray1D_Set                ( CBooleanArray1D  *self          ,
                                                              CBoolean          value         )
    cdef void             BooleanArray1D_SetItem            ( CBooleanArray1D  *self          ,
                                                              CInteger          index         ,
                                                              CBoolean          value         ,
                                                              CStatus          *status        )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class BooleanArray1D ( BaseArray1D ):

    cdef CBooleanArray1D *cObject
