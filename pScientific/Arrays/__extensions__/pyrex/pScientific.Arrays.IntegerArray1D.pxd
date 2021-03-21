from pCore.CPrimitiveTypes          cimport CBoolean      , \
                                            CFalse        , \
                                            CInteger      , \
                                            CReal         , \
                                            CTrue
from pCore.IntegerBlock             cimport CIntegerBlock , \
                                            IntegerBlock
from pCore.Status                   cimport CStatus       , \
                                            CStatus_OK
from pScientific.Arrays.BaseArray1D cimport BaseArray1D   , \
                                            CView1D

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "IntegerArray1D.h":

    # . The array type.
    ctypedef struct CIntegerArray1D "IntegerArray1D":
        pass

    # . Functions.
    cdef CIntegerArray1D *IntegerArray1D_Allocate           ( CStatus          *status        )
    cdef CIntegerArray1D *IntegerArray1D_AllocateWithExtent ( CInteger          extent        ,
                                                              CStatus          *status        )  
    cdef void             IntegerArray1D_AssignBlock        ( CIntegerArray1D  *self          ,
                                                              CIntegerBlock    *block         ,
                                                              CBoolean          withReference ,
                                                              CStatus          *status        )
    cdef void             IntegerArray1D_CopyTo             ( CIntegerArray1D  *self          ,
                                                              CIntegerArray1D  *other         ,
                                                              CStatus          *status        )
    cdef void             IntegerArray1D_Deallocate         ( CIntegerArray1D **self          )
    cdef CInteger         IntegerArray1D_GetItem            ( CIntegerArray1D  *self          ,
                                                              CInteger          index         ,
                                                              CStatus          *status        )
    cdef CInteger        *IntegerArray1D_PointerToData      ( CIntegerArray1D  *self          )
    cdef void             IntegerArray1D_Set                ( CIntegerArray1D  *self          ,
                                                              CInteger          value         )
    cdef void             IntegerArray1D_SetItem            ( CIntegerArray1D  *self          ,
                                                              CInteger          index         ,
                                                              CInteger          value         ,
                                                              CStatus          *status        )
    cdef void             IntegerArray1D_Sort               ( CIntegerArray1D  *self          )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class IntegerArray1D ( BaseArray1D ):

    cdef CIntegerArray1D *cObject
