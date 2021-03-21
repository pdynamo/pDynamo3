from pCore.CPrimitiveTypes          cimport CBoolean    , \
                                            CFalse      , \
                                            CInteger    , \
                                            CReal       , \
                                            CTrue
from pCore.RealBlock                cimport CRealBlock  , \
                                            RealBlock
from pCore.Status                   cimport CStatus     , \
                                            CStatus_OK
from pScientific.Arrays.BaseArray1D cimport BaseArray1D , \
                                            CView1D

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "RealArray1D.h":

    # . The array type.
    ctypedef struct CRealArray1D "RealArray1D":
        pass

    # . Functions.
    cdef CRealArray1D *RealArray1D_Allocate           ( CStatus       *status        )
    cdef CRealArray1D *RealArray1D_AllocateWithExtent ( CInteger       extent        ,
                                                        CStatus       *status        )
    cdef void          RealArray1D_AssignBlock        ( CRealArray1D  *self          ,
                                                        CRealBlock    *block         ,
                                                        CBoolean       withReference ,
                                                        CStatus       *status        )
    cdef void          RealArray1D_CopyTo             ( CRealArray1D  *self          ,
                                                        CRealArray1D  *other         ,
                                                        CStatus       *status        )
    cdef void          RealArray1D_Deallocate         ( CRealArray1D **self          )
    cdef CReal         RealArray1D_GetItem            ( CRealArray1D  *self          ,
                                                        CInteger       index         ,
                                                        CStatus       *status        )
    cdef CReal        *RealArray1D_PointerToData      ( CRealArray1D  *self          )
    cdef void          RealArray1D_Scale              ( CRealArray1D  *self          ,
                                                        CReal          value         )
    cdef void          RealArray1D_Set                ( CRealArray1D  *self          ,
                                                        CReal          value         )
    cdef void          RealArray1D_SetItem            ( CRealArray1D  *self          ,
                                                        CInteger       index         ,
                                                        CReal          value         ,
                                                        CStatus       *status        )
    cdef void          RealArray1D_Sort               ( CRealArray1D  *self          )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class RealArray1D ( BaseArray1D ):

    cdef CRealArray1D *cObject
