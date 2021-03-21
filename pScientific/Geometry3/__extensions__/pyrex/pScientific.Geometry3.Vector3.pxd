from pCore.CPrimitiveTypes          cimport CBoolean        , \
                                            CFalse          , \
                                            CInteger        , \
                                            CReal           , \
                                            CTrue
from pCore.Status                   cimport CStatus         , \
                                            CStatus_OK
from pScientific.Arrays.RealArray1D cimport CRealArray1D    , \
                                            RealArray1D

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "Vector3.h":

#    ctypedef struct CVector3 "Vector3":
#        pass

    cdef void Vector3_CrossProduct ( CRealArray1D *self, CRealArray1D *other )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class Vector3 ( RealArray1D ):
    pass
