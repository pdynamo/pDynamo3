from pCore.CPrimitiveTypes          cimport CBoolean            , \
                                            CFalse              , \
                                            CInteger            , \
                                            CReal               , \
                                            CTrue
from pCore.Status                   cimport CStatus             , \
                                            CStatus_OK
from pScientific.Arrays.RealArray1D cimport CRealArray1D        , \
                                            RealArray1D         , \
                                            RealArray1D_GetItem
from pScientific.Arrays.RealArray2D cimport CRealArray2D        , \
                                            RealArray2D
from pScientific.Geometry3.Vector3  cimport Vector3

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "Matrix33.h":

#    ctypedef struct CRealArray2D "Matrix33":
#        pass

    cdef void     Matrix33_ApplyToVector3         ( CRealArray2D  *self    ,
                                                    CRealArray1D  *vector3 )
    cdef void     Matrix33_Invert                 ( CRealArray2D  *self    ,
                                                    CRealArray2D  *other   )
    cdef CBoolean Matrix33_IsEqual                ( CRealArray2D  *self    ,
                                                    CRealArray2D  *other   )
    cdef CBoolean Matrix33_IsIdentity             ( CRealArray2D  *self    )
    cdef CBoolean Matrix33_IsImproperRotation     ( CRealArray2D  *self    )
    cdef CBoolean Matrix33_IsOrthogonal           ( CRealArray2D  *self    )
    cdef CBoolean Matrix33_IsProperRotation       ( CRealArray2D  *self    )
    cdef void     Matrix33_PostMultiplyBy         ( CRealArray2D  *self    ,
                                                    CRealArray2D  *other   )
    cdef void     Matrix33_PreMultiplyBy          ( CRealArray2D  *self    ,
                                                    CRealArray2D  *other   )
    cdef void     Matrix33_Reflection             ( CRealArray2D **self    ,
                                                    CRealArray1D  *normal  )
    cdef CStatus  Matrix33_RotationAboutAxis      ( CRealArray2D **self    ,
                                                    CReal          angle   ,
                                                    CReal          x       ,
                                                    CReal          y       ,
                                                    CReal          z       )
    cdef CStatus  Matrix33_RotationFromQuaternion ( CRealArray2D **self    ,
                                                    CReal          q0      ,
                                                    CReal          q1      ,
                                                    CReal          q2      ,
                                                    CReal          q3      )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class Matrix33 ( RealArray2D ):
    pass
