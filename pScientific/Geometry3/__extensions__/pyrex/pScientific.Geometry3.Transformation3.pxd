from pCore.CPrimitiveTypes          cimport CBoolean     , \
                                            CFalse       , \
                                            CInteger     , \
                                            CReal        , \
                                            CTrue
from pCore.Status                   cimport CStatus      , \
                                            CStatus_OK
from pScientific.Arrays.RealArray1D cimport CRealArray1D , \
                                            RealArray1D
from pScientific.Arrays.RealArray2D cimport CRealArray2D , \
                                            RealArray2D
from pScientific.Geometry3.Matrix33 cimport Matrix33
from pScientific.Geometry3.Vector3  cimport Vector3

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "Transformation3.h":

    ctypedef struct CTransformation3 "Transformation3":
        CRealArray2D *rotation
        CRealArray1D *translation

    cdef CTransformation3 *Transformation3_Allocate                ( CStatus           *status      )
    cdef CTransformation3 *Transformation3_Clone                   ( CTransformation3  *self        ,
                                                                     CStatus           *status      )
    cdef void              Transformation3_CopyTo                  ( CTransformation3  *self        ,
                                                                     CTransformation3  *other       )
    cdef void              Transformation3_Deallocate              ( CTransformation3 **self        )
    cdef CTransformation3 *Transformation3_FromRotationTranslation ( CRealArray2D      *rotation    ,
                                                                     CRealArray1D      *translation ,
                                                                     CStatus           *status      )
    cdef CBoolean          Transformation3_IsEqual                 ( CTransformation3  *self        ,
                                                                     CTransformation3  *other       )
    cdef CBoolean          Transformation3_IsIdentity              ( CTransformation3  *self        )
    cdef void              Transformation3_Orthogonalize           ( CTransformation3  *self        ,
                                                                     CRealArray2D      *A           ,
                                                                     CRealArray2D      *B           )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class Transformation3:

    cdef CTransformation3 *cObject
    cdef public object     isOwner
    cdef public object     rotation
    cdef public object     translation
