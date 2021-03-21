from pCore.CPrimitiveTypes                 cimport CBoolean, CFalse, CInteger, CReal, CTrue
from pCore.Status                          cimport CStatus, CStatus_OK
from pScientific.Geometry3.Transformation3 cimport CTransformation3           , \
                                                   Transformation3            , \
                                                   Transformation3_IsEqual    , \
                                                   Transformation3_IsIdentity

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "Transformation3Container.h":

    ctypedef struct CTransformation3Container "Transformation3Container":
        CBoolean           isOwner
        CInteger           capacity
        CInteger           identity
        CInteger          *inverses
        CTransformation3 **items

    cdef CTransformation3Container *Transformation3Container_Allocate                ( CInteger                    capacity ,
                                                                                       CStatus                    *status   )
    cdef void                       Transformation3Container_Deallocate              ( CTransformation3Container **self     )
    cdef void                       Transformation3Container_FindIdentity            ( CTransformation3Container  *self     )
    cdef void                       Transformation3Container_FindInverses            ( CTransformation3Container  *self     ,
                                                                                       CStatus                    *status   )
    cdef CInteger                   Transformation3Container_NumberOfNonSelfInverses ( CTransformation3Container  *self     )
    cdef CInteger                   Transformation3Container_NumberOfSelfInverses    ( CTransformation3Container  *self     )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class Transformation3Container:

    cdef CTransformation3Container *cObject
    cdef public object              isOwner
    cdef public object              items
