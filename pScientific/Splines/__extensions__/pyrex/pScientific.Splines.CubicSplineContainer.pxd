from pCore.CPrimitiveTypes           cimport CBoolean     , \
                                             CFalse       , \
                                             CInteger     , \
                                             CReal        , \
                                             CTrue           
from pCore.Status                    cimport CStatus      , \
                                             CStatus_OK      
from pScientific.Splines.CubicSpline cimport CCubicSpline , \
                                             CubicSpline

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "CubicSplineContainer.h":

    ctypedef struct CCubicSplineContainer "CubicSplineContainer":
        CBoolean       isOwner 
        CInteger       capacity
        CCubicSpline **entries 

    cdef CCubicSplineContainer *CubicSplineContainer_Allocate   ( CInteger                capacity ,
                                                                  CStatus                *status   )
    cdef CCubicSplineContainer *CubicSplineContainer_Clone      ( CCubicSplineContainer  *self     ,
                                                                  CStatus                *status   )
    cdef void                   CubicSplineContainer_Deallocate ( CCubicSplineContainer **self     )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class CubicSplineContainer:

    cdef CCubicSplineContainer *cObject
    cdef public object          isOwner
    cdef public object          keys
    cdef public object          label
    cdef public object          uniqueEntries
