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

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "CubicSpline.h":

    # . The cubic spline type.
    ctypedef struct CCubicSpline "CubicSpline":
        CRealArray1D *x
        CRealArray2D *y
        CRealArray2D *h

    # . Functions.
    cdef CCubicSpline *CubicSpline_Allocate            ( CStatus       *status          )
    cdef CCubicSpline *CubicSpline_AllocateWithExtents ( CInteger       points          ,
                                                         CInteger       splines         ,
                                                         CStatus       *status          )
    cdef  void         CubicSpline_AssignArrays        ( CCubicSpline  *self            ,
                                                         CRealArray1D  *x               ,
                                                         CRealArray2D  *y               ,
                                                         CRealArray2D  *h               ,
                                                         CStatus       *status          )
    cdef void          CubicSpline_CheckXYH            ( CCubicSpline  *self            ,
                                                         CStatus       *status          )
    cdef CCubicSpline *CubicSpline_Clone               ( CCubicSpline  *self            ,
                                                         CStatus       *status          )
    cdef void          CubicSpline_Deallocate          ( CCubicSpline **self            )
    cdef void          CubicSpline_DeassignArrays      ( CCubicSpline  *self            )
    cdef void          CubicSpline_Evaluate            ( CCubicSpline  *self            ,
                                                         CInteger       spline          ,
                                                         CReal          x               ,
                                                         CReal         *f               ,
                                                         CReal         *g               ,
                                                         CReal         *h               ,
                                                         CStatus       *status          )
    cdef void          CubicSpline_EvaluateLUDST       ( CCubicSpline  *self            ,
                                                         CReal          x               ,
                                                         CInteger      *l               ,
                                                         CInteger      *u               ,
                                                         CReal         *d               ,
                                                         CReal         *s               ,
                                                         CReal         *t               )
    cdef void          CubicSpline_FindExtrema         ( CCubicSpline  *self            ,
                                                         CInteger       spline          ,
                                                         CRealArray1D  *maxima          ,
                                                         CRealArray1D  *minima          ,
                                                         CInteger      *nMaxima         ,
                                                         CInteger      *nMinima         ,
                                                         CStatus       *status          )
    cdef CCubicSpline *CubicSpline_FromRealArrays      ( CRealArray1D  *x               ,
                                                         CRealArray2D  *y               ,
                                                         CRealArray2D  *h               ,
                                                         CInteger       lowerDerivative ,
                                                         CReal          lowerValue      ,
                                                         CInteger       upperDerivative ,
                                                         CReal          upperValue      ,
                                                         CStatus       *status          )
    cdef void          CubicSpline_Initialize          ( CCubicSpline  *self            )
    cdef CReal         CubicSpline_Integrate           ( CCubicSpline  *self            ,
                                                         CInteger       spline          ,
                                                         CReal          a               ,
                                                         CReal          b               ,
                                                         CStatus       *status          )
    cdef CReal         CubicSpline_IntegrateFull       ( CCubicSpline  *self            ,
                                                         CInteger       spline          ,
                                                         CStatus       *status          )
    cdef void          CubicSpline_SetUpSpline         ( CCubicSpline  *self            ,
                                                         CInteger       spline          ,
                                                         CInteger       lowerDerivative ,
                                                         CReal          lowerValue      ,
                                                         CInteger       upperDerivative ,
                                                         CReal          upperValue      ,
                                                         CStatus       *status          )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class CubicSpline:

    cdef CCubicSpline      *cObject
    cdef public object      isOwner
    cdef public RealArray1D x
    cdef public RealArray2D y
    cdef public RealArray2D h 
