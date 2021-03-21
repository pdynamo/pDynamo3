from pCore.CPrimitiveTypes          cimport CBoolean, CFalse, CInteger, CReal, CTrue
from pCore.Status                   cimport CStatus
from pScientific.Arrays.RealArray1D cimport CRealArray1D
from pScientific.Arrays.RealArray2D cimport CRealArray2D
from pScientific.Arrays.RealArrayND cimport CRealArrayND

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "CubicSpline.h":

    # . Spline boundary conditions.
    ctypedef enum BicubicSplineType:
        BicubicSplineType_Clamped  = 0
        BicubicSplineType_Natural  = 1
        BicubicSplineType_NotAKnot = 2
        BicubicSplineType_Periodic = 3

    # . The bicubic spline type.
    ctypedef struct CBicubicSpline "BicubicSpline":
        BicubicSplineType type
        CInteger           lengthX
        CInteger           lengthY
        CRealArray1D      *x
        CRealArray1D      *y
        CRealArray2D      *f
        CRealArrayND      *coefficients

    # . Functions.
    cdef CBicubicSpline *BicubicSpline_Allocate            ( CInteger lengthx, CInteger lengthy, CBoolean doX, CBoolean doY, CBoolean doF, CStatus *status )
    cdef CBicubicSpline *BicubicSpline_Clone               ( CBicubicSpline  *self, CStatus *status )
    cdef void            BicubicSpline_Deallocate          ( CBicubicSpline **self )
    cdef void            BicubicSpline_Evaluate            ( CBicubicSpline  *self, CReal x, CReal y, CReal *f, CReal *g1, CReal *g2, CStatus *status )
    cdef CBicubicSpline *BicubicSpline_MakeFromRealArray2D ( CRealArray1D **x, CRealArray1D **y, CRealArray2D **f, BicubicSplineType type, CStatus *status )
