# ifndef _BICUBICSPLINE
# define _BICUBICSPLINE

# include "Boolean.h"
# include "Integer.h"
# include "Real.h"
# include "RealArray1D.h"
# include "RealArray2D.h"
# include "RealArrayND.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Spline boundary conditions. */
typedef enum {
    BicubicSplineType_Clamped  = 0,
    BicubicSplineType_Natural  = 1,
    BicubicSplineType_NotAKnot = 2,
    BicubicSplineType_Periodic = 3
} BicubicSplineType ;

/* . The bicubic spline type. */
typedef struct {
    BicubicSplineType type    ;
    Integer           lengthX ; /* . Number of points >= 2. Number of intervals is length - 1. */
    Integer           lengthY ;
    RealArray1D      *x       ; /* . X-values. */
    RealArray1D      *y       ; /* . Y-values. */
    RealArray2D      *f       ; /* . Function values. */
    RealArrayND      *coefficients ; /* . Coefficients required for evaluating the spline. */
} BicubicSpline ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern BicubicSpline *BicubicSpline_Allocate            ( const Integer lengthx, const Integer lengthy, const Boolean doX, const Boolean doY, const Boolean doF, Status *status ) ;
extern BicubicSpline *BicubicSpline_Clone               ( const BicubicSpline  *self, Status *status ) ;
extern void           BicubicSpline_Deallocate          (       BicubicSpline **self ) ;
extern void           BicubicSpline_Evaluate            ( const BicubicSpline  *self, const Real x, const Real y, Real *f, Real *g1, Real *g2, Status *status ) ;
extern BicubicSpline *BicubicSpline_MakeFromRealArray2D ( RealArray1D **x, RealArray1D **y, RealArray2D **f, const BicubicSplineType type, Status *status ) ;

# endif
