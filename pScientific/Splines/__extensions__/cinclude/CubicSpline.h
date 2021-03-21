# ifndef _CUBICSPLINE
# define _CUBICSPLINE

# include "Boolean.h"
# include "Integer.h"
# include "Real.h"
# include "RealArray1D.h"
# include "RealArray2D.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The (multi-) cubic spline type. */
/* . L = number of points in spline >= 2 ; N = number of splines. */
typedef struct {
    RealArray1D *x ; /* . X-values           (L x 1). */
    RealArray2D *y ; /* . Y-values           (L x N). */
    RealArray2D *h ; /* . Second derivatives (L x N). */
} CubicSpline ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define CubicSpline_FastEvaluateFGN( self, n, l, u, d, s, t, f, g ) \
    { \
        auto Real hl, hu, yl, yu ; \
        hl = Array2D_Item ( self->h, l, n ) * d / 6.0e+00 ; \
        hu = Array2D_Item ( self->h, u, n ) * d / 6.0e+00 ; \
        yl = Array2D_Item ( self->y, l, n ) ; \
        yu = Array2D_Item ( self->y, u, n ) ; \
        f  = t * yl + s * yu + d * ( t * ( t * t - 1.0e+00 ) * hl + s * ( s * s - 1.0e+00 ) * hu ) ; \
        g  = ( yu - yl ) / d + ( - ( 3.0e+00 * t * t - 1.0e+00 ) * hl + ( 3.0e+00 * s * s - 1.0e+00 ) * hu ) ; \
    }

# define CubicSpline_FastEvaluateLUDST( self, x0, l, u, d, s, t ) \
    { \
        auto Integer i ; \
        auto RealArray1D *abscissa = self->x ; \
        l = 0 ; \
        u = View1D_Extent ( abscissa ) - 1 ; \
        while ( ( u - l ) > 1 ) \
        { \
            i = ( u + l ) >> 1 ; \
            if ( Array1D_Item ( abscissa, i ) > x0 ) u = i ; \
            else                                     l = i ; \
        } \
        d = ( Array1D_Item ( abscissa, u ) - Array1D_Item ( abscissa, l ) ) ; \
        s = ( x0 - Array1D_Item ( abscissa, l ) ) / d ; \
        t = ( Array1D_Item ( abscissa, u ) - x0 ) / d ; \
    }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern CubicSpline *CubicSpline_Allocate            (       Status       *status          ) ;
extern CubicSpline *CubicSpline_AllocateWithExtents ( const Integer       points          ,
                                                      const Integer       splines         ,
                                                            Status       *status          ) ;
extern  void        CubicSpline_AssignArrays        (       CubicSpline  *self            ,
                                                            RealArray1D  *x               ,
                                                            RealArray2D  *y               ,
                                                            RealArray2D  *h               ,
                                                            Status       *status          ) ;
extern void         CubicSpline_CheckXYH            (       CubicSpline  *self            ,
                                                            Status       *status          ) ;
extern CubicSpline *CubicSpline_Clone               ( const CubicSpline  *self            ,
                                                            Status       *status          ) ;
extern void         CubicSpline_Deallocate          (       CubicSpline **self            ) ;
extern void         CubicSpline_DeassignArrays      (       CubicSpline  *self            ) ;
extern void         CubicSpline_Evaluate            ( const CubicSpline  *self            ,
                                                      const Integer       spline          ,
                                                      const Real          x               ,
                                                            Real         *f               ,
                                                            Real         *g               ,
                                                            Real         *h               ,
                                                            Status       *status          ) ;
extern void         CubicSpline_EvaluateLUDST       ( const CubicSpline  *self            ,
                                                      const Real          x               ,
                                                            Integer      *l               ,
                                                            Integer      *u               ,
                                                            Real         *d               ,
                                                            Real         *s               ,
                                                            Real         *t               ) ;
extern void         CubicSpline_FindExtrema         ( const CubicSpline  *self            ,
                                                      const Integer       spline          ,
                                                            RealArray1D  *maxima          ,
                                                            RealArray1D  *minima          ,
                                                            Integer      *nMaxima         ,
                                                            Integer      *nMinima         ,
                                                            Status       *status          ) ;
extern CubicSpline *CubicSpline_FromRealArrays      (       RealArray1D  *x               ,
                                                            RealArray2D  *y               ,
                                                            RealArray2D  *h               ,
                                                      const Integer       lowerDerivative ,
                                                      const Real          lowerValue      ,
                                                      const Integer       upperDerivative ,
                                                      const Real          upperValue      ,
                                                            Status       *status          ) ;
extern void         CubicSpline_Initialize          (       CubicSpline  *self            ) ;
extern Real         CubicSpline_Integrate           ( const CubicSpline  *self            ,
                                                      const Integer       spline          ,
                                                      const Real          a               ,
                                                      const Real          b               ,
                                                            Status       *status          ) ;
extern Real         CubicSpline_IntegrateFull       ( const CubicSpline  *self            ,
                                                      const Integer       spline          ,
                                                            Status       *status          ) ;
extern void         CubicSpline_SetUpSpline         (       CubicSpline  *self            ,  
                                                      const Integer       spline          ,  
                                                      const Integer       lowerDerivative ,  
                                                      const Real          lowerValue      ,  
                                                      const Integer       upperDerivative ,  
                                                      const Real          upperValue      ,  
                                                            Status       *status          ) ;
# endif

