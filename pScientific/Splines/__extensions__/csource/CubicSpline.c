/*==================================================================================================================================
! . Cubic splines.
! . Follows the symmetrical formulation in Numerical Recipes.
!=================================================================================================================================*/

# include <stdlib.h>

# include "CubicSpline.h"
# include "f2clapack.h"
# include "math.h"
# include "Memory.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Basic allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
CubicSpline *CubicSpline_Allocate ( Status *status )
{
    CubicSpline *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( CubicSpline ) ;
        if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
        else                CubicSpline_Initialize ( self ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation with extents.
!---------------------------------------------------------------------------------------------------------------------------------*/
CubicSpline *CubicSpline_AllocateWithExtents ( const Integer points  ,
                                               const Integer splines ,
                                                     Status *status  )
{
    CubicSpline *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        if ( ( points >= 2 ) && ( splines > 0 ) )
        {
            self    = CubicSpline_Allocate ( status ) ;
            self->x = RealArray1D_AllocateWithExtent  ( points,          status ) ;
            self->y = RealArray2D_AllocateWithExtents ( points, splines, status ) ;
            self->h = RealArray2D_AllocateWithExtents ( points, splines, status ) ;
            if ( ! Status_IsOK ( status ) ) CubicSpline_Deallocate ( &self ) ;
        }    
        else Status_Set ( status, Status_InvalidArgument ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Assign arrays - with dimension checks.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CubicSpline_AssignArrays ( CubicSpline *self    ,
                                RealArray1D *x       ,
                                RealArray2D *y       ,
                                RealArray2D *h       ,
                                Status      *status  )
{
    if ( Status_IsOK ( status ) )
    {
        if ( ( x != NULL ) &&
             ( y != NULL ) &&
             ( h != NULL ) &&
             ( View1D_Extent  ( x ) >  1                    ) &&
             ( View2D_Columns ( y ) >  0                    ) &&
             ( View1D_Extent  ( x ) == View2D_Rows    ( y ) ) &&
             ( View1D_Extent  ( x ) == View2D_Rows    ( h ) ) &&
             ( View2D_Columns ( y ) == View2D_Columns ( h ) ) )
        {
            self->x = x ;
            self->y = y ;
            self->h = h ;
        }
        else Status_Set ( status, Status_InvalidArgument ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check the assigned arrays - their dimensions should already be OK.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CubicSpline_CheckXYH ( CubicSpline *self, Status *status )
{
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer n = View1D_Extent ( self->x ), i ;
        /* . Reverse x (and y) if necessary so that x is increasing. */
        if ( Array1D_Item ( self->x, 0 ) > Array1D_Item ( self->x, n-1 ) )
        {
            auto Integer     s ;
            auto RealArray1D view ;
            RealArray1D_Reverse ( self->x ) ;
            for ( s = 0 ; s < View2D_Columns ( self->y ) ; s++ )
            {
                RealArray2D_ColumnView ( self->y, s, False, &view, NULL ) ;
                RealArray1D_Reverse    ( &view ) ;
            }
        }
        /* . Initialize H. */
        RealArray2D_Set ( self->h, 0.0e+00 ) ;
        /* . Check to see if x is always increasing and the intervals are non-zero. */
        for ( i = 0 ; i < n-1 ; i++ )
        {
            if ( Array1D_Item ( self->x, i+1 ) <= Array1D_Item ( self->x, i ) )
            {
                Status_Set ( status, Status_InvalidArgument ) ; break ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
CubicSpline *CubicSpline_Clone ( const CubicSpline *self, Status *status )
{
    CubicSpline *clone = NULL ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        clone = CubicSpline_AllocateWithExtents ( View2D_Rows ( self->y ), View2D_Columns ( self->y ), status ) ;
        if ( clone != NULL )
        {
            RealArray1D_CopyTo ( self->x, clone->x, NULL ) ;
            RealArray2D_CopyTo ( self->h, clone->h, NULL ) ;
            RealArray2D_CopyTo ( self->y, clone->y, NULL ) ;
        }
        else Status_Set ( status, Status_OutOfMemory ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CubicSpline_Deallocate ( CubicSpline **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        RealArray1D_Deallocate ( &((*self)->x) ) ;
        RealArray2D_Deallocate ( &((*self)->y) ) ;
        RealArray2D_Deallocate ( &((*self)->h) ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deassign arrays.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CubicSpline_DeassignArrays ( CubicSpline *self ) { CubicSpline_Initialize ( self ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Evaluation - with checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CubicSpline_Evaluate ( const CubicSpline *self   ,
                            const Integer      spline ,
                            const Real         x      ,
                                  Real        *f      ,
                                  Real        *g      ,
                                  Real        *h      ,
                                  Status      *status )
{
    if ( f != NULL ) (*f) = 0.0e+00 ;
    if ( g != NULL ) (*g) = 0.0e+00 ;
    if ( h != NULL ) (*h) = 0.0e+00 ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto Boolean      isOK ;
        auto Integer      n = View1D_Extent ( self->x ) ;
        auto RealArray1D *abscissa = self->x ;
        isOK = ( ( spline >= 0                              ) &&
                 ( spline <  View2D_Columns ( self->y       ) ) &&
                 ( x      >= Array1D_Item   ( abscissa, 0   ) ) &&
                 ( x      <= Array1D_Item   ( abscissa, n-1 ) ) ) ;
        if ( isOK )
        {
            auto Integer i, l, u ;
            auto Real    hl, hu, d, s, t, yl, yu ;
            /* . Locate the index. */
            l = 0     ;
            u = n - 1 ;
            while ( ( u - l ) > 1 )
            {
                i = ( u + l ) >> 1 ;
                if ( Array1D_Item ( abscissa, i ) > x ) u = i ;
                else                                    l = i ;
            }
            /* . Factors. */
            d  = ( Array1D_Item ( abscissa, u ) - Array1D_Item ( abscissa, l ) ) ;
            s  = ( x - Array1D_Item ( abscissa, l ) ) / d ;
            t  = ( Array1D_Item ( abscissa, u ) - x ) / d ;
            hl = Array2D_Item ( self->h, l, spline ) * d / 6.0e+00 ;
            hu = Array2D_Item ( self->h, u, spline ) * d / 6.0e+00 ;
            yl = Array2D_Item ( self->y, l, spline ) ;
            yu = Array2D_Item ( self->y, u, spline ) ;
            /* . Calculate the function and its derivatives. */
            if ( f != NULL ) (*f) = t * yl + s * yu + d * ( t * ( t * t - 1.0e+00 ) * hl + s * ( s * s - 1.0e+00 ) * hu ) ;
            if ( g != NULL ) (*g) = ( yu - yl ) / d + ( - ( 3.0e+00 * t * t - 1.0e+00 ) * hl + ( 3.0e+00 * s * s - 1.0e+00 ) * hu ) ;
            if ( h != NULL ) (*h) = 6.0e+00 * ( t * hl + s * hu ) / d ;
        }
        else Status_Set ( status, Status_InvalidArgument ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Evaluation of quantities required to evaluate the spline - there is no checking!
!---------------------------------------------------------------------------------------------------------------------------------*/
void CubicSpline_EvaluateLUDST ( const CubicSpline *self ,
                                 const Real         x    ,
                                       Integer     *l    ,
                                       Integer     *u    ,
                                       Real        *d    ,
                                       Real        *s    ,
                                       Real        *t    )
{
    auto Integer      i, l0, u0 ;
    auto RealArray1D *abscissa = self->x ;
    /* . Locate the index. */
    l0 = 0 ;
    u0 = View1D_Extent ( abscissa ) - 1 ;
    while ( ( u0 - l0 ) > 1 )
    {
        i = ( u0 + l0 ) >> 1 ;
        if ( Array1D_Item ( abscissa, i ) > x ) u0 = i ;
        else                                    l0 = i ;
    }
    /* . Save data. */
    (*l) = l0 ;
    (*u) = u0 ;
    (*d) = ( Array1D_Item ( abscissa, u0 ) - Array1D_Item ( abscissa, l0 ) ) ;
    (*s) = ( x - Array1D_Item ( abscissa, l0 ) ) / (*d) ;
    (*t) = ( Array1D_Item ( abscissa, u0 ) - x ) / (*d) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find the positions of any extrema (maxima and minima only).
!---------------------------------------------------------------------------------------------------------------------------------*/
void CubicSpline_FindExtrema ( const CubicSpline *self    ,
                               const Integer      spline  ,
                                     RealArray1D *maxima  ,
                                     RealArray1D *minima  ,
                                     Integer     *nMaxima ,
                                     Integer     *nMinima ,
                                     Status      *status  )
{
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        if ( ( ( maxima  == NULL ) || ( ( maxima  != NULL ) && ( nMaxima != NULL ) ) ) &&
             ( ( minima  == NULL ) || ( ( minima  != NULL ) && ( nMinima != NULL ) ) ) &&
             ( ( nMaxima != NULL ) ||   ( nMinima != NULL ) ) &&
               ( spline >= 0                        ) &&
               ( spline <  View2D_Columns ( self->y ) ) )
        {
            auto Boolean      doMaxima, doMinima ;
            auto Integer      i, iTry, mMaxima = 0, mMinima = 0, n = View1D_Extent ( self->x ) , nTry = 0 ;
            auto Real         a, b, c, d, factor, h, hl, hu, s, sign, x, xl, xu, yl, yu ;
            auto RealArray1D *abscissa = self->x ;
            doMaxima = ( maxima != NULL ) ;
            doMinima = ( minima != NULL ) ;
            /* . Find the extrema. */
            for ( i = 0 ; i < n - 1 ; i++ )
            {
                /* . The interval range. */
                xl = Array1D_Item ( abscissa, i     ) ;
                xu = Array1D_Item ( abscissa, i + 1 ) ;
                d  = xu - xl ;
                /* . Other factors. */
                hl = Array2D_Item ( self->h, i    , spline ) * d / 6.0e+00 ;
                hu = Array2D_Item ( self->h, i + 1, spline ) * d / 6.0e+00 ;
                yl = Array2D_Item ( self->y, i    , spline ) ;
                yu = Array2D_Item ( self->y, i + 1, spline ) ;
                /* . Factors for a quadratic equation in terms of s. */
                a = 3.0e+00 * ( hu - hl ) ;
                b = 6.0e+00 * hl ;
                c = ( yu - yl ) / d - 2.0e+00 * hl - hu ;
                factor = b * b - 4.0e+00 * a * c ;
                /* . Solve the first derivative quadratic. */
                if ( factor >= 0.0e+00 )
                {
                    if ( factor == 0.0e+00 ) nTry = 1 ;
                    else { factor = sqrt ( factor ) ; nTry = 2 ; }
                    for ( iTry = 0, sign = -1.0e+00 ; iTry < nTry ; iTry++ )
                    {
                        s = ( - b + sign * factor ) / ( 2.0e+00 * a ) ;
                        /* . Only check for lower boundary extrema except for the last interval when both are checked. */
                        if ( ( s >= 0.0e+00 ) && ( ( ( s < 1.0e+00 ) && ( i < n - 2 ) ) || ( ( s <= 1.0e+00 ) && ( i == n - 2 ) ) ) )
                        {
                            x = d * s + xl ;
                            /* . Apply the second derivative test. */
                            /* . The more general extremum test is unnecessary for a cubic function because if the second derivative is zero
                            ! .  the point is an inflection point (when the third derivative is non-zero) or the function is a constant
                            ! . (as all derivatives are zero). */
                            CubicSpline_Evaluate ( self, spline, x, NULL, NULL, &h, status ) ;
                            if ( h > 0.0e+00 )
                            {
                                if ( doMinima )
                                {
                                    Array1D_Item ( minima, mMinima ) = x ;
                                    if ( mMinima >= View1D_Extent ( minima ) - 1 ) doMinima = False ;
                                }
                                mMinima ++ ;
                            }
                            else if ( h < 0.0e+00 )
                            {
                                if ( doMaxima )
                                {
                                    Array1D_Item ( maxima, mMaxima ) = x ;
                                    if ( mMaxima >= View1D_Extent ( maxima ) - 1 ) doMaxima = False ;
                                }
                                mMaxima ++ ;
                            }
                        }
                        sign *= -1.0e00 ;
                    }
                }
            }
            if ( nMaxima != NULL ) (*nMaxima) = mMaxima ;
            if ( nMinima != NULL ) (*nMinima) = mMinima ;
        }
        else Status_Set ( status, Status_InvalidArgument ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Spline constructor given appropriate arrays and boundary conditions (the same for all splines).
!---------------------------------------------------------------------------------------------------------------------------------*/
CubicSpline *CubicSpline_FromRealArrays (       RealArray1D  *x               ,
                                                RealArray2D  *y               ,
                                                RealArray2D  *h               ,
                                          const Integer       lowerDerivative ,
                                          const Real          lowerValue      ,
                                          const Integer       upperDerivative ,
                                          const Real          upperValue      ,
                                                Status       *status          )
{
    Integer s ;
    CubicSpline *self = CubicSpline_Allocate ( status ) ;
    CubicSpline_AssignArrays ( self, x, y, h, status ) ;
    CubicSpline_CheckXYH     ( self         , status ) ;
    for ( s = 0 ; s < View2D_Columns ( self->y ) ; s++ )
    {
        CubicSpline_SetUpSpline ( self, s, lowerDerivative, lowerValue, upperDerivative, upperValue, status ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CubicSpline_Initialize ( CubicSpline *self )
{
    if ( self != NULL )
    {
        self->h = NULL ;
        self->x = NULL ;
        self->y = NULL ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the integral of the spline in the range [a,b].
!---------------------------------------------------------------------------------------------------------------------------------*/
Real CubicSpline_Integrate ( const CubicSpline *self   ,
                             const Integer      spline ,
                             const Real         a      ,
                             const Real         b      ,
                                   Status      *status )
{
    Real integral = 0.0e+00 ;
    if ( ( self != NULL ) && ( a != b ) && Status_IsOK ( status ) )
    {
        auto Integer n = View1D_Extent ( self->x ) ;
        /* . Check the range of the integral and the spline index.*/
        if ( ( a < Array1D_Item ( self->x, 0 ) ) || ( a > Array1D_Item ( self->x, n-1 ) ) ||
             ( b < Array1D_Item ( self->x, 0 ) ) || ( b > Array1D_Item ( self->x, n-1 ) ) ||
             ( a > b ) || ( spline < 0 ) || ( spline >= View2D_Columns ( self->y ) ) ) Status_Set ( status, Status_InvalidArgument ) ;
        /* . OK. */
        else
        {
            Integer l, lA, lB, u ;
            Real    d, local, sA, sB, sF2, sF4, tA, tB, tF2, tF4 ;
            /* . Get data for a and b. */
            CubicSpline_EvaluateLUDST ( self, a, &lA, &u, &d, &sA, &tA ) ;
            CubicSpline_EvaluateLUDST ( self, b, &lB, &u, &d, &sB, &tB ) ;
            /* . Loop over intervals. */
            for ( l = lA ; l <= lB ; l++ )
            {
                u   = l + 1 ;
                d   = Array1D_Item ( self->x, u         ) - Array1D_Item ( self->x, l ) ;
                sF4 = Array2D_Item ( self->h, u, spline ) * d * d / 6.0e+00 ;
                sF2 = 2.0e+00 * ( Array2D_Item ( self->y, u, spline ) - sF4 ) ;
                tF4 = Array2D_Item ( self->h, l, spline ) * d * d / 6.0e+00 ;
                tF2 = 2.0e+00 * ( Array2D_Item ( self->y, l, spline ) - tF4 ) ;
                if ( l == lA ) local  = tA * tA * ( tF2 + tA * tA * tF4 ) - sA * sA * ( sF2 + sA * sA * sF4 ) ;
                else           local  = tF2 + tF4 ;
                if ( l == lB ) local += sB * sB * ( sF2 + sB * sB * sF4 ) - tB * tB * ( tF2 + tB * tB * tF4 ) ;
                else           local += sF2 + sF4 ;
                integral += 0.25e+00 * d * local ;
            }
        }
    }
    return integral ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the integral of the full spline.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real CubicSpline_IntegrateFull ( const CubicSpline *self   ,
                                 const Integer      spline ,
                                       Status      *status )
{
    Real integral = 0.0e+00 ;
    if ( self != NULL )
    {
        integral = CubicSpline_Integrate ( self, spline, Array1D_Item ( self->x, 0 ), Array1D_Item ( self->x, View1D_Extent ( self->x ) - 1 ), status ) ;
    }
    return integral ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a spline with given boundary conditions.
! . Periodic splines requires solving a cyclic tridiagonal system (to be added if necessary).
!---------------------------------------------------------------------------------------------------------------------------------*/
void CubicSpline_SetUpSpline (       CubicSpline *self            ,
                               const Integer      spline          ,
                               const Integer      lowerDerivative ,
                               const Real         lowerValue      ,
                               const Integer      upperDerivative ,
                               const Real         upperValue      ,
                                     Status      *status          )
{
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        /* . Check arguments. */
        if ( ( lowerDerivative < 1 ) ||
             ( lowerDerivative > 2 ) ||
             ( upperDerivative < 1 ) ||
             ( upperDerivative > 2 ) ||
             ( spline          < 0 ) ||
             ( spline          >= View2D_Columns ( self->y ) ) ) Status_Set ( status, Status_InvalidArgument ) ;
        /* . OK. */
        else
        {
            auto Integer      n = View1D_Extent ( self->x ) ;
            auto RealArray1D *diagonal = NULL, *rhs = NULL, *subDiagonal = NULL, *superDiagonal = NULL ;
            /* . Allocate space. */
            diagonal      = RealArray1D_AllocateWithExtent ( n    , status ) ;
            rhs           = RealArray1D_AllocateWithExtent ( n    , status ) ;
            subDiagonal   = RealArray1D_AllocateWithExtent ( n - 1, status ) ;
            superDiagonal = RealArray1D_AllocateWithExtent ( n - 1, status ) ;
            if ( Status_IsOK ( status ) )
            {
                auto Integer i ;
                auto Real    dl, du ;
                /* . Set up the tridiagonal system. */
                RealArray1D_Set ( rhs, 0.0e+00 ) ;
                /* . Lower boundary. */
                if ( lowerDerivative == 1 )
                {
                    dl = ( Array1D_Item ( self->x, 1 ) - Array1D_Item ( self->x, 0 ) ) ;
                    Array1D_Item ( diagonal     , 0 ) = dl / 3.0e+00 ;
                    Array1D_Item ( superDiagonal, 0 ) = dl / 6.0e+00 ;
                    Array1D_Item ( rhs          , 0 ) = ( Array2D_Item ( self->y, 1, spline ) - Array2D_Item ( self->y, 0, spline ) ) / dl - lowerValue ;
                }
                else if ( lowerDerivative == 2 )
                {
                    Array1D_Item ( diagonal     , 0 ) = 1.0e+00    ;
                    Array1D_Item ( superDiagonal, 0 ) = 0.0e+00    ;
                    Array1D_Item ( rhs          , 0 ) = lowerValue ;
                }
                /* . Interior conditions. */
                for ( i = 1 ; i < n - 1 ; i++ )
                {
                    dl = ( Array1D_Item ( self->x, i   ) - Array1D_Item ( self->x, i-1 ) ) ;
                    du = ( Array1D_Item ( self->x, i+1 ) - Array1D_Item ( self->x, i   ) ) ;
                    Array1D_Item (   subDiagonal, i-1 ) =   dl        / 6.0e+00 ;
                    Array1D_Item (      diagonal, i   ) = ( dl + du ) / 3.0e+00 ;
                    Array1D_Item ( superDiagonal, i   ) =        du   / 6.0e+00 ;
                    Array1D_Item ( rhs          , i   ) = ( Array2D_Item ( self->y, i+1, spline ) - Array2D_Item ( self->y, i, spline ) ) / du +
                                                          ( Array2D_Item ( self->y, i-1, spline ) - Array2D_Item ( self->y, i, spline ) ) / dl ;
                }
                /* . Upper boundary. */
                if ( upperDerivative == 1 )
                {
                    du = ( Array1D_Item ( self->x, n-1 ) - Array1D_Item ( self->x, n-2 ) ) ;
                    Array1D_Item ( subDiagonal, n-2 ) = du / 6.0e+00 ;
                    Array1D_Item ( diagonal   , n-1 ) = du / 3.0e+00 ;
                    Array1D_Item ( rhs        , n-1 ) = upperValue - ( Array2D_Item ( self->y, n-1, spline ) - Array2D_Item ( self->y, n-2, spline ) ) / du ;
                }
                else if ( upperDerivative == 2 )
                {
                    Array1D_Item ( subDiagonal, n-2 ) = 0.0e+00    ;
                    Array1D_Item ( diagonal   , n-1 ) = 1.0e+00    ;
                    Array1D_Item ( rhs        , n-1 ) = upperValue ;
                }
                /* . Solve the linear equations for the tridiagonal matrix. */
                {
                    auto integer info   = 0 ;
                    auto integer length = n ;
                    auto integer nRHS   = 1 ;
                    dgtsv_ ( &length, &nRHS, Array1D_Data ( subDiagonal ), Array1D_Data ( diagonal ), Array1D_Data ( superDiagonal ), Array1D_Data ( rhs ), &length, &info ) ;
                    /* . Set the status. */
                    if ( info != 0 ) Status_Set ( status, Status_AlgorithmError ) ;
                }
                /* . Copy rhs to h. */
                {
                    auto RealArray1D view ;
                    RealArray2D_ColumnView ( self->h, spline, False, &view, NULL ) ;
                    RealArray1D_CopyTo     ( rhs, &view, NULL ) ;
                }
            }
            /* . Deallocate space. */
            RealArray1D_Deallocate ( &diagonal      ) ;
            RealArray1D_Deallocate ( &rhs           ) ;
            RealArray1D_Deallocate ( &subDiagonal   ) ;
            RealArray1D_Deallocate ( &superDiagonal ) ;
        }
    }
}
