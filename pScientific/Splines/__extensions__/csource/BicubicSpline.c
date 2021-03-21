/*==================================================================================================================================
! . Bicubic splines.
! . Abscissae should be in strictly ascending order.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "BicubicSpline.h"
# include "Memory.h"
# include "NumericalMacros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void    BicubicSpline_EvaluateCoefficientTable     ( BicubicSpline *self, const RealArray2D *p, const RealArray2D *q, const RealArray2D *r ) ;
static void    BicubicSpline_Get1DDerivatives             ( const RealArray1D *x, const RealArray1D *u, RealArray1D *d, BicubicSplineType type, RealArray1D *Ad, RealArray1D *Asd, RealArray1D *qdy, RealArray1D *lll ) ;
static Integer BicubicSpline_Locate                       ( const RealArray1D *abscissa, const Real x, const Boolean isPeriodic ) ;
static void    BicubicSpline_Setup                        ( BicubicSpline *self, Status *status ) ;
static void    BicubicSpline_TridiagonalLDLtSolve         ( const Integer n, RealArray1D *d, RealArray1D *l, RealArray1D *b ) ;
static void    BicubicSpline_TridiagonalLDLtSolvePeriodic ( const Integer n, RealArray1D *d, RealArray1D *lsd, RealArray1D *lll, RealArray1D *b ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
BicubicSpline *BicubicSpline_Allocate ( const Integer lengthX, const Integer lengthY, const Boolean doX, const Boolean doY, const Boolean doF, Status *status )
{
    BicubicSpline *self = NULL ;
    if ( ( lengthX >= 2 ) && ( lengthY >= 2 ) )
    {
        self = Memory_AllocateType ( BicubicSpline ) ;
        if ( self != NULL )
        {
            auto Integer lengths[4] ;
            auto Status  localStatus = Status_OK ;
            self->type         = BicubicSplineType_Natural ;
            self->lengthX      = lengthX ;
            self->lengthY      = lengthY ;
            self->x            = NULL    ;
            self->y            = NULL    ;
            self->f            = NULL    ;
            self->coefficients = NULL    ;
            /* . Array allocation. */
            lengths[0] = lengthX - 1 ; lengths[1] = lengthY - 1 ; lengths[2] = 4 ; lengths[3] = 4 ;
            self->coefficients = RealArrayND_AllocateWithShape   ( 4, lengths, &localStatus ) ;
/*
{
auto Integer i ;
auto ViewND *view = self->coefficients->view ;
printf ( "Array Data %d %d %d\n", view->rank, view->offset, view->size ) ;
for ( i = 0 ; i < view->rank ; i++ )
{
    printf ( "%d %d\n", view->extents[i], view->strides[i] ) ;
}
fflush ( stdout ) ;
}
*/
            if ( doX ) self->x = RealArray1D_AllocateWithExtent  ( lengthX, &localStatus ) ;
            if ( doY ) self->y = RealArray1D_AllocateWithExtent  ( lengthY, &localStatus ) ;
            if ( doF ) self->f = RealArray2D_AllocateWithExtents ( lengthX, lengthY, &localStatus ) ;
            if ( ! Status_IsValueOK ( localStatus ) ) BicubicSpline_Deallocate ( &self ) ;
        }
        if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    }
    else Status_Set ( status, Status_InvalidArgument ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
BicubicSpline *BicubicSpline_Clone ( const BicubicSpline *self, Status *status )
{
    BicubicSpline *clone = NULL ;
    if ( self != NULL )
    {
        clone = BicubicSpline_Allocate ( self->lengthX, self->lengthY, True, True, True, status ) ;
        if ( clone != NULL )
        {
            clone->type = self->type ;
            RealArray1D_CopyTo ( self->x, clone->x, status ) ;
            RealArray1D_CopyTo ( self->y, clone->y, status ) ;
            RealArray2D_CopyTo ( self->f, clone->f, status ) ;
            RealArrayND_CopyTo ( self->coefficients, clone->coefficients, status ) ;
        }
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void BicubicSpline_Deallocate ( BicubicSpline **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        RealArray1D_Deallocate ( &((*self)->x) ) ;
        RealArray1D_Deallocate ( &((*self)->y) ) ;
        RealArray2D_Deallocate ( &((*self)->f) ) ;
        RealArrayND_Deallocate ( &((*self)->coefficients) ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Evaluation (function and first derivatives only).
!---------------------------------------------------------------------------------------------------------------------------------*/
void BicubicSpline_Evaluate ( const BicubicSpline *self, const Real x, const Real y, Real *f, Real *g1, Real *g2, Status *status )
{
    if ( f  != NULL ) (*f)  = 0.0e+00 ;
    if ( g1 != NULL ) (*g1) = 0.0e+00 ;
    if ( g2 != NULL ) (*g2) = 0.0e+00 ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer ix, iy ;
        ix = BicubicSpline_Locate ( self->x, x, ( self->type == BicubicSplineType_Periodic ) ) ;
        iy = BicubicSpline_Locate ( self->y, y, ( self->type == BicubicSplineType_Periodic ) ) ;
        if ( ( ix >= 0 ) && ( iy >= 0 ) )
        {
            auto Integer     i, indices[2] ;
            auto Real        dudx, dudy, dx, dy, u ;
            auto RealArray2D C ;
            indices[0] = ix ; indices[1] = iy ;
            RealArrayND_ViewTail2D ( self->coefficients, indices, False, &C, status ) ;
            dx   = x - Array1D_Item ( self->x, ix ) ;
            dy   = y - Array1D_Item ( self->y, iy ) ;
            u    = 0.0e+00 ;
            dudx = 0.0e+00 ;
            dudy = 0.0e+00 ;
            for ( i = 3 ; i >= 0 ; i-- )
            {
                u    = Array2D_Item ( &C, i, 0 ) + dy * ( Array2D_Item ( &C, i, 1 ) + dy * ( Array2D_Item ( &C, i, 2 ) + dy * Array2D_Item ( &C, i, 3 ) ) ) + u * dx ;
                dudx = Array2D_Item ( &C, 1, i ) + dx * ( 2.0e+00 * Array2D_Item ( &C, 2, i ) + 3.0e+00 * dx * Array2D_Item ( &C, 3, i ) ) + dudx * dy ;
                dudy = Array2D_Item ( &C, i, 1 ) + dy * ( 2.0e+00 * Array2D_Item ( &C, i, 2 ) + 3.0e+00 * dy * Array2D_Item ( &C, i, 3 ) ) + dudy * dx ;
            }
            if ( f  != NULL ) (*f)  = u ;
            if ( g1 != NULL ) (*g1) = dudx ;
            if ( g2 != NULL ) (*g2) = dudy ;
        }
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Evaluate the coefficient table.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void BicubicSpline_EvaluateCoefficientTable ( BicubicSpline *self, const RealArray2D *p, const RealArray2D *q, const RealArray2D *r )
{
    auto Integer     i, j, indices[2] ;
    auto Real        a, b, c, d, dx, dy ;
    auto RealArray2D C, *u = self->f ;
    for ( j = 0 ; j < self->lengthY - 1 ; j++ )
    {
        dy = 1.0e+00 / ( Array1D_Item ( self->y, j+1 ) - Array1D_Item ( self->y, j ) ) ;
        for ( i = 0 ; i < self->lengthX - 1 ; i++ )
        {
            dx = 1.0e+00 / ( Array1D_Item ( self->x, i+1 ) - Array1D_Item ( self->x, i ) ) ;
            indices[0] = i ; indices[1] = j ;
            RealArrayND_ViewTail2D ( self->coefficients, indices, False, &C, NULL ) ;

            Array2D_Item ( &C, 0, 0 ) = Array2D_Item ( u, i, j ) ;
            Array2D_Item ( &C, 1, 0 ) = Array2D_Item ( p, i, j ) ;
            Array2D_Item ( &C, 0, 1 ) = Array2D_Item ( q, i, j ) ;
            Array2D_Item ( &C, 1, 1 ) = Array2D_Item ( r, i, j ) ;

            a = ( Array2D_Item ( u, i+1, j ) - Array2D_Item ( u, i, j ) ) * dx ;
            Array2D_Item ( &C, 2, 0 ) = ( 3.0e+00 * a - 2.0e+00 * Array2D_Item ( p, i, j ) - Array2D_Item ( p, i+1, j ) ) * dx ;
            Array2D_Item ( &C, 3, 0 ) = ( Array2D_Item ( p, i+1, j ) + Array2D_Item ( p, i, j ) - 2.0e+00 * a ) * ( dx * dx ) ;

            a = ( Array2D_Item ( u, i, j+1 ) - Array2D_Item ( u, i, j ) ) * dy ;
            Array2D_Item ( &C, 0, 2 ) = ( 3.0e+00 * a - 2.0e+00 * Array2D_Item ( q, i, j ) - Array2D_Item ( q, i, j+1 ) ) * dy ;
            Array2D_Item ( &C, 0, 3 ) = ( Array2D_Item ( q, i, j+1 ) + Array2D_Item ( q, i, j ) - 2.0e+00 * a ) * ( dy * dy ) ;

            a = ( Array2D_Item ( q, i+1, j ) - Array2D_Item ( q, i, j ) ) * dx ;
            Array2D_Item ( &C, 2, 1 ) = ( 3.0e+00 * a - Array2D_Item ( r, i+1, j ) - 2.0e+00 * Array2D_Item ( r, i, j ) ) * dx ;
            Array2D_Item ( &C, 3, 1 ) = ( Array2D_Item ( r, i+1, j ) + Array2D_Item ( r, i, j ) - 2.0e+00 * a ) * ( dx * dx ) ;

            a = ( Array2D_Item ( p, i, j+1 ) - Array2D_Item ( p, i, j ) ) * dy ;
            Array2D_Item ( &C, 1, 2 ) = ( 3.0e+00 * a - Array2D_Item ( r, i, j+1 ) - 2.0e+00 * Array2D_Item ( r, i, j ) ) * dy ;
            Array2D_Item ( &C, 1, 3 ) = ( Array2D_Item ( r, i, j+1 ) + Array2D_Item ( r, i, j ) - 2.0e+00 * a ) * ( dy * dy ) ;

            a = ( Array2D_Item ( u, i+1, j+1 ) + Array2D_Item ( u, i, j ) - Array2D_Item ( u, i+1, j ) - Array2D_Item ( u, i, j+1 ) ) * ( dx * dx * dy * dy ) -
                ( Array2D_Item ( p, i, j+1 ) - Array2D_Item ( p, i, j ) ) * ( dx * dy * dy ) - ( Array2D_Item ( q, i+1, j ) - Array2D_Item ( q, i, j ) ) * ( dx * dx * dy ) + Array2D_Item ( r, i, j ) * ( dx * dy ) ;
            b = ( Array2D_Item ( p, i+1, j+1 ) + Array2D_Item ( p, i, j ) - Array2D_Item ( p, i+1, j ) - Array2D_Item ( p, i, j+1 ) ) * ( dx * dy * dy ) - ( Array2D_Item ( r, i+1, j ) - Array2D_Item ( r, i, j ) ) * ( dx * dy ) ;
            c = ( Array2D_Item ( q, i+1, j+1 ) + Array2D_Item ( q, i, j ) - Array2D_Item ( q, i+1, j ) - Array2D_Item ( q, i, j+1 ) ) * ( dx * dx * dy ) - ( Array2D_Item ( r, i, j+1 ) - Array2D_Item ( r, i, j ) ) * ( dx * dy ) ;
            d = ( Array2D_Item ( r, i+1, j+1 ) + Array2D_Item ( r, i, j ) - Array2D_Item ( r, i+1, j ) - Array2D_Item ( r, i, j+1 ) ) * ( dx * dy ) ;
            Array2D_Item ( &C, 2, 2 ) =    9.0e+00 * a - 3.0e+00 * b - 3.0e+00 * c + d ;
            Array2D_Item ( &C, 2, 3 ) = ( -6.0e+00 * a + 2.0e+00 * b + 3.0e+00 * c - d ) * dy ;
            Array2D_Item ( &C, 3, 2 ) = ( -6.0e+00 * a + 3.0e+00 * b + 2.0e+00 * c - d ) * dx ;
            Array2D_Item ( &C, 3, 3 ) = (  4.0e+00 * a - 2.0e+00 * b - 2.0e+00 * c + d ) * dx * dy ;
/*
printf ( "\nCoefficients for (%d,%d):\n", i, j ) ;
RealArray2D_Print ( &C ) ;
*/
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the 1-D derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void BicubicSpline_Get1DDerivatives ( const RealArray1D *x, const RealArray1D *y, RealArray1D *d, BicubicSplineType type, RealArray1D *Ad, RealArray1D *Asd, RealArray1D *qdy, RealArray1D *lll )
{
    Integer     i, n ;
    Real        r ;
    RealArray1D a, b, c ;

    /* . Initialization. */
    n = View1D_Extent ( x ) ;

    /* . Setup. */
    for ( i = 0 ; i < n-1 ; i++ )
    {
        Array1D_Item ( Asd, i ) = 1.0e+00 / ( Array1D_Item ( x, i+1 ) - Array1D_Item ( x, i ) ) ;
        Array1D_Item ( qdy, i ) = ( Array1D_Item ( y, i+1 ) - Array1D_Item ( y, i ) ) * pow ( Array1D_Item ( Asd, i ), 2 ) ;
    }
    for ( i = 1 ; i < n-1 ; i++ )
    {
        Array1D_Item ( Ad, i ) = 2.0e+00 * ( Array1D_Item ( Asd, i-1 ) + Array1D_Item ( Asd, i ) ) ;
        Array1D_Item (  d, i ) = 3.0e+00 * ( Array1D_Item ( qdy, i-1 ) + Array1D_Item ( qdy, i ) ) ;
    }

    /* . Branch on the type. */
    switch ( type )
    {
       case BicubicSplineType_Clamped:
           Array1D_Item ( d, 1   ) -= Array1D_Item ( d, 0   ) * Array1D_Item ( Asd, 0   ) ;
           Array1D_Item ( d, n-2 ) -= Array1D_Item ( d, n-1 ) * Array1D_Item ( Asd, n-2 ) ;
           RealArray1D_View ( Ad ,  1, n-2, 1, False, &a, NULL ) ;
           RealArray1D_View ( Asd,  1, n-2, 1, False, &b, NULL ) ;
           RealArray1D_View ( d  ,  1, n-2, 1, False, &c, NULL ) ;
           BicubicSpline_TridiagonalLDLtSolve ( n-2, &a, &b, &c ) ;
           break ;
       case BicubicSplineType_Natural:
           Array1D_Item ( Ad, 0   ) = 2.0e+00 * Array1D_Item ( Asd, 0   ) ;
           Array1D_Item (  d, 0   ) = 3.0e+00 * Array1D_Item ( qdy, 0   ) ;
           Array1D_Item ( Ad, n-1 ) = 2.0e+00 * Array1D_Item ( Asd, n-2 ) ;
           Array1D_Item (  d, n-1 ) = 3.0e+00 * Array1D_Item ( qdy, n-2 ) ;
           BicubicSpline_TridiagonalLDLtSolve ( n, Ad, Asd, d ) ;
           break ;
       case BicubicSplineType_NotAKnot:
           r = Array1D_Item ( Asd, 1   ) / Array1D_Item ( Asd, 0   ) ;
           Array1D_Item ( Ad, 0   ) = Array1D_Item ( Asd, 0   ) / ( 1.0e+00 + r ) ;
           Array1D_Item ( d,  0   ) = ( ( 3.0e+00 * r + 2.0e+00 ) * Array1D_Item ( qdy, 0   ) + r * Array1D_Item ( qdy, 1   ) ) / pow ( ( 1.0e+00 + r ), 2 ) ;
           r = Array1D_Item ( Asd, n-3 ) / Array1D_Item ( Asd, n-2 ) ;
           Array1D_Item ( Ad, n-1 ) = Array1D_Item ( Asd, n-2 ) / ( 1.0e+00 + r ) ;
           Array1D_Item ( d,  n-1 ) = ( ( 3.0e+00 * r + 2.0e+00 ) * Array1D_Item ( qdy, n-2 ) + r * Array1D_Item ( qdy, n-3 ) ) / pow ( ( 1.0e+00 + r ), 2 ) ;
           BicubicSpline_TridiagonalLDLtSolve ( n, Ad, Asd, d ) ;
           break ;
       case BicubicSplineType_Periodic:
           RealArray1D_Set ( lll, 0.0e+00 ) ;
           Array1D_Item ( Ad,  0   ) = 2.0e+00 * ( Array1D_Item ( Asd, 0 ) + Array1D_Item ( Asd, n-2 ) ) ;
           Array1D_Item ( d,   0   ) = 3.0e+00 * ( Array1D_Item ( qdy, 0 ) + Array1D_Item ( qdy, n-2 ) ) ;
           Array1D_Item ( lll, 0   ) = Array1D_Item ( Asd, n-2 ) ;
           Array1D_Item ( lll, n-3 ) = Array1D_Item ( Asd, n-3 ) ;
           BicubicSpline_TridiagonalLDLtSolvePeriodic ( n-1, Ad, Asd, lll, d ) ;
           Array1D_Item ( d, n-1 ) = Array1D_Item ( d, 0 ) ;
           break ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Locate the index of a point such that point lies between i and i+1.
! . Points outside the range return -1.
! . This should never happen for periodic abscissa.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Integer BicubicSpline_Locate ( const RealArray1D *abscissa, const Real x, const Boolean isPeriodic )
{
    Integer index = -1 ;
    if ( abscissa != NULL )
    {
        auto Integer n ;
        auto Real    a, al, au ;

        /* . Initialization. */
        n  = View1D_Extent ( abscissa ) ;
        a  = x ;
        al = Array1D_Item ( abscissa, 0     ) ;
        au = Array1D_Item ( abscissa, n - 1 ) ;

        /* . Adjust x for a periodic abscissa. */
        if ( ( isPeriodic ) && ( ( a < al ) || ( a > au ) ) )
        {
            auto Real dx, rF, rI ;
            dx = au - al ;
            rF = modf ( ( a - al ) / dx, &rI );
            if ( rF >= 0.0e+00 ) a = al + rF * dx ;
            else                 a = au + rF * dx ;
            if      ( a < al ) a = al ;
            else if ( a > au ) a = au ;
        }

        /* . Locate the index. */
        if ( ( a >= al ) && ( a <= au ) )
        {
            auto Integer i, l, u ;
            l = 0     ;
            u = n - 1 ;
            while ( ( u - l ) > 1 )
            {
                i = ( u + l ) >> 1 ;
                if ( Array1D_Item ( abscissa, i ) >= a ) u = i ;
                else                                     l = i ;
            }
            index = l ;
        }
    }
    return index ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a spline given (x,y,f) and a type.
!---------------------------------------------------------------------------------------------------------------------------------*/
BicubicSpline *BicubicSpline_MakeFromRealArray2D ( RealArray1D **x, RealArray1D **y, RealArray2D **f, const BicubicSplineType type, Status *status )
{
    BicubicSpline *self = NULL ;
    Boolean        isOK ;
    isOK = ( x != NULL ) && ( (*x) != NULL ) && ( y != NULL ) && ( (*y) != NULL ) && ( f != NULL ) && ( (*f) != NULL ) ;
    if ( isOK )
    {
        auto Integer lengthX, lengthY ;

        /* . Checks. */
        lengthX = View1D_Extent ( (*x) ) ;
        lengthY = View1D_Extent ( (*y) ) ;
        isOK = ( lengthX > 1 ) && ( lengthY > 1 ) && ( View2D_Rows ( (*f) ) == lengthX ) && ( View2D_Columns ( (*f) ) == lengthY ) ;
        if ( isOK )
        {
            auto Integer i ;
            for ( i = 1 ; i < lengthX ; i++ ) { if ( Array1D_Item ( (*x), i ) <=  Array1D_Item ( (*x), i-1 ) ) { isOK = False ; break ; } }
            for ( i = 1 ; i < lengthY ; i++ ) { if ( Array1D_Item ( (*y), i ) <=  Array1D_Item ( (*y), i-1 ) ) { isOK = False ; break ; } }
        }

        /* . OK so far. */
        if ( isOK )
        {
            /* . Allocate space. */
            self = BicubicSpline_Allocate ( lengthX, lengthY, False, False, False, status ) ;
            if ( self != NULL )
            {
                auto Status localStatus = Status_OK ;
                /* . Finish construction. */
                self->type = type ;
                self->x = (*x) ; (*x) = NULL ;
                self->y = (*y) ; (*y) = NULL ;
                self->f = (*f) ; (*f) = NULL ;
                BicubicSpline_Setup ( self, &localStatus ) ;
                if ( ! Status_IsValueOK ( localStatus ) )
                {
                    BicubicSpline_Deallocate ( &self ) ;
                    Status_Set ( status, localStatus ) ;
                }
            }
            else Status_Set ( status, Status_OutOfMemory ) ;
        }
    }
    if ( ! isOK ) Status_Set ( status, Status_InvalidArgument ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Setup the bicubic spline.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void BicubicSpline_Setup ( BicubicSpline *self, Status *status )
{
    Integer      n ;
    RealArray1D *Ad, *Asd, *d, *ll = NULL, *qdu, t, u ;
    RealArray2D *p, *q, *r ;
    Status       localStatus = Status_OK ;

    /* . Initialization . */
    n = Maximum ( self->lengthX, self->lengthY ) ;

    /* . Allocate space. */
    Ad  = RealArray1D_AllocateWithExtent  ( n, &localStatus ) ;
    Asd = RealArray1D_AllocateWithExtent  ( n, &localStatus ) ;
    d   = RealArray1D_AllocateWithExtent  ( self->lengthY, &localStatus ) ;
    if ( self->type == BicubicSplineType_Periodic ) ll = RealArray1D_AllocateWithExtent ( n, &localStatus ) ;
    qdu = RealArray1D_AllocateWithExtent  ( n, &localStatus ) ;
    p   = RealArray2D_AllocateWithExtents ( self->lengthX, self->lengthY, &localStatus ) ;
    q   = RealArray2D_AllocateWithExtents ( self->lengthX, self->lengthY, &localStatus ) ;
    r   = RealArray2D_AllocateWithExtents ( self->lengthX, self->lengthY, &localStatus ) ;
    if ( Status_IsValueOK ( localStatus ) )
    {
        auto Integer i ;

        /* . du/dx. */
        for ( i = 0 ; i < self->lengthY ; i++ )
        {
            RealArray2D_ColumnView ( p      , i, False, &t, &localStatus ) ;
            RealArray2D_ColumnView ( self->f, i, False, &u, &localStatus ) ;
            BicubicSpline_Get1DDerivatives ( self->x, &u, &t, self->type, Ad, Asd, qdu, ll ) ;
        }
        /* . du/dy. */
        for ( i = 0 ; i < self->lengthX ; i++ )
        {
            RealArray2D_RowView ( q      , i, False, &t, &localStatus ) ;
            RealArray2D_RowView ( self->f, i, False, &u, &localStatus ) ;
            BicubicSpline_Get1DDerivatives ( self->y, &u, d, self->type, Ad, Asd, qdu, ll ) ;
            RealArray1D_CopyTo ( d, &t, &localStatus ) ;
        }
        /* . d2u/dxdy. */
        RealArray2D_ColumnView ( q, 0, False, &u, &localStatus ) ;
        RealArray2D_ColumnView ( r, 0, False, &t, &localStatus ) ;
        BicubicSpline_Get1DDerivatives ( self->x, &u, &t, self->type, Ad, Asd, qdu, ll ) ;
        RealArray2D_ColumnView ( q, self->lengthY-1, False, &u, &localStatus ) ;
        RealArray2D_ColumnView ( r, self->lengthY-1, False, &t, &localStatus ) ;
        BicubicSpline_Get1DDerivatives ( self->x, &u, &t, self->type, Ad, Asd, qdu, ll ) ;
        for ( i = 0 ; i < self->lengthX ; i++ )
        {
            RealArray2D_RowView ( p, i, False, &u, &localStatus ) ;
            Array1D_Item ( d, 0               ) = Array2D_Item ( r, i, 0               ) ;
            Array1D_Item ( d, self->lengthY-1 ) = Array2D_Item ( r, i, self->lengthY-1 ) ;
            BicubicSpline_Get1DDerivatives ( self->y, &u, d, BicubicSplineType_Clamped, Ad, Asd, qdu, ll ) ;
            RealArray1D_View   ( d,       1, self->lengthY-2, 1, False, &u, &localStatus ) ;
            RealArray2D_View1D ( r, 1, i, 1, self->lengthY-2, 1, False, &t, &localStatus ) ;
            RealArray1D_CopyTo ( &u, &t, &localStatus ) ;
        }

        /* . Coefficient table. */
        BicubicSpline_EvaluateCoefficientTable ( self, p, q, r ) ;
    }
    else Status_Set ( status, localStatus ) ;

    /* . Finish up. */
    RealArray1D_Deallocate ( &Ad  ) ;
    RealArray1D_Deallocate ( &Asd ) ;
    RealArray1D_Deallocate ( &d   ) ;
    RealArray1D_Deallocate ( &ll  ) ;
    RealArray1D_Deallocate ( &qdu ) ;
    RealArray2D_Deallocate ( &p   ) ;
    RealArray2D_Deallocate ( &q   ) ;
    RealArray2D_Deallocate ( &r   ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Solution of a tridiagonal system.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void BicubicSpline_TridiagonalLDLtSolve ( const Integer n, RealArray1D *d, RealArray1D *l, RealArray1D *b )
{
    auto Integer i ;
    auto Real temp ;
    for ( i = 1 ; i < n ; i++ )
    {
        temp = Array1D_Item ( l, i-1 ) ;
        Array1D_Item ( l, i-1 ) /= Array1D_Item ( d, i-1 ) ;
        Array1D_Item ( d, i   ) -= temp * Array1D_Item ( l, i-1 ) ;
        Array1D_Item ( b, i   ) -= Array1D_Item ( l, i-1 ) * Array1D_Item ( b, i-1 ) ;
    }
    Array1D_Item ( b, n-1 ) /= Array1D_Item ( d, n-1 ) ;
    for ( i = n-2 ; i >= 0 ; i-- ) Array1D_Item ( b, i ) = ( Array1D_Item ( b, i ) / Array1D_Item ( d, i ) - Array1D_Item ( l, i ) * Array1D_Item ( b, i+1 ) ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Solution of a periodic tridiagonal system.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void BicubicSpline_TridiagonalLDLtSolvePeriodic ( const Integer n, RealArray1D *d, RealArray1D *lsd, RealArray1D *lll, RealArray1D *b )
{
    auto Integer i ;
    auto Real temp1, temp2 ;
    for ( i = 0 ; i < n-2 ; i++ )
    {
       temp1 = Array1D_Item ( lsd, i ) ;
       temp2 = Array1D_Item ( lll, i ) ;
       Array1D_Item ( lsd, i   ) /= Array1D_Item ( d, i ) ;
       Array1D_Item ( lll, i   ) /= Array1D_Item ( d, i ) ;
       Array1D_Item (   d, i+1 ) -= Array1D_Item ( lsd, i ) * temp1 ;
       Array1D_Item ( lll, i+1 ) -= Array1D_Item ( lll, i ) * temp1 ;
       Array1D_Item (   d, n-1 ) -= Array1D_Item ( lll, i ) * temp2 ;
    }
    temp2 = Array1D_Item ( lll, n-2 ) ;
    Array1D_Item ( lll, n-2 ) /= Array1D_Item (   d, n-2 ) ;
    Array1D_Item (   d, n-1 ) -= Array1D_Item ( lll, n-2 ) * temp2 ;
    for ( i = 1 ; i < n-1 ; i++ ) Array1D_Item ( b, i   ) -= Array1D_Item ( lsd, i-1 ) * Array1D_Item ( b, i-1 ) ;
    for ( i = 0 ; i < n-1 ; i++ ) Array1D_Item ( b, n-1 ) -= Array1D_Item ( lll, i   ) * Array1D_Item ( b, i   ) ;
    for ( i = 0 ; i < n   ; i++ ) Array1D_Item ( b, i   ) /= Array1D_Item (   d, i   ) ;
    Array1D_Item ( b, n-2 ) -= Array1D_Item ( lll, n-2 ) * Array1D_Item ( b, n-1 ) ;
    for ( i = n-3 ; i >= 0 ; i-- ) Array1D_Item ( b, i ) -= ( Array1D_Item ( lsd, i ) * Array1D_Item ( b, i+1 ) + Array1D_Item ( lll, i ) * Array1D_Item ( b, n-1 ) ) ;
}

