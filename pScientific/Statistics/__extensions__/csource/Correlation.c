/*==================================================================================================================================
! . Functions for calculating correlation functions.
!=================================================================================================================================*/

/* . For cross-correlation the function calculated is X(t)Y(0) so need a second call to get X(0)Y(t). */

# include <math.h>
# include <stdio.h>

# include "Correlation.h"
# include "NumericalMacros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define DEFAULTSMALL 1.0e-10

/*----------------------------------------------------------------------------------------------------------------------------------
! . Private procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void Correlation_DotProductDirect      ( const RealArray2D *x, const RealArray2D *y, RealArray1D *c, Status *status ) ;
static void Correlation_DotProductNormalize   ( const RealArray2D *x, const RealArray2D *y, const Real small, RealArray1D *c, Status *status ) ;
static void Correlation_DotProductRemoveMeans ( RealArray2D *x, Status *status ) ;

static void Correlation_SimpleDirect          ( const RealArray1D *x, const RealArray1D *y, RealArray1D *c ) ;
static void Correlation_SimpleNormalize       ( const RealArray1D *x, const RealArray1D *y, const Real small, RealArray1D *c ) ;
static void Correlation_SimpleRemoveMean      ( RealArray1D *x ) ;

/*==================================================================================================================================
! . Public procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Dot-product auto or cross-correlation function.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Correlation_MakeDotProduct (       RealArray2D *x            ,
                                        RealArray2D *y            ,
                                  const Boolean      useFFT       ,
                                  const Boolean      normalize    ,
                                  const Boolean      removeMean   ,
                                  const Integer      tCorrelation ,
                                  const Real        *tolerance    ,
                                        RealArray1D *f            ,
                                        Status      *status       )
{
    if ( ( x != NULL ) && ( f != NULL ) && Status_IsOK ( status ) )
    {
        if ( View1D_Extent ( f ) >= tCorrelation )
        {
            /* . Remove means. */
            if ( removeMean )
            {
                Correlation_DotProductRemoveMeans ( x, status ) ;
                if ( y != x ) Correlation_DotProductRemoveMeans ( y, status ) ;
            }
            /* . Calculate the correlation function. */
            Correlation_DotProductDirect ( x, y, f, status ) ;
/*
            if ( useFFT ) Correlation_DotProductFFT    ( x, y, f, status ) ;
            else          Correlation_DotProductDirect ( x, y, f, status ) ;
*/
            /* . Normalize. */
            if ( normalize )
            {
                auto Real small ;
                if ( tolerance == NULL ) small = DEFAULTSMALL ;
                else                     small = fabs ( *tolerance ) ;
                Correlation_DotProductNormalize ( x, y, small, f, status ) ;
            }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Simple auto or cross-correlation function.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Correlation_MakeSimple (       RealArray1D *x            ,
                                    RealArray1D *y            ,
                              const Boolean      useFFT       ,
                              const Boolean      normalize    ,
                              const Boolean      removeMean   ,
                              const Integer      tCorrelation ,
                              const Real        *tolerance    ,
                                    RealArray1D *f            ,
                                    Status      *status       )
{
    if ( ( x != NULL ) && ( f != NULL ) && Status_IsOK ( status ) )
    {
        if ( View1D_Extent ( f ) >= tCorrelation )
        {
            /* . Remove means. */
            if ( removeMean )
            {
                Correlation_SimpleRemoveMean ( x ) ;
                if ( y != x ) Correlation_SimpleRemoveMean ( y ) ;
            }
            /* . Calculate the correlation function. */
            Correlation_SimpleDirect ( x, y, f ) ;
/*
            if ( useFFT ) Correlation_SimpleFFT    ( x, y, f, status ) ;
            else          Correlation_SimpleDirect ( x, y, f ) ;
*/
            /* . Normalize. */
            if ( normalize )
            {
                auto Real small ;
                if ( tolerance == NULL ) small = DEFAULTSMALL ;
                else                     small = fabs ( *tolerance ) ;
                Correlation_SimpleNormalize ( x, y, small, f ) ;
            }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*==================================================================================================================================
! . Private procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Dot-product auto or cross-correlation function via a direct method.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void Correlation_DotProductDirect ( const RealArray2D *x, const RealArray2D *y, RealArray1D *c, Status *status )
{
    if ( ( c != NULL ) && ( x != NULL ) )
    {
        auto       Integer      i, t, tCorrelation, tMaximum, tRun ;
        auto       Real         sum   ;
        auto       RealArray1D  aRow, bRow ;
        auto const RealArray2D *a, *b ;
        /* . Assign array aliases. */
        a = x ;
        if ( y == NULL ) b = x ;
        else             b = y ;
        /* . Get problem size. */
        tRun         = Minimum ( View2D_Rows ( a ), View2D_Rows ( b ) ) ;
        tCorrelation = Minimum ( View1D_Extent ( c ) - 1, tRun - 1 ) ;
        /* . Initialization. */
        RealArray1D_Set ( c, 0.0e+00 ) ;
        /* . Loop over the elements in the correlation function. */
        for ( t = 0 ; t <= tCorrelation ; t++ )
        {
            tMaximum = tRun - t ;
            sum      = 0.0e+00  ;
            for ( i = 0 ; i < tMaximum ; i++ )
            {
                RealArray2D_RowView ( a, t+i, False, &aRow, status ) ;
                RealArray2D_RowView ( b,   i, False, &bRow, status ) ;
                sum += RealArray1D_Dot ( &aRow, &bRow, status ) ;
            }
            Array1D_Item ( c, t ) = sum / ( ( Real ) tMaximum ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Normalize a dot-product correlation function.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void Correlation_DotProductNormalize ( const RealArray2D *x, const RealArray2D *y, const Real small, RealArray1D *c, Status *status )
{
    if ( ( x != NULL ) && ( c != NULL ) )
    {
        auto Real scale ;
        /* . Auto-correlation. */
        if ( y == NULL ) scale = Array1D_Item ( c, 0 ) ;
        /* . Cross-correlation. */
        else
        {
            auto Integer     i ;
            auto Real        x2, y2 ;
            auto RealArray1D column ;
            x2 = 0.0e+00 ;
            y2 = 0.0e+00 ;
            for ( i = 0 ; i < View2D_Columns ( x ) ; i++ )
            {
                RealArray2D_ColumnView ( x, i, False, &column, status ) ;
                x2 += RealArray1D_Dot  ( &column, &column, status ) ;
            }
            for ( i = 0 ; i < View2D_Columns ( y ); i++ )
            {
                RealArray2D_ColumnView ( y, i, False, &column, status ) ;
                y2 += RealArray1D_Dot  ( &column, &column, status ) ;
            }
            scale = sqrt ( fabs ( x2 * y2 ) ) ;
        }
        /* . Scaling. */
        if ( fabs ( scale ) > small ) RealArray1D_Scale ( c, 1.0 / scale ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Remove the means from a 2-D array of data.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void Correlation_DotProductRemoveMeans ( RealArray2D *x, Status *status )
{
    auto Integer n ;
    n = View2D_Rows ( x ) ;
    if ( ( x != NULL ) && ( n > 0 ) )
    {
        auto Integer     i      ;
        auto Real        mean   ;
        auto RealArray1D column ;
        for ( i = 0 ; i < View2D_Columns ( x ) ; i++ )
        {
            RealArray2D_ColumnView ( x, i, False, &column, status ) ;
            mean = RealArray1D_Sum ( &column ) / ( ( Real ) n ) ;
            RealArray1D_Increment  ( &column, -mean ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Simple auto or cross-correlation function via a direct method.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void Correlation_SimpleDirect ( const RealArray1D *x, const RealArray1D *y, RealArray1D *c )
{
    if ( ( c != NULL ) && ( x != NULL ) )
    {
        auto       Integer      i, t, tCorrelation, tMaximum, tRun ;
        auto       Real         sum   ;
        auto const RealArray1D *a, *b ;
        /* . Assign array aliases. */
        a = x ;
        if ( y == NULL ) b = x ;
        else             b = y ;
        /* . Get problem size. */
        tRun         = Minimum ( View1D_Extent ( a ), View1D_Extent ( b ) ) ;
        tCorrelation = Minimum ( View1D_Extent ( c ) - 1, tRun - 1 ) ;
        /* . Initialization. */
        RealArray1D_Set ( c, 0.0e+00 ) ;
        /* . Loop over the elements in the correlation function. */
        for ( t = 0 ; t <= tCorrelation ; t++ )
        {
            tMaximum = tRun - t ;
            sum      = 0.0e+00  ;
            for ( i = 0 ; i < tMaximum ; i++ ) sum += ( Array1D_Item ( a, t+i ) * Array1D_Item ( b, i ) ) ;
            Array1D_Item ( c, t ) = sum / ( ( Real ) tMaximum ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Normalize a simple correlation function.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void Correlation_SimpleNormalize ( const RealArray1D *x, const RealArray1D *y, const Real small, RealArray1D *c )
{
    if ( ( x != NULL ) && ( c != NULL ) )
    {
        auto Real scale ;
        /* . Auto-correlation. */
        if ( y == NULL ) scale = Array1D_Item ( c, 0 ) ;
        /* . Cross-correlation. */
        else             scale = sqrt ( fabs ( RealArray1D_Dot ( x, x, NULL ) * RealArray1D_Dot ( y, y, NULL ) ) ) ;
        /* . Scaling. */
        if ( fabs ( scale ) > small ) RealArray1D_Scale ( c, 1.0 / scale ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Remove the mean from an array of data.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void Correlation_SimpleRemoveMean ( RealArray1D *x )
{
    auto Integer n ;
    n = View1D_Extent ( x ) ;
    if ( ( x != NULL ) && ( n > 0 ) )
    {
        auto Real mean ;
        mean = RealArray1D_Sum ( x ) / ( ( Real ) n ) ;
        RealArray1D_Increment ( x, -mean ) ;
    }
}
