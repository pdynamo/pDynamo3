/*==================================================================================================================================
! . Gaussian basis set normalization.
!=================================================================================================================================*/

# include <math.h>
# include <stdarg.h>
# include <stdlib.h>

# include "GaussianBasisIntegrals_b2e1n0.h"
# include "GaussianBasisNormalize.h"
# include "NumericalMacros.h"
# include "OrthogonalizingTransformation.h"
# include "SymmetricMatrix.h"

/*# define _DebugPrint*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check normalization.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real CheckNormalization ( GaussianBasis *self, RealArray2D *S, Status *status )
{
    Real maximumDeviation = 0.0e+00 ;
    if ( ( self != NULL ) && ( S != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer      i, o, w ;
        auto RealArray2D *a, *b ;
        /* . Form (c->o)^T * S * (c->o). */
        w = View2D_Rows    ( self->c2o ) ;
        o = View2D_Columns ( self->c2o ) ;
        a = RealArray2D_AllocateWithExtents ( w, o, NULL ) ;
        b = RealArray2D_AllocateWithExtents ( o, o, NULL ) ;
        RealArray2D_MatrixMultiply ( False, False, 1.0e+00, S        , self->c2o, 0.0e+00, a, NULL ) ;
        RealArray2D_MatrixMultiply ( True , False, 1.0e+00, self->c2o, a        , 0.0e+00, b, NULL ) ;
        /* . Check the result. */
        for ( i = 0 ; i < o ; i++ ) Array2D_Item ( b, i, i ) -= 1.0e+00 ;
        maximumDeviation = RealArray2D_AbsoluteMaximum ( b ) ;
        /* . Finish up. */
        RealArray2D_Deallocate ( &a ) ;
        RealArray2D_Deallocate ( &b ) ;
    }
    return maximumDeviation ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Normalize the basis.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _DiagonalTolerance   1.0e-10
# define _EigenValueTolerance 1.0e-30
Real GaussianBasis_Normalize (       GaussianBasis *self               ,
                               const Boolean        checkNormalization ,
                                     Status        *status             )
{
    Real deviation = 0.0e+00 ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer          d, N, u, v ;
        auto Real             r[3] = { 0.0e+00, 0.0e+00, 0.0e+00 } ;
        auto RealArray2D     *c2s = NULL, *S  = NULL, *s2c = NULL, *X = NULL, *Y = NULL ;
        auto SymmetricMatrix *Ss  = NULL, *Sw = NULL ;
        /* . Deallocate the previous Xs. */
        RealArray2D_Deallocate ( &(self->c2o) ) ;
        RealArray2D_Deallocate ( &(self->o2c) ) ;
        /* . Fill the primitive ccbf for the basis. */
        GaussianBasis_FillPrimitiveCCBF ( self ) ;
        /* . Calculate the normalization S in the working basis. */
        S  = RealArray2D_AllocateWithExtents ( self->nbasisw, self->nbasisw, status ) ;
        Sw = SymmetricMatrix_AllocateWithExtent ( self->nbasisw, status ) ;
        if ( self->basisType == GaussianBasisType_Coulomb ) GaussianBasisIntegrals_2Coulomb ( self, r, self, r, S ) ;
        else                                                GaussianBasisIntegrals_2Overlap ( self, r, self, r, S ) ;
        for ( u = 0 ; u < self->nbasisw ; u++ )
        {
            for ( v = 0 ; v <= u ; v++ ) { SymmetricMatrix_Item ( Sw, u, v ) = Array2D_Item ( S, u, v ) ; }
        }
# ifdef _DebugPrint
printf ( "\nIntegrals (%u):\n", self->QSPHERICAL ) ; fflush ( stdout ) ;
SymmetricMatrix_Print ( Sw ) ; fflush ( stdout ) ;
# endif
        /* . Transform Sw to the spherical harmonic basis if necessary. */
        if ( self->QSPHERICAL )
        {
            auto Integer      i, iShell, j ;
            auto RealArray2D *mc, *ms ;
            /* . Create the forwards and backwards Xs. */
            c2s = RealArray2D_AllocateWithExtents ( self->nbasisw, self->nbasis, status ) ; RealArray2D_Set ( c2s, 0.0e+00 ) ;
            s2c = RealArray2D_AllocateWithExtents ( self->nbasisw, self->nbasis, status ) ; RealArray2D_Set ( s2c, 0.0e+00 ) ;
            for ( iShell = 0 ; iShell < self->nshells ; iShell++ )
            {
                mc = self->shells[iShell].c2s ;
                ms = self->shells[iShell].s2c ;
                if ( mc != NULL )
                {
                    for ( i = 0 ; i < View2D_Rows ( mc ) ; i++ )
                    {
                        for ( j = 0 ; j < View2D_Columns ( mc ) ; j++ )
                        {
                            Array2D_Item ( c2s, i + self->shells[iShell].nstartw, j + self->shells[iShell].nstart ) =  Array2D_Item ( mc, i, j ) ;
                        }
                    }
                }
                else
                {
                    for ( i = 0 ; i < self->shells[iShell].nbasisw ; i++ )
                    {
                        Array2D_Item ( c2s, i + self->shells[iShell].nstartw, i + self->shells[iShell].nstart ) = 1.0e+00 ;
                    }
                }
                if ( ms != NULL )
                {
                    for ( i = 0 ; i < View2D_Rows ( ms ) ; i++ )
                    {
                        for ( j = 0 ; j < View2D_Columns ( ms ) ; j++ )
                        {
                            Array2D_Item ( s2c, i + self->shells[iShell].nstartw, j + self->shells[iShell].nstart ) =  Array2D_Item ( ms, i, j ) ;
                        }
                    }
                }
                else
                {
                    for ( i = 0 ; i < self->shells[iShell].nbasisw ; i++ )
                    {
                        Array2D_Item ( s2c, i + self->shells[iShell].nstartw, i + self->shells[iShell].nstart ) = 1.0e+00 ;
                    }
                }
            }
            /*# define _CHECKCTOS*/
            # ifdef _CHECKCTOS
            { auto Real e = 0.0e+00 ; e = CheckOrthogonalization ( c2s, s2c, &e, NULL ) ; printf ( "\nC->S deviation %f %f\n", e ) ; }
            # endif
            /* . Transform Sw to Ss. */
            Ss = SymmetricMatrix_AllocateWithExtent ( self->nbasis, status ) ;
            SymmetricMatrix_Transform ( Sw, c2s, False, Ss, status ) ;
        }
        else { Ss = Sw ; Sw = NULL ; }
        /* . Allocate X. */
        d = SymmetricMatrix_Extent ( Ss ) ;
        X = RealArray2D_AllocateWithExtents ( d, d, status ) ;
        N = d ;
        /* . Diagonal normalization. */
        if ( self->normalizationType == NormalizationType_Diagonal )
        {
            auto Boolean isOK ;
            auto Real    dTolerance = _DiagonalTolerance ;
            /* . The basis should already be diagonal.
            !    In contrast to the general case, this option preserves basis function order.
            !    It is necessary for cases, such as MNDO, where this order is imposed elsewhere,
            !    and not by the basis itself.
            !    This should be improved - e.g. by normalizing Ss first, then checking its form.
            */
            isOK = SymmetricMatrix_IsDiagonal ( Ss, dTolerance ) ;
            if ( isOK )
            {
                auto Integer i ;
                auto Real    v ;
                RealArray2D_Set ( X, 0.0e+00 ) ;
                for ( i = 0 ; i < d ; i++ )
                {
                    v = SymmetricMatrix_Item ( Ss, i, i ) ;
                    if ( fabs ( v ) > _EigenValueTolerance ) { Array2D_Item ( X, i, i ) = 1.0e+00 / sqrt ( v ) ; }
                    else                                     { Status_Set ( status, Status_AlgorithmError ) ; break ; }
                }
            }
            else Status_Set ( status, Status_AlgorithmError ) ;
        }
        /* . Determine the X in the general case. */
        else
        {
            auto Boolean doCanonical = ( self->normalizationType == NormalizationType_Canonical ) ;
            auto Real    eTolerance  = _EigenValueTolerance ;
# ifdef _DebugPrint
doCanonical = True ;
# endif
            N = OrthogonalizingTransformation ( Ss, doCanonical, True, &eTolerance, NULL, NULL, X, status ) ;
            if ( N < d )
            {
                auto RealArray2D *Xt = NULL ;
                Xt = RealArray2D_Allocate ( status ) ;
                RealArray2D_View ( X, 0, 0, d, N, 1, 1, True, Xt, status ) ;
                RealArray2D_Deallocate ( &X ) ;
                X = Xt ;
            }
        }
        if ( Status_IsOK ( status ) )
        {
            /* . Create Y = S * X. */
            Y = RealArray2D_AllocateWithExtents ( SymmetricMatrix_Extent ( Ss ), N, status ) ;
            SymmetricMatrix_PostMatrixMultiply ( Ss, X, False, Y, status ) ;
            /* . Set the Xs. */
            if ( self->QSPHERICAL )
            {
                self->c2o = RealArray2D_AllocateWithExtents ( self->nbasisw, N, NULL ) ;
                self->o2c = RealArray2D_AllocateWithExtents ( self->nbasisw, N, NULL ) ;
                RealArray2D_MatrixMultiply ( False, False, 1.0e+00, c2s, X, 0.0e+00, self->c2o, status ) ;
                RealArray2D_MatrixMultiply ( False, False, 1.0e+00, s2c, Y, 0.0e+00, self->o2c, status ) ;
/*
printf ( "\nSc:"  ) ; SymmetricMatrix_Print ( Sw  ) ;
printf ( "\nSs:"  ) ; SymmetricMatrix_Print ( Ss  ) ;
printf ( "\nc2s:" ) ; RealArray2D_Print ( c2s ) ;
printf ( "\ns2c:" ) ; RealArray2D_Print ( s2c ) ;
printf ( "\nX:"   ) ; RealArray2D_Print ( X   ) ;
printf ( "\nY:"   ) ; RealArray2D_Print ( Y   ) ;
printf ( "\nc2o:" ) ; RealArray2D_Print ( self->c2o ) ;
printf ( "\no2c:" ) ; RealArray2D_Print ( self->o2c ) ;
*/
            }
            else { self->c2o = X ; self->o2c = Y ; }
/*
# define _TEST
# ifdef _TEST
if ( ( self->normalizationType != NormalizationType_Diagonal ) && ( self->QSPHERICAL ) )
{
RealArray2D_CopyTo ( c2s, self->c2o, status ) ;
RealArray2D_CopyTo ( s2c, self->o2c, status ) ;
}
# endif
*/
# ifdef _DebugPrint
printf ( "\nAtomic number and types = %d %u %u", self->atomicNumber, self->basisType, self->normalizationType ) ; fflush ( stdout ) ;
printf ( "\nc2o:" ) ; RealArray2D_Print ( self->c2o ) ; fflush ( stdout ) ;
printf ( "\no2c:" ) ; RealArray2D_Print ( self->o2c ) ; fflush ( stdout ) ;
# endif
            /* . Do the normalization check. */
            if ( checkNormalization )
            {
                auto Integer u, v ;
                auto Real    deviation1, deviation2 ;
                for ( u = 0 ; u < self->nbasisw ; u++ )
                {
                    for ( v = 0 ; v < u ; v++ ) { Array2D_Item ( S, v, u ) = Array2D_Item ( S, u, v ) ; }
                }
                deviation1 = CheckNormalization     ( self     , S        , status ) ;
                deviation2 = CheckOrthogonalization ( self->c2o, self->o2c, status ) ;
                deviation  = Maximum ( deviation1, deviation2 ) ;
            }
        }
        /* . Clean up. */
        RealArray2D_Deallocate     ( &c2s ) ;
        RealArray2D_Deallocate     ( &S   ) ;
        RealArray2D_Deallocate     ( &s2c ) ;
        SymmetricMatrix_Deallocate ( &Ss  ) ;
        SymmetricMatrix_Deallocate ( &Sw  ) ;
        if ( self->QSPHERICAL )
        {
            RealArray2D_Deallocate ( &Y ) ;
            RealArray2D_Deallocate ( &X ) ;
        }
    }
    return deviation ;
}
# undef _DiagonalTolerance
# undef _EigenValueTolerance
