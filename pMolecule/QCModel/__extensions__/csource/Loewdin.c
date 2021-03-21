/*==================================================================================================================================
! . Loewdin charge analysis.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>

# include "Integer.h"
# include "Loewdin.h"
# include "NumericalMacros.h"

/* . basisIndices refers to the orthogonal basis. */

/*----------------------------------------------------------------------------------------------------------------------------------
! . Atomic charges.
! . Charges incremented here.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Loewdin_AtomicCharges ( const IntegerArray1D  *basisIndices , /* . Orthogonal basis. */
                             const RealArray2D     *loewdinT     , /* . Nwork x Northogonal. */
                             const SymmetricMatrix *density      ,
                                   RealArray1D     *charges      )
{
    if ( ( basisIndices != NULL ) &&
         ( charges      != NULL ) &&
         ( density      != NULL ) &&
         ( loewdinT     != NULL ) )
    {
        auto Integer i, n, u, u0, u1, v, w ;
        auto Real    f, g ;
        n = View2D_Rows ( loewdinT ) ;
        for ( i = 0 ; i < View1D_Extent ( charges ) ; i++ )
        {
            u0 = Array1D_Item ( basisIndices, i   ) ;
            u1 = Array1D_Item ( basisIndices, i+1 ) ;
            for ( u = u0, f = 0.0e+00 ; u < u1 ; u++ )
            {
                for ( v = 0 ; v < n ; v++ )
                {
                    for ( w = 0, g = 0.0e+00 ; w < v ; w++ ) g += ( SymmetricMatrix_Item ( density, v, w ) * Array2D_Item ( loewdinT, w, u ) ) ;
                    for ( w = v              ; w < n ; w++ ) g += ( SymmetricMatrix_Item ( density, w, v ) * Array2D_Item ( loewdinT, w, u ) ) ;
                    f += ( Array2D_Item ( loewdinT, v, u ) * g ) ;
                }
            }
            Array1D_Item ( charges, i ) -= f ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Bond orders.
! . Bond orders incremented here.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Loewdin_BondOrders ( const IntegerArray1D  *basisIndices ,
                          const RealArray2D     *loewdinT     ,
                          const SymmetricMatrix *density      ,
                                SymmetricMatrix *bondOrders   ,
                                Status          *status       )
{
    if ( ( basisIndices != NULL ) &&
         ( bondOrders   != NULL ) &&
         ( density      != NULL ) &&
         ( loewdinT     != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto SymmetricMatrix *ps = NULL ;
        ps = SymmetricMatrix_AllocateWithExtent ( View2D_Columns ( loewdinT ), status ) ;
        SymmetricMatrix_Transform ( density, loewdinT, False, ps, status ) ;
        if ( Status_IsOK ( status ) )
        {
            auto Integer i, j, u, u0, u1, v, v0, v1 ;
            auto Real    sum ;
            for ( i = 0 ; i < SymmetricMatrix_Extent ( bondOrders ) ; i++ )
            {
                u0 = Array1D_Item ( basisIndices, i   ) ;
                u1 = Array1D_Item ( basisIndices, i+1 ) ;
                /* . Off-diagonal. */
                for ( j = 0 ; j < i ; j++ )
                {
                    v0 = Array1D_Item ( basisIndices, j   ) ;
                    v1 = Array1D_Item ( basisIndices, j+1 ) ;
                    for ( u = u0, sum = 0.0e+00 ; u < u1 ; u++ )
                    {
                        for ( v = v0 ; v < v1 ; v++ ) sum += pow ( SymmetricMatrix_Item ( ps, u, v ), 2 ) ;
                    }
                    SymmetricMatrix_Item ( bondOrders, i, j ) += sum ;
                }
                /* . Diagonal. */
                for ( u = u0, sum = 0.0e+00 ; u < u1 ; u++ )
                {
                    for ( v = u0 ; v < u ; v++ ) sum += ( 2.0e+00 * pow ( SymmetricMatrix_Item ( ps, u, v ), 2 ) ) ;
                    sum += pow ( SymmetricMatrix_Item ( ps, u, u ), 2 ) ;
                }
                SymmetricMatrix_Item ( bondOrders, i, i ) += sum ;
            }
        }
        SymmetricMatrix_Deallocate ( &ps ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Charge density derivatives.
! . Fock incremented here only.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Loewdin_ChargeDensityDerivatives ( const IntegerArray1D  *basisIndices   ,
                                        const RealArray1D     *potentials     , /* . = dXdQ. */
                                        const RealArray2D     *loewdinT       ,
                                              SymmetricMatrix *fock           )
{
    if ( ( basisIndices != NULL ) &&
         ( fock         != NULL ) &&
         ( loewdinT     != NULL ) &&
         ( potentials   != NULL ) )
    {
        auto Integer i, u, u0, u1, v, w ;
        auto Real    f, p ;
        for ( v = 0 ; v < SymmetricMatrix_Extent ( fock ) ; v++ )
        {
            for ( w = 0 ; w <= v ; w++ )
            {
                for ( i = 0, f = 0.0e+00 ; i < View1D_Extent ( potentials ) ; i++ )
                {
                    p  = Array1D_Item ( potentials  , i   ) ;
                    u0 = Array1D_Item ( basisIndices, i   ) ;
                    u1 = Array1D_Item ( basisIndices, i+1 ) ;
                    for ( u = u0 ; u < u1 ; u++ ) f += ( p * Array2D_Item ( loewdinT, v, u ) * Array2D_Item ( loewdinT, w, u ) ) ;
                }
                SymmetricMatrix_Item ( fock, v, w ) -= f ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Weighted density matrix.
! . wDensity incremented here only.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _EigenValueTolerance 1.0e-10
void Loewdin_WeightedDensity ( const IntegerArray1D  *basisIndices        ,
                               const RealArray1D     *potentials          , /* . = dXdQ. */
                               const RealArray1D     *eigenValues         ,
                               const RealArray2D     *eigenVectors        ,
                               const RealArray2D     *loewdinT            ,
                               const RealArray2D     *o2C                 ,
                               const RealArray2D     *c2O                 ,
                               const SymmetricMatrix *density             ,
                                     Real            *eigenValueTolerance ,
                                     SymmetricMatrix *wDensity            ,
                                     Status          *status              )
{
    if ( ( basisIndices != NULL ) &&
         ( c2O          != NULL ) &&
         ( density      != NULL ) &&
         ( eigenValues  != NULL ) &&
         ( eigenVectors != NULL ) &&
         ( loewdinT     != NULL ) &&
         ( potentials   != NULL ) &&
         ( o2C          != NULL ) &&
         ( wDensity     != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Boolean isOK ;
        auto Integer nA, nC, nE, nO ;
        /* . Get and check dimensions. */
        nA   = View1D_Extent          ( potentials  ) ;
        nC   = SymmetricMatrix_Extent ( density     ) ;
        nE   = View1D_Extent          ( eigenValues ) ;
        nO   = View2D_Columns         ( loewdinT    ) ;
        isOK = ( View1D_Extent          ( basisIndices     ) == ( nA + 1 ) ) &&
               ( View2D_Columns         ( c2O              ) ==   nO       ) &&
               ( View2D_Rows            ( c2O              ) ==   nC       ) &&
               ( View2D_Columns         ( eigenVectors     ) ==   nE       ) &&
               ( View2D_Rows            ( eigenVectors     ) ==   nO       ) &&
               ( View2D_Rows            ( loewdinT         ) ==   nC       ) &&
               ( View2D_Columns         ( o2C              ) ==   nO       ) &&
               ( View2D_Rows            ( o2C              ) ==   nC       ) &&
               ( Array1D_Item           ( basisIndices, nA ) ==   nO       ) &&
               ( SymmetricMatrix_Extent ( wDensity         ) ==   nC       ) ;
        if ( ! isOK ) Status_Set ( status, Status_NonConformableArrays ) ;
        if ( Status_IsOK ( status ) )
        {
            auto RealArray1D     *tempEV ;
            auto RealArray2D     *pO2c ;
            auto SymmetricMatrix *tempNC, *tempNE, *tempNO ;
            /* . Allocate space. */
            tempEV = RealArray1D_AllocateWithExtent     ( nE    , status ) ;
            pO2c   = RealArray2D_AllocateWithExtents    ( nC, nO, status ) ;
            tempNC = SymmetricMatrix_AllocateWithExtent ( nC    , status ) ;
            tempNE = SymmetricMatrix_AllocateWithExtent ( nE    , status ) ;
            tempNO = SymmetricMatrix_AllocateWithExtent ( nO    , status ) ;
            if ( Status_IsOK ( status ) )
            {
                auto Integer i, j, u, u0, u1, v, v0, v1, w ;
                auto Real    a, ab, b, f, tolerance ;
                if ( eigenValueTolerance == NULL ) tolerance =  _EigenValueTolerance  ;
                else                               tolerance = (*eigenValueTolerance) ;
                /* . P * o2C. */
/*# define _DebugPrint*/
# ifdef _DebugPrint
printf ( "\nDensity:\n" ) ; fflush ( stdout ) ;
SymmetricMatrix_Print ( density ) ; fflush ( stdout ) ;
if ( o2C != NULL )
{
printf ( "O->C:\n" ) ; fflush ( stdout ) ;
RealArray2D_Print ( o2C ) ; fflush ( stdout ) ;
}
if ( c2O != NULL )
{
printf ( "C->O:\n" ) ; fflush ( stdout ) ;
RealArray2D_Print ( c2O ) ; fflush ( stdout ) ;
}
printf ( "Potentials:\n" ) ; fflush ( stdout ) ;
RealArray1D_Print ( potentials ) ; fflush ( stdout ) ;
printf ( "Loewdin Transformation:\n" ) ; fflush ( stdout ) ;
RealArray2D_Print ( loewdinT ) ; fflush ( stdout ) ;
printf ( "Eigenvalues\n" ) ; fflush ( stdout ) ;
RealArray1D_Print ( eigenValues ) ; fflush ( stdout ) ;
printf ( "Eigenvectors\n" ) ; fflush ( stdout ) ;
RealArray2D_Print ( eigenVectors ) ; fflush ( stdout ) ;
# endif
                SymmetricMatrix_PostMatrixMultiply ( density, o2C, False, pO2c, status ) ;
# ifdef _DebugPrint
printf ( "pO2c:\n" ) ; fflush ( stdout ) ;
RealArray2D_Print ( pO2c ) ; fflush ( stdout ) ;
# endif
                /* . Symmetrized core matrix. */
                for ( i = 0 ; i < View1D_Extent ( potentials ) ; i++ )
                {
                    a  = Array1D_Item ( potentials  , i   ) ;
                    u0 = Array1D_Item ( basisIndices, i   ) ;
                    u1 = Array1D_Item ( basisIndices, i+1 ) ;
                    for ( u = u0 ; u < u1 ; u++ )
                    {
                         for ( j = 0 ; j <= i ; j++ )
                        {
                            b  = Array1D_Item ( potentials  , j ) ;
                            v0 = Array1D_Item ( basisIndices, j ) ;
                            v1 = Minimum ( u, Array1D_Item ( basisIndices, j+1 ) - 1 ) ;
                            for ( v = v0 ; v <= v1 ; v++ )
                            {
                                for ( w = 0, f = 0.0e+00 ; w < nC ; w++ )
                                {
                                    f += b * Array2D_Item ( pO2c, w, u ) * Array2D_Item ( loewdinT, w, v ) +
                                         a * Array2D_Item ( pO2c, w, v ) * Array2D_Item ( loewdinT, w, u ) ;
                                }
                                SymmetricMatrix_Item ( tempNO, u, v ) = f ;
                            }
                        }
                    }
                }
# ifdef _DebugPrint
printf ( "\nSymmetrized Core Matrix:\n" ) ; fflush ( stdout ) ;
SymmetricMatrix_Print ( tempNO ) ; fflush ( stdout ) ;
# endif
                /* . First transformation. */
                SymmetricMatrix_Transform ( tempNO, eigenVectors, False, tempNE, status ) ;
# ifdef _DebugPrint
printf ( "\nTransformed Symmetrized Core Matrix:\n" ) ; fflush ( stdout ) ;
SymmetricMatrix_Print ( tempNE ) ; fflush ( stdout ) ;
# endif
                /* . Get the square roots of the eigenvalues. */
                for ( u = 0 ; u < nE ; u++ )
                {
                    a = Array1D_Item ( eigenValues, u ) ;
                    Array1D_Item ( tempEV, u ) = ( a < 0.0e+00 ? 0.0e+00 : sqrt ( a ) ) ;
                }
                /* . Scale by the eigenvalue factors. */
                for ( u = 0 ; u < nE ; u++ )
                {
                    a = Array1D_Item ( tempEV, u ) ;
                    for ( v = 0 ; v <= u ; v++ )
                    {
                        b  = Array1D_Item ( tempEV, v ) ;
                        ab = a + b ;
                        if ( ab > tolerance ) SymmetricMatrix_Item ( tempNE, u, v ) /= ab      ;
                        else                  SymmetricMatrix_Item ( tempNE, u, v )  = 0.0e+00 ;
                    }
                }
# ifdef _DebugPrint
printf ( "\nScaled Symmetrized Core Matrix:\n" ) ; fflush ( stdout ) ;
SymmetricMatrix_Print ( tempNE ) ; fflush ( stdout ) ;
# endif
                /* . Second transformation. */
                SymmetricMatrix_Transform ( tempNE, eigenVectors, True, tempNO, status ) ;
# ifdef _DebugPrint
printf ( "\nSecond Transformed Symmetrized Core Matrix:\n" ) ; fflush ( stdout ) ;
SymmetricMatrix_Print ( tempNO ) ; fflush ( stdout ) ;
# endif
                /* . Back transform temp by c2O * temp * c2O^T. */
                SymmetricMatrix_Transform ( tempNO, c2O, True, tempNC, status ) ;
                /* . Add in the contributions to the density. */
                SymmetricMatrix_Add ( wDensity, -2.0e+00, tempNC, status ) ;
# ifdef _DebugPrint
printf ( "WDM Contribution\n" ) ; fflush ( stdout ) ;
SymmetricMatrix_Print ( tempNC ) ; fflush ( stdout ) ;
printf ( "WDM Final\n" ) ; fflush ( stdout ) ;
SymmetricMatrix_Print ( wDensity ) ; fflush ( stdout ) ;
# endif
            }
            /* . Deallocate space. */
            RealArray1D_Deallocate     ( &tempEV ) ;
            RealArray2D_Deallocate     ( &pO2c   ) ;
            SymmetricMatrix_Deallocate ( &tempNC ) ;
            SymmetricMatrix_Deallocate ( &tempNE ) ;
            SymmetricMatrix_Deallocate ( &tempNO ) ;
        }
    }
}
# undef _EigenValueTolerance
