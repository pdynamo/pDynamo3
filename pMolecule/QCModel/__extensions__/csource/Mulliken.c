/*==================================================================================================================================
! . Mulliken charge analysis.
!=================================================================================================================================*/

# include "Integer.h"
# include "Mulliken.h"
# include "Real.h"

/* . Basis indices refer to whatever representation the density, overlap and Fock matrices are in. */

/*----------------------------------------------------------------------------------------------------------------------------------
! . Atomic charges.
! . Charges incremented here.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Mulliken_AtomicCharges ( const IntegerArray1D  *basisIndices ,
                              const SymmetricMatrix *density      ,
                              const SymmetricMatrix *overlap      ,
                                    RealArray1D     *charges      )
{
    if ( ( basisIndices != NULL ) &&
         ( charges      != NULL ) &&
         ( density      != NULL ) &&
         ( overlap      != NULL ) )
    {
        auto Integer i, n, u, u0, u1, v ;
        auto Real    ps ;
        n = SymmetricMatrix_Extent ( density ) ;
        for ( i = 0 ; i < View1D_Extent ( charges ) ; i++ )
        {
            u0 = Array1D_Item ( basisIndices, i   ) ;
            u1 = Array1D_Item ( basisIndices, i+1 ) ;
            for ( u = u0, ps = 0.0e+00 ; u < u1 ; u++ )
            {
                for ( v = 0 ; v < u ; v++ ) ps += ( SymmetricMatrix_Item ( density, u, v ) * SymmetricMatrix_Item ( overlap, u, v ) ) ;
                for ( v = u ; v < n ; v++ ) ps += ( SymmetricMatrix_Item ( density, v, u ) * SymmetricMatrix_Item ( overlap, v, u ) ) ;
            }
            Array1D_Item ( charges, i ) -= ps ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Bond orders.
! . Bond orders incremented here.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Mulliken_BondOrders ( const IntegerArray1D  *basisIndices ,
                           const SymmetricMatrix *density      ,
                           const SymmetricMatrix *overlap      ,
                                 SymmetricMatrix *bondOrders   ,
                                 Status          *status       )
{
    if ( ( basisIndices != NULL ) &&
         ( bondOrders   != NULL ) &&
         ( density      != NULL ) &&
         ( overlap      != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer      n  = SymmetricMatrix_Extent ( overlap ) ;
        auto RealArray2D *ps = NULL ;
        ps = RealArray2D_AllocateWithExtents    ( n, n, status ) ;
        SymmetricMatrix_SymmetricMatrixMultiply ( density, overlap, ps, status ) ;
        if ( Status_IsOK ( status ) )
        {
            auto Integer i, j, u, u0, u1, v, v0, v1 ;
            auto Real    sum ;
            for ( i = 0 ; i < SymmetricMatrix_Extent ( bondOrders ) ; i++ )
            {
                u0 = Array1D_Item ( basisIndices, i   ) ;
                u1 = Array1D_Item ( basisIndices, i+1 ) ;
                /* . Off-diagonal and diagonal. */
                for ( j = 0 ; j <= i ; j++ )
                {
                    v0 = Array1D_Item ( basisIndices, j   ) ;
                    v1 = Array1D_Item ( basisIndices, j+1 ) ;
                    for ( u = u0, sum = 0.0e+00 ; u < u1 ; u++ )
                    {
                        for ( v = v0 ; v < v1 ; v++ ) sum += ( Array2D_Item ( ps, u, v ) * Array2D_Item ( ps, v, u ) ) ;
                    }
                    SymmetricMatrix_Item ( bondOrders, i, j ) += sum ;
                }
            }
        }
        RealArray2D_Deallocate ( &ps ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Charge density derivatives.
! . Fock incremented here only.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Mulliken_ChargeDensityDerivatives ( const IntegerArray1D  *basisIndices   ,
                                         const RealArray1D     *potentials     , /* . = dXdQ. */
                                         const SymmetricMatrix *overlap        ,
                                               SymmetricMatrix *fock           )
{
    if ( ( basisIndices != NULL ) &&
         ( fock         != NULL ) &&
         ( overlap      != NULL ) &&
         ( potentials   != NULL ) )
    {
        auto Integer i, j, u, u0, u1, v, v0, v1 ;
        auto Real    p, pq, q ;
        for ( i = 0 ; i < View1D_Extent ( potentials ) ; i++ )
        {
            p  = Array1D_Item ( potentials  , i   ) ;
            u0 = Array1D_Item ( basisIndices, i   ) ;
            u1 = Array1D_Item ( basisIndices, i+1 ) ;
            /* . Off-diagonal. */
            for ( j = 0 ; j < i ; j++ )
            {
                q  = Array1D_Item ( potentials  , j   ) ;
                v0 = Array1D_Item ( basisIndices, j   ) ;
                v1 = Array1D_Item ( basisIndices, j+1 ) ;
                pq = 0.5e+00 * ( p + q ) ;
                for ( u = u0 ; u < u1 ; u++ )
                {
                    for ( v = v0 ; v < v1 ; v++ ) SymmetricMatrix_Item ( fock, u, v ) -= ( pq * SymmetricMatrix_Item ( overlap, u, v ) ) ;
                }
            }
            /* . Diagonal. */
            for ( u = u0 ; u < u1 ; u++ )
            {
                for ( v = u0 ; v <= u ; v++ ) SymmetricMatrix_Item ( fock, u, v ) -= ( p * SymmetricMatrix_Item ( overlap, u, v ) ) ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Charge restraint W-matrix and core term.
! . This is hugely wasteful for MNDO and potentially very wasteful for Mulliken methods as W is sparse!
! . However it is done for the moment to simplify the charge restraint code, in particular for those methods,
! . such as DFT with Loewdin charges, for which W is dense.
! . Only basic checking is done.
! . The input W matrix is initialized on entry.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Mulliken_ChargeRestraintMatrix ( const IntegerArray1D  *basisIndices   ,
                                      const RealArray1D     *nuclearCharges ,
                                      const IntegerArray1D  *crIndices      ,
                                      const RealArray1D     *crWeights      ,
                                      const Boolean          isSpin         ,
                                      const SymmetricMatrix *overlap        ,
                                            SymmetricMatrix *W              )
{
    Real core = 0.0e+00 ;
    if ( ( basisIndices   != NULL ) &&
         ( crIndices      != NULL ) &&
         ( crWeights      != NULL ) &&
         ( nuclearCharges != NULL ) &&
         ( overlap        != NULL ) &&
         ( W              != NULL ) )
    {
        auto Integer a, b, i, u, u0, u1, v, v0, v1 ;
        auto Real    w ;
        SymmetricMatrix_Set ( W, 0.0e+00 ) ;
        for ( i = 0 ; i < View1D_Extent ( crIndices ) ; i++ )
        {
            a = Array1D_Item ( crIndices, i ) ;
            w = Array1D_Item ( crWeights, i ) ;
            if ( ! isSpin )
            {    
                 core += ( w * Array1D_Item ( nuclearCharges, a ) ) ;
                 w    *= -1.0e+00 ; /* . w is -1.0 for electrons in this case. */
            }
            u0 = Array1D_Item ( basisIndices, a   ) ;
            u1 = Array1D_Item ( basisIndices, a+1 ) ;
            /* . AA block. */
            for ( u = u0 ; u < u1 ; u++ )
            {
                for ( v = u0 ; v < u ; v++ )
                {
                    SymmetricMatrix_Item ( W, u, v ) += ( w * SymmetricMatrix_Item ( overlap, u, v ) ) ;
                }
                SymmetricMatrix_Item ( W, u, u ) += ( w * SymmetricMatrix_Item ( overlap, u, u ) ) ;
            }
            /* . AB blocks ( A > B ). */
            for ( b = 0 ; b < a ; b++ )
            {
                v0 = Array1D_Item ( basisIndices, b   ) ;
                v1 = Array1D_Item ( basisIndices, b+1 ) ;
                for ( u = u0 ; u < u1 ; u++ )
                {
                    for ( v = v0 ; v < v1 ; v++ )
                    {
                        SymmetricMatrix_Item ( W, u, v ) += ( 0.5 * w * SymmetricMatrix_Item ( overlap, u, v ) ) ;
                    }
                }
            }
            /* . AB blocks ( A < B ). */
            for ( b = a+1 ; b < View1D_Extent ( basisIndices ) - 1 ; b++ )
            {
                v0 = Array1D_Item ( basisIndices, b   ) ;
                v1 = Array1D_Item ( basisIndices, b+1 ) ;
                for ( v = v0 ; v < v1 ; v++ )
                {
                    for ( u = u0 ; u < u1 ; u++ )
                    {
                        SymmetricMatrix_Item ( W, v, u ) += ( 0.5 * w * SymmetricMatrix_Item ( overlap, v, u ) ) ;
                    }
                }
            }
        }
    }
    return core ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Charge restraint weighted density.
! . This method needs to be called for each restraint separately together with the appropriate derivative of the restraint
! . energy model with respect to the restraint, dRdL, and the density matrix, density. The latter will be the total density
! . for a charge restraint (isSpin = False) or the spin density for a spin restraint (isSpin = True).
! . The weighted density is incremented only and so should be initialized before entry.
! . This is essentially the same as the previous method with the density in place of the overlap
! . except that the diagonal blocks are ignored as the diagonal overlap integrals are zero.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _Factor 2.0e+00 /* . The appropriate weight factor. */
void Mulliken_ChargeRestraintWeightedDensity ( const IntegerArray1D  *basisIndices ,
                                               const IntegerArray1D  *crIndices    ,
                                               const RealArray1D     *crWeights    ,
                                               const Boolean          isSpin       ,
                                               const Real             dRdL         ,
                                               const SymmetricMatrix *density      ,
                                                     SymmetricMatrix *wdm          )
{
    if ( ( basisIndices != NULL    ) &&
         ( crIndices    != NULL    ) &&
         ( crWeights    != NULL    ) &&
         ( density      != NULL    ) &&
         ( wdm          != NULL    ) &&
         ( dRdL         != 0.0e+00 ) ) /* . Do nothing if dRdL is zero! */
    {
        auto Integer a, b, i, u, u0, u1, v, v0, v1 ;
        auto Real    w ;
        for ( i = 0 ; i < View1D_Extent ( crIndices ) ; i++ )
        {
            a =                        Array1D_Item ( crIndices, i ) ;
            w = 0.5 * _Factor * dRdL * Array1D_Item ( crWeights, i ) ; /* . 0.5 here as diagonal blocks zero. */
            if ( ! isSpin ) w *= -1.0e+00 ;
            u0 = Array1D_Item ( basisIndices, a   ) ;
            u1 = Array1D_Item ( basisIndices, a+1 ) ;
            /* . AB blocks ( A > B ). */
            for ( b = 0 ; b < a ; b++ )
            {
                v0 = Array1D_Item ( basisIndices, b   ) ;
                v1 = Array1D_Item ( basisIndices, b+1 ) ;
                for ( u = u0 ; u < u1 ; u++ )
                {
                    for ( v = v0 ; v < v1 ; v++ )
                    {
                        SymmetricMatrix_Item ( wdm, u, v ) += ( w * SymmetricMatrix_Item ( density, u, v ) ) ;
                    }
                }
            }
            /* . AB blocks ( A < B ). */
            for ( b = a+1 ; b < View1D_Extent ( basisIndices ) - 1 ; b++ )
            {
                v0 = Array1D_Item ( basisIndices, b   ) ;
                v1 = Array1D_Item ( basisIndices, b+1 ) ;
                for ( v = v0 ; v < v1 ; v++ )
                {
                    for ( u = u0 ; u < u1 ; u++ )
                    {
                        SymmetricMatrix_Item ( wdm, v, u ) += ( w * SymmetricMatrix_Item ( density, v, u ) ) ;
                    }
                }
            }
        }
    }
}
# undef _Factor

/*----------------------------------------------------------------------------------------------------------------------------------
! . Weighted density matrix.
! . wDensity incremented here only.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Mulliken_WeightedDensity ( const IntegerArray1D  *basisIndices   ,
                                const RealArray1D     *potentials     , /* . = dXdQ. */
                                const SymmetricMatrix *density        ,
                                      SymmetricMatrix *wDensity       )
{
    if ( ( basisIndices != NULL ) &&
         ( density      != NULL ) &&
         ( potentials   != NULL ) &&
         ( wDensity     != NULL ) )
    {
        auto Integer i, j, u, u0, u1, v, v0, v1 ;
        auto Real    p, pq, q ;
        for ( i = 0 ; i < View1D_Extent ( potentials ) ; i++ )
        {
            p  = Array1D_Item ( potentials  , i   ) ;
            u0 = Array1D_Item ( basisIndices, i   ) ;
            u1 = Array1D_Item ( basisIndices, i+1 ) ;
            /* . Off-diagonal only. */
            for ( j = 0 ; j < i ; j++ )
            {
                q  = Array1D_Item ( potentials  , j   ) ;
                v0 = Array1D_Item ( basisIndices, j   ) ;
                v1 = Array1D_Item ( basisIndices, j+1 ) ;
                pq = - ( p + q ) ;
                for ( u = u0 ; u < u1 ; u++ )
                {
                    for ( v = v0 ; v < v1 ; v++ ) SymmetricMatrix_Item ( wDensity, u, v ) += ( pq * SymmetricMatrix_Item ( density, u, v ) ) ;
                }
            }
        }
    }
}
