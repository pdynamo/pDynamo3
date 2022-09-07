/*==================================================================================================================================
! . MNDO atomic multipoles up to quadrupole.
!=================================================================================================================================*/

# include <math.h>

# include "Integer.h"
# include "MNDODefinitions.h"
# include "MNDOMultipoles.h"
# include "Real.h"

/*
  Conversions MNDOD / pDynamo:

    dA = DD(2)          = ddp[1]
    qA = DD(3)**2       = ddp[2]**2 / 2.0
    qB = DD(4)**2 / 2.0 = ddp[3]**2 / 2.0
    dB = DD(5)          = ddp[4]
    qC = DD(6)**2 / 2.0 = ddp[5]**2 / 2.0
*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Atomic charges, dipoles and quadrupoles.
! . Note nuclear charges are excluded!
! . Multipoles is a flattened Nmult x Nqc row-major matrix.
! . Multipoles set here.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDO_AtomicMultipoles ( const MNDOParametersContainer *parameters              ,
                             const IntegerArray1D          *basisIndices            ,
                             const SymmetricMatrix         *density                 ,
                             const MultipoleRepresentation  multipoleRepresentation ,
                             const Integer                  multipoleOrder          ,
                                   RealArray1D             *multipoles              )
{
    RealArray1D_Set ( multipoles, 0.0e+00 ) ;
    if ( ( parameters   != NULL ) &&
         ( basisIndices != NULL ) &&
         ( density      != NULL ) &&
         ( multipoles   != NULL ) )
    {
        auto Boolean         doDipoles = ( multipoleOrder > 0 ), doQuadrupoles = ( multipoleOrder > 1 ) ;
        auto Integer         i, i0, nI, nP, u ;
        auto MNDOParameters *iData ;
        auto Real            iSqrt3, q, sqrt3 ;
        sqrt3  = sqrt ( 3.0e+00 ) ;
        iSqrt3 = 1.0e+00 / sqrt3 ;
        nP     = parameters->capacity ;
        for ( i = 0 ; i < nP ; i++ )
        {
            iData = parameters->entries[i] ;
            i0    = Array1D_Item ( basisIndices, i ) ;
            nI    = iData->norbitals ;
            /* . Charge. */
            q = 0.0e+00 ;
            for ( u = 0 ; u < nI ; u++ ) q -= SymmetricMatrix_Item ( density, i0 + u, i0 + u ) ;
            Array1D_Item ( multipoles, i ) = q ;
            /* . p-orbital contributions. */
            if ( doDipoles && ( iData->norbitals >= 4 ) )
            {
                /* . Orbital orders S < PZ < PX < PY. */
                auto Integer _PX = i0 + PX, _PY = i0 + PY, _PZ = i0 + PZ, _S = i0 ;
                auto Real dA, dX, dY, dZ, qA, qT, qXX = 0.0e+00 , qXY = 0.0e+00 , qXZ = 0.0e+00 ,
                                                  qYY = 0.0e+00 , qYZ = 0.0e+00 , qZZ = 0.0e+00 ;
                dA = iData->ddp[1] ;
                qA = iData->ddp[2] * iData->ddp[2] * 0.5e+00 ;
                /* . Dipole - sp. */
                dX = dA * SymmetricMatrix_Item ( density, _PX, _S ) ;
                dY = dA * SymmetricMatrix_Item ( density, _PY, _S ) ;
                dZ = dA * SymmetricMatrix_Item ( density, _PZ, _S ) ;
                if ( doQuadrupoles )
                {
                    /* . Quadrupole - pp. */
                    qXX = qA * SymmetricMatrix_Item ( density, _PX, _PX ) ;
                    qXY = qA * SymmetricMatrix_Item ( density, _PY, _PX ) ;
                    qXZ = qA * SymmetricMatrix_Item ( density, _PX, _PZ ) ;
                    qYY = qA * SymmetricMatrix_Item ( density, _PY, _PY ) ;
                    qYZ = qA * SymmetricMatrix_Item ( density, _PY, _PZ ) ;
                    qZZ = qA * SymmetricMatrix_Item ( density, _PZ, _PZ ) ;
                }
                /* . d-orbital contributions. */
                if ( iData->norbitals >= 9 )
                {
                    /* . Orbital orders DZ2 < DXZ < DYZ < DX2Y2 < DXY. */
                    auto Integer _DX2Y2 = i0 + DX2Y2, _DXY = i0 + DXY, _DXZ = i0 + DXZ, _DYZ = i0 + DYZ, _DZ2 = i0 + DZ2 ;
                    auto Real    dB, qB, qC ;
                    dB = iData->ddp[4] ;
                    qB = iData->ddp[3] * iData->ddp[3] * 0.5e+00 ;
                    qC = iData->ddp[5] * iData->ddp[5] * 0.5e+00 ;
                    /* . Dipole - pd. */
                    dX += dB * (             SymmetricMatrix_Item ( density, _DXZ,   _PZ ) +
                                             SymmetricMatrix_Item ( density, _DX2Y2, _PX ) +
                                             SymmetricMatrix_Item ( density, _DXY,   _PY ) -
                                    iSqrt3 * SymmetricMatrix_Item ( density, _DZ2,   _PX ) ) ;
                    dY += dB * (             SymmetricMatrix_Item ( density, _DYZ,   _PZ ) -
                                             SymmetricMatrix_Item ( density, _DX2Y2, _PY ) +
                                             SymmetricMatrix_Item ( density, _DXY,   _PX ) -
                                    iSqrt3 * SymmetricMatrix_Item ( density, _DZ2,   _PY ) ) ;
                    dZ += dB * (             SymmetricMatrix_Item ( density, _DXZ,   _PX ) +
                                             SymmetricMatrix_Item ( density, _DYZ,   _PY ) +
                          2.0e+00 * iSqrt3 * SymmetricMatrix_Item ( density, _DZ2,   _PZ ) ) ;
                    if ( doQuadrupoles )
                    {
                        /* . Quadrupole - sd. */
                        qXX += qB         * SymmetricMatrix_Item ( density, _DX2Y2, _S ) ;
                        qXY += qB         * SymmetricMatrix_Item ( density, _DXY  , _S ) ;
                        qXZ += qB         * SymmetricMatrix_Item ( density, _DXZ  , _S ) ;
                        qYY -= qB         * SymmetricMatrix_Item ( density, _DX2Y2, _S ) ;
                        qYZ += qB         * SymmetricMatrix_Item ( density, _DYZ  , _S ) ;
                        qZZ += qB * sqrt3 * SymmetricMatrix_Item ( density, _DZ2  , _S ) ;
                        /* . Quadrupole - dd. */
                        qXX -= qC * (                    SymmetricMatrix_Item ( density, _DYZ  , _DYZ   ) +
                                      2.0e+00 * iSqrt3 * SymmetricMatrix_Item ( density, _DX2Y2, _DZ2   ) ) ;
                        qXY += qC * (                    SymmetricMatrix_Item ( density, _DYZ  , _DXZ   ) -
                                      2.0e+00 * iSqrt3 * SymmetricMatrix_Item ( density, _DXY  , _DZ2   ) ) ;
                        qXZ += qC * (                    SymmetricMatrix_Item ( density, _DX2Y2, _DXZ   ) +
                                                iSqrt3 * SymmetricMatrix_Item ( density, _DXZ  , _DZ2   ) +
                                                         SymmetricMatrix_Item ( density, _DXY  , _DYZ   ) ) ;
                        qYY -= qC * (                    SymmetricMatrix_Item ( density, _DXZ  , _DXZ   ) -
                                      2.0e+00 * iSqrt3 * SymmetricMatrix_Item ( density, _DX2Y2, _DZ2   ) ) ;
                        qYZ += qC * (                    SymmetricMatrix_Item ( density, _DXY  , _DXZ   ) -
                                                         SymmetricMatrix_Item ( density, _DX2Y2, _DYZ   ) +
                                                iSqrt3 * SymmetricMatrix_Item ( density, _DYZ  , _DZ2   ) ) ;
                        qZZ -= qC * (                    SymmetricMatrix_Item ( density, _DX2Y2, _DX2Y2 ) +
                                                         SymmetricMatrix_Item ( density, _DXY  , _DXY   ) -
                                                         SymmetricMatrix_Item ( density, _DZ2  , _DZ2   ) ) ;
                    }
                }
                Array1D_Item ( multipoles, i +   nP ) = dX ;
                Array1D_Item ( multipoles, i + 2*nP ) = dY ;
                Array1D_Item ( multipoles, i + 3*nP ) = dZ ;
                if ( doQuadrupoles )
                {
                    qT = qXX + qYY + qZZ ;
                    if ( multipoleRepresentation == MultipoleRepresentation_Buckingham )
                    {
                        Array1D_Item ( multipoles, i + 4*nP ) = 0.5e+00 * ( 3.0e+00 * qXX - qT ) ;
                        Array1D_Item ( multipoles, i + 5*nP ) = 1.5e+00 * qXY ;
                        Array1D_Item ( multipoles, i + 6*nP ) = 1.5e+00 * qXZ ;
                        Array1D_Item ( multipoles, i + 7*nP ) = 0.5e+00 * ( 3.0e+00 * qYY - qT ) ;
                        Array1D_Item ( multipoles, i + 8*nP ) = 1.5e+00 * qYZ ;
                        Array1D_Item ( multipoles, i + 9*nP ) = 0.5e+00 * ( 3.0e+00 * qZZ - qT ) ;
                    }
                    else if ( multipoleRepresentation == MultipoleRepresentation_Spherical )
                    {
                        Array1D_Item ( multipoles, i + 4*nP ) = 0.5e+00 * ( 3.0e+00 * qZZ - qT ) ;
                        Array1D_Item ( multipoles, i + 5*nP ) = sqrt3 * qXZ ;
                        Array1D_Item ( multipoles, i + 6*nP ) = sqrt3 * qYZ ;
                        Array1D_Item ( multipoles, i + 7*nP ) = 0.5e+00 * sqrt3 * ( qXX - qYY ) ;
                        Array1D_Item ( multipoles, i + 8*nP ) = sqrt3 * qXY ;
                    }
                    else
                    {
                        Array1D_Item ( multipoles, i + 4*nP ) = qXX ;
                        Array1D_Item ( multipoles, i + 5*nP ) = qXY ;
                        Array1D_Item ( multipoles, i + 6*nP ) = qXZ ;
                        Array1D_Item ( multipoles, i + 7*nP ) = qYY ;
                        Array1D_Item ( multipoles, i + 8*nP ) = qYZ ;
                        Array1D_Item ( multipoles, i + 9*nP ) = qZZ ;
                    }
                }
            }
        }
        /* . Scaling of higher multipoles. */
        if ( doDipoles )
        {
            RealArray1D notCharges ;
            RealArray1D_View  ( multipoles, nP, ( View1D_Extent ( multipoles ) - nP ), 1, False, &notCharges, NULL ) ;
            RealArray1D_Scale ( &notCharges, -2.0e+00 ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Fock contribution of atomic charges, dipoles and quadrupoles to energy - Cartesian multipoles only.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Energy is C^T * P, Fock is initialized on entry, and potentials stored as flattened Nmult * Natoms matrix. */
/* . Require factors of two for diagonal terms, but not off-diagonal ones. */
void MNDO_AtomicMultipolesFock ( const MNDOParametersContainer *parameters     ,
                                 const IntegerArray1D          *basisIndices   ,
                                 const RealArray1D             *potentials     ,
                                 const Integer                  multipoleOrder ,
                                       SymmetricMatrix         *fock           )
{
    if ( ( parameters   != NULL ) &&
         ( basisIndices != NULL ) &&
         ( potentials   != NULL ) &&
         ( fock         != NULL ) )
    {
        auto Boolean         doDipoles = ( multipoleOrder > 0 ), doQuadrupoles = ( multipoleOrder > 1 ) ;
        auto Integer         i, i0, nI, nP, u ;
        auto MNDOParameters *iData ;
        auto Real            iSqrt3, sqrt3 ;
        sqrt3  = sqrt ( 3.0e+00 ) ;
        iSqrt3 = 1.0e+00 / sqrt3 ;
        nP     = parameters->capacity ;
        for ( i = 0 ; i < nP ; i++ )
        {
            iData = parameters->entries[i] ;
            i0    = Array1D_Item ( basisIndices, i ) ;
            nI    = iData->norbitals ;
            /* . Charge. */
            for ( u = 0 ; u < nI ; u++ ) SymmetricMatrix_Item ( fock, i0 + u, i0 + u ) -= Array1D_Item ( potentials, i ) ;
            /* . p-orbital contributions. */
            if ( doDipoles && ( iData->norbitals >= 4 ) )
            {
                /* . Orbital orders S < PZ < PX < PY. */
                auto Integer _PX = i0 + PX, _PY = i0 + PY, _PZ = i0 + PZ, _S = i0 ;
                auto Real    dA, qA, tX, tXX, tXY, tXZ, tY, tYY, tYZ, tZ, tZZ ;
                dA = iData->ddp[1] ;
                qA = iData->ddp[2] * iData->ddp[2] * 0.5e+00 ;
                /* . Potentials. */
                tX  = Array1D_Item ( potentials, i +   nP ) ;
                tY  = Array1D_Item ( potentials, i + 2*nP ) ;
                tZ  = Array1D_Item ( potentials, i + 3*nP ) ;
                /* . Dipole - sp. */
                SymmetricMatrix_Item ( fock, _PX, _S ) -= dA * tX ;
                SymmetricMatrix_Item ( fock, _PY, _S ) -= dA * tY ;
                SymmetricMatrix_Item ( fock, _PZ, _S ) -= dA * tZ ;
                if ( doQuadrupoles )
                {
                    /* . Potentials. */
                    tXX = Array1D_Item ( potentials, i + 4*nP ) ;
                    tXY = Array1D_Item ( potentials, i + 5*nP ) ;
                    tXZ = Array1D_Item ( potentials, i + 6*nP ) ;
                    tYY = Array1D_Item ( potentials, i + 7*nP ) ;
                    tYZ = Array1D_Item ( potentials, i + 8*nP ) ;
                    tZZ = Array1D_Item ( potentials, i + 9*nP ) ;
                    /* . Quadrupole - pp. */
                    SymmetricMatrix_Item ( fock, _PX, _PX ) -= qA * tXX * 2.0e+00 ;
                    SymmetricMatrix_Item ( fock, _PY, _PX ) -= qA * tXY ;
                    SymmetricMatrix_Item ( fock, _PX, _PZ ) -= qA * tXZ ;
                    SymmetricMatrix_Item ( fock, _PY, _PY ) -= qA * tYY * 2.0e+00 ;
                    SymmetricMatrix_Item ( fock, _PY, _PZ ) -= qA * tYZ ;
                    SymmetricMatrix_Item ( fock, _PZ, _PZ ) -= qA * tZZ * 2.0e+00 ;
                }
                /* . d-orbital contributions. */
                if ( iData->norbitals >= 9 )
                {
                    /* . Orbital orders DZ2 < DXZ < DYZ < DX2Y2 < DXY. */
                    auto Integer _DX2Y2 = i0 + DX2Y2, _DXY = i0 + DXY, _DXZ = i0 + DXZ, _DYZ = i0 + DYZ, _DZ2 = i0 + DZ2 ;
                    auto Real    dB, qB, qC ;
                    dB = iData->ddp[4] ;
                    qB = iData->ddp[3] * iData->ddp[3] * 0.5e+00 ;
                    qC = iData->ddp[5] * iData->ddp[5] * 0.5e+00 ;
                    /* . Dipole - pd. */
                    SymmetricMatrix_Item ( fock, _DXZ,   _PZ ) -= dB * tX ;
                    SymmetricMatrix_Item ( fock, _DX2Y2, _PX ) -= dB * tX ;
                    SymmetricMatrix_Item ( fock, _DXY,   _PY ) -= dB * tX ;
                    SymmetricMatrix_Item ( fock, _DZ2,   _PX ) += dB * tX * iSqrt3 ;
                    SymmetricMatrix_Item ( fock, _DYZ,   _PZ ) -= dB * tY ;
                    SymmetricMatrix_Item ( fock, _DX2Y2, _PY ) += dB * tY ;
                    SymmetricMatrix_Item ( fock, _DXY,   _PX ) -= dB * tY ;
                    SymmetricMatrix_Item ( fock, _DZ2,   _PY ) += dB * tY * iSqrt3 ;
                    SymmetricMatrix_Item ( fock, _DXZ,   _PX ) -= dB * tZ ;
                    SymmetricMatrix_Item ( fock, _DYZ,   _PY ) -= dB * tZ ;
                    SymmetricMatrix_Item ( fock, _DZ2,   _PZ ) -= dB * tZ * 2.0e+00 * iSqrt3 ;
                    if ( doQuadrupoles )
                    {
                        /* . Quadrupole - sd. */
                        SymmetricMatrix_Item ( fock, _DX2Y2, _S ) -= qB * tXX ;
                        SymmetricMatrix_Item ( fock, _DXY  , _S ) -= qB * tXY ;
                        SymmetricMatrix_Item ( fock, _DXZ  , _S ) -= qB * tXZ ;
                        SymmetricMatrix_Item ( fock, _DX2Y2, _S ) += qB * tYY ;
                        SymmetricMatrix_Item ( fock, _DYZ  , _S ) -= qB * tYZ ;
                        SymmetricMatrix_Item ( fock, _DZ2  , _S ) -= qB * tZZ * sqrt3 ;
                        /* . Quadrupole - dd. */
                        SymmetricMatrix_Item ( fock, _DYZ  , _DYZ   ) += qC * tXX * 2.0e+00 ;
                        SymmetricMatrix_Item ( fock, _DX2Y2, _DZ2   ) += qC * tXX * 2.0e+00 * iSqrt3 ;
                        SymmetricMatrix_Item ( fock, _DYZ  , _DXZ   ) -= qC * tXY ;
                        SymmetricMatrix_Item ( fock, _DXY  , _DZ2   ) += qC * tXY * 2.0e+00 * iSqrt3 ;
                        SymmetricMatrix_Item ( fock, _DX2Y2, _DXZ   ) -= qC * tXZ ;
                        SymmetricMatrix_Item ( fock, _DXZ  , _DZ2   ) -= qC * tXZ * iSqrt3 ;
                        SymmetricMatrix_Item ( fock, _DXY  , _DYZ   ) -= qC * tXZ ;
                        SymmetricMatrix_Item ( fock, _DXZ  , _DXZ   ) += qC * tYY * 2.0e+00 ;
                        SymmetricMatrix_Item ( fock, _DX2Y2, _DZ2   ) -= qC * tYY * 2.0e+00 * iSqrt3 ;
                        SymmetricMatrix_Item ( fock, _DXY  , _DXZ   ) -= qC * tYZ ;
                        SymmetricMatrix_Item ( fock, _DX2Y2, _DYZ   ) += qC * tYZ ;
                        SymmetricMatrix_Item ( fock, _DYZ  , _DZ2   ) -= qC * tYZ * iSqrt3 ;
                        SymmetricMatrix_Item ( fock, _DX2Y2, _DX2Y2 ) += qC * tZZ * 2.0e+00 ;
                        SymmetricMatrix_Item ( fock, _DXY  , _DXY   ) += qC * tZZ * 2.0e+00 ;
                        SymmetricMatrix_Item ( fock, _DZ2  , _DZ2   ) -= qC * tZZ * 2.0e+00 ;
                    }
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Bond orders.
! . Bond orders incremented here.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDO_BondOrders ( const IntegerArray1D  *basisIndices ,
                       const SymmetricMatrix *density      ,
                             SymmetricMatrix *bondOrders   )
{
    if ( ( basisIndices != NULL ) &&
         ( bondOrders   != NULL ) &&
         ( density      != NULL ) )
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
                    for ( v = v0 ; v < v1 ; v++ ) sum += pow ( SymmetricMatrix_Item ( density, u, v ), 2 ) ;
                }
                SymmetricMatrix_Item ( bondOrders, i, j ) += sum ;
            }
            /* . Diagonal. */
            for ( u = u0, sum = 0.0e+00 ; u < u1 ; u++ )
            {
                for ( v = u0 ; v < u ; v++ ) sum += ( 2.0e+00 * pow ( SymmetricMatrix_Item ( density, u, v ), 2 ) ) ;
                sum += pow ( SymmetricMatrix_Item ( density, u, u ), 2 ) ;
            }
            SymmetricMatrix_Item ( bondOrders, i, i ) += sum ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Charge restraint W-matrix and core term.
! . This is hugely wasteful for MNDO methods as W is diagonally sparse!
! . However it is done for the moment to simplify the charge restraint code, in particular for those methods,
! . such as DFT with Loewdin charges, for which W is dense.
! . Only basic checking is done.
! . The input W matrix is initialized on entry.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real MNDO_ChargeRestraintMatrix ( const IntegerArray1D  *basisIndices   ,
                                  const RealArray1D     *nuclearCharges ,
                                  const IntegerArray1D  *crIndices      ,
                                  const RealArray1D     *crWeights      ,
                                  const Boolean          isSpin         ,
                                        SymmetricMatrix *W              )
{
    Real core = 0.0e+00 ;
    if ( ( basisIndices   != NULL ) &&
         ( crIndices      != NULL ) &&
         ( crWeights      != NULL ) &&
         ( nuclearCharges != NULL ) &&
         ( W              != NULL ) )
    {
        auto Integer a, i, u, u0, u1 ;
        auto Real    scale, w ;
        if ( isSpin ) scale =  1.0e+00 ;
        else          scale = -1.0e+00 ;
        SymmetricMatrix_Set ( W, 0.0e+00 ) ;
        for ( i = 0 ; i < View1D_Extent ( crIndices ) ; i++ )
        {
            a = Array1D_Item ( crIndices, i ) ;
            w = Array1D_Item ( crWeights, i ) ;
            if ( ! isSpin ) core += ( w * Array1D_Item ( nuclearCharges, a ) ) ;
            u0 = Array1D_Item ( basisIndices, a   ) ;
            u1 = Array1D_Item ( basisIndices, a+1 ) ;
            for ( u = u0 ; u < u1 ; u++ ) SymmetricMatrix_Item ( W, u, u ) += ( scale * w ) ;
        }
    }
    return core ;
}
