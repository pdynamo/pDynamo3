/*==================================================================================================================================
! . Functions for solving the CPHF equations.
!=================================================================================================================================*/

# include <math.h>
# include <stdlib.h>

# include "CICPHF.h"
# include "FockConstruction.h"
# include "Real.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate A * B where A is the CPHF TEI matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CICPHF_ApplyCPHFMatrix ( const Integer          n1                   ,
                              const IntegerArray2D  *in1                  ,
                              const Integer          n2                   ,
                              const IntegerArray2D  *in2                  ,
                              const RealArray1D     *aDiagonal            ,
                              const RealArray1D     *b                    ,
                              const RealArray2D     *orbitals             ,
                                    BlockStorage    *twoElectronIntegrals ,
                                    SymmetricMatrix *work1                ,
                                    SymmetricMatrix *work2                ,
                                    RealArray1D     *x                    )
{
    if ( ( in1      != NULL ) &&
         ( in2      != NULL ) &&
         ( b        != NULL ) &&
         ( orbitals != NULL ) &&
         ( work1    != NULL ) &&
         ( work2    != NULL ) &&
         ( x        != NULL ) )
    {
        auto Integer  i, j, n ;
        /* . Initialization. */
        RealArray1D_Set ( x, 0.0e+00 ) ;
        /* . Transform B to the A.O. basis - in work2. */
        CICPHF_Transform ( n2, in2, b, 0, in2, b, orbitals, True, work1, work2 ) ;
        /* . Build Y in the A.O. basis in work1. */
        Fock_MakeFromTEIs ( twoElectronIntegrals, work2, NULL, 1.0e+00, work1, NULL ) ;
        /* . Transform Y to the M.O. basis - in work2.*/
        SymmetricMatrix_Transform ( work1, orbitals, False, work2, NULL ) ;
        /* . Fill X and scale. */
        for ( n = 0 ; n < n1 ; n++ )
        {
            i = Array2D_Item ( in1, n, 0 ) ;
            j = Array2D_Item ( in1, n, 1 ) ;
            Array1D_Item ( x, n ) = SymmetricMatrix_Item ( work2, j, i ) ;
        }
        RealArray1D_Scale ( x, 4.0e+00 ) ;
        /* . Add in the diagonal terms. */
        if ( aDiagonal != NULL )
        {
            for ( i = 0 ; i < n1 ; i++ ) Array1D_Item ( x, i ) += ( Array1D_Item ( aDiagonal, i ) * Array1D_Item ( b, i ) ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the vectors required for solution of the CPHF equations.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define OccupancyTolerance         1.0e-06
# define OrbitalDegeneracyTolerance 1.0e-06
# define PreconditionerTolerance    1.0e-06
# define ZeroTolerance              1.0e-12
void CICPHF_CalculateCPHFVectors ( const Integer                nActive                   ,
                                   const Integer                nCore                     ,
                                   const Integer                nOrbitals                 ,
                                         BlockStorage          *twoElectronIntegrals      ,
                                   const DoubleSymmetricMatrix *twoPDM                    ,
                                   const RealArray1D           *energies                  ,
                                   const RealArray1D           *occupancies               ,
                                   const RealArray2D           *orbitals                  ,
                                   const RealArrayND           *moTEI234                  ,
                                   const SymmetricMatrix       *fCore                     ,
                                   const SymmetricMatrix       *onePDM                    ,
                                   const SymmetricMatrix       *onePDMMO                  ,
                                         SymmetricMatrix       *work1                     ,
                                         SymmetricMatrix       *work2                     ,
                                         Integer               *numberDegenerateRedundant ,
                                         Integer               *numberNonRedundant        ,
                                         Integer               *numberRedundant           ,
                                         IntegerArray2D        *indicesNR                 ,
                                         IntegerArray2D        *indicesR                  ,
                                         RealArray1D           *aDiagonal                 ,
                                         RealArray1D           *qNR                       ,
                                         RealArray1D           *qR                        ,
                                         RealArray1D           *preconditioner            ,
                                         Status                *status                    )
{
    if ( ( twoElectronIntegrals      != NULL ) &&
         ( twoPDM                    != NULL ) &&
         ( energies                  != NULL ) &&
         ( occupancies               != NULL ) &&
         ( orbitals                  != NULL ) &&
         ( moTEI234                  != NULL ) &&
         ( fCore                     != NULL ) &&
         ( onePDM                    != NULL ) &&
         ( onePDMMO                  != NULL ) &&
         ( work1                     != NULL ) &&
         ( work2                     != NULL ) &&
         ( numberDegenerateRedundant != NULL ) &&
         ( numberNonRedundant        != NULL ) &&
         ( numberRedundant           != NULL ) &&
         ( indicesNR                 != NULL ) &&
         ( indicesR                  != NULL ) &&
         ( aDiagonal                 != NULL ) &&
         ( qNR                       != NULL ) &&
         ( qR                        != NULL ) &&
         ( preconditioner            != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer          a, i, nCoreAndActive, nDegenerateRedundant, nNonRedundant, nRedundant, p, q, r ;
        auto Real             energyDifference, f, iEnergy, iOccupancy, occupancyDifference, pEnergy, pOccupancy ;
        auto RealArray2D     *twoXY ;
        auto SymmetricMatrix *fTransformed, *gGamma ;
        /* . Initialization. */
        nCoreAndActive       = nCore + nActive ;
        nDegenerateRedundant = 0 ;
        nNonRedundant        = 0 ;
        nRedundant           = 0 ;
        RealArray1D_Set ( qNR, 0.0e+00 ) ;
        RealArray1D_Set ( qR , 0.0e+00 ) ;
        /* . Calculate gGamma - use fTransformed as scratch. */
        fTransformed = work1 ;
        gGamma       = work2 ;
        if ( nCore > 0 )
        {
            Fock_MakeFromTEIs         ( twoElectronIntegrals, onePDM, NULL, 1.0e+00, fTransformed, NULL ) ;
            SymmetricMatrix_Transform ( fTransformed, orbitals, False, gGamma, NULL ) ;
            SymmetricMatrix_Scale     ( gGamma, 2.0e+00 ) ;
        }
        /* . Transform fCore to the M.O. basis. */
        SymmetricMatrix_Transform ( fCore, orbitals, False, fTransformed, NULL ) ;
        /* . Gamma terms. */
        twoXY = CICPHF_CalculateKPA ( nActive, nOrbitals, twoPDM, orbitals, moTEI234, status ) ;
        if ( ! Status_IsOK ( status ) ) goto FinishUp ;
        /* . Fill the elements in order - non-redundant and redundant at the same time. */
        /* . Core-active. */
        for ( i = 0 ; i < nCore ; i++ )
        {
            iEnergy    = Array1D_Item ( energies   , i ) ;
            iOccupancy = Array1D_Item ( occupancies, i ) ;
            for ( p = nCore ; p < nCoreAndActive ; p++ )
            {
                energyDifference    = Array1D_Item ( energies, p ) - iEnergy ;
                occupancyDifference = 0.5e+00 * ( iOccupancy - Array1D_Item ( occupancies, p ) ) ;
                for ( f = 0.0e+00, r = 0 ; r < nActive ; r++ ) f += ( SymmetricMatrix_GetItem ( onePDMMO, p - nCore, r, NULL ) * SymmetricMatrix_Item ( fTransformed, r + nCore, i ) ) ;
                f += ( Array2D_Item ( twoXY, i, p - nCore ) - 2.0e+00 * SymmetricMatrix_Item ( fTransformed, p, i ) - SymmetricMatrix_Item ( gGamma, p, i ) ) ;
                if ( fabs ( occupancyDifference ) > OccupancyTolerance )
                {
                    Array2D_Item ( indicesNR, nNonRedundant, 0 ) = i ;
                    Array2D_Item ( indicesNR, nNonRedundant, 1 ) = p ;
                    Array1D_Item     ( aDiagonal, nNonRedundant    ) = energyDifference / occupancyDifference ;
                    Array1D_Item     ( qNR      , nNonRedundant    ) = f                / occupancyDifference ;
                    nNonRedundant += 1 ;
                }
                else if ( fabs ( f ) > ZeroTolerance )
                {
                    if ( fabs ( energyDifference ) > OrbitalDegeneracyTolerance )
                    {
                        Array2D_Item ( indicesR, nRedundant, 0 ) = i ;
                        Array2D_Item ( indicesR, nRedundant, 1 ) = p ;
                        Array1D_Item     ( qR      , nRedundant    ) = f / energyDifference ;
                        nRedundant += 1 ;
                    }
                    else nDegenerateRedundant += 1 ;
                }
            }
        }
        /* . Core-virtual (only non-redundant). */
        for ( i = 0 ; i < nCore ; i++ )
        {
            iEnergy    = Array1D_Item ( energies   , i ) ;
            iOccupancy = Array1D_Item ( occupancies, i ) ;
            for ( a = nCoreAndActive ; a < nOrbitals ; a++ )
            {
                energyDifference    = Array1D_Item ( energies, a ) - iEnergy ;
                occupancyDifference = 0.5e+00 * ( iOccupancy - Array1D_Item ( occupancies, a ) ) ;
                f = - 2.0e+00 * SymmetricMatrix_Item ( fTransformed, a, i ) - SymmetricMatrix_Item ( gGamma, a, i ) ;
                Array2D_Item ( indicesNR, nNonRedundant, 0 ) = i ;
                Array2D_Item ( indicesNR, nNonRedundant, 1 ) = a ;
                Array1D_Item     ( aDiagonal, nNonRedundant    ) = energyDifference / occupancyDifference ;
                Array1D_Item     ( qNR      , nNonRedundant    ) = f                / occupancyDifference ;
                nNonRedundant += 1 ;
            }
        }
        /* . Active-active. */
        for ( p = nCore ; p < nCoreAndActive ; p++ )
        {
            pEnergy    = Array1D_Item ( energies   , p ) ;
            pOccupancy = Array1D_Item ( occupancies, p ) ;
            for ( q = ( p + 1 ) ; q < nCoreAndActive ; q++ )
            {
                energyDifference    = Array1D_Item ( energies, q ) - pEnergy ;
                occupancyDifference = 0.5e+00 * ( pOccupancy - Array1D_Item ( occupancies, q ) ) ;
                for ( f = 0.0e+00, r = 0 ; r < nActive ; r++ ) f += ( SymmetricMatrix_GetItem ( onePDMMO, q - nCore, r, NULL ) * SymmetricMatrix_GetItem ( fTransformed, p, r + nCore, NULL ) -
                                                                      SymmetricMatrix_GetItem ( onePDMMO, p - nCore, r, NULL ) * SymmetricMatrix_GetItem ( fTransformed, q, r + nCore, NULL ) ) ;
                f += ( Array2D_Item ( twoXY, p, q - nCore ) - Array2D_Item ( twoXY, q, p - nCore ) ) ;
                if ( fabs ( occupancyDifference ) > OccupancyTolerance )
                {
                    Array2D_Item ( indicesNR, nNonRedundant, 0 ) = p ;
                    Array2D_Item ( indicesNR, nNonRedundant, 1 ) = q ;
                    Array1D_Item     ( aDiagonal, nNonRedundant    ) = energyDifference / occupancyDifference ;
                    Array1D_Item     ( qNR      , nNonRedundant    ) = f                / occupancyDifference ;
                    nNonRedundant += 1 ;
                }
                else if ( fabs ( f ) > ZeroTolerance )
                {
                    if ( fabs ( energyDifference ) > OrbitalDegeneracyTolerance )
                    {
                        Array2D_Item ( indicesR, nRedundant, 0 ) = p ;
                        Array2D_Item ( indicesR, nRedundant, 1 ) = q ;
                        Array1D_Item     ( qR      , nRedundant    ) = f / energyDifference ;
                        nRedundant += 1 ;
                    }
                    else nDegenerateRedundant += 1 ;
                }
            }
        }
        /* . Active-virtual. */
        for ( p = nCore ; p < nCoreAndActive ; p++ )
        {
            pEnergy    = Array1D_Item ( energies   , p ) ;
            pOccupancy = Array1D_Item ( occupancies, p ) ;
            for ( a = nCoreAndActive ; a < nOrbitals ; a++ )
            {
                energyDifference    = Array1D_Item ( energies, a ) - pEnergy ;
                occupancyDifference = 0.5e+00 * ( pOccupancy - Array1D_Item ( occupancies, a ) ) ;
                for ( f = 0.0e+00, r = 0 ; r < nActive ; r++ ) f -= ( SymmetricMatrix_GetItem ( onePDMMO, p - nCore, r, NULL ) * SymmetricMatrix_Item ( fTransformed, a, r + nCore ) ) ;
                f -= Array2D_Item ( twoXY, a, p - nCore ) ;
                if ( fabs ( occupancyDifference ) > OccupancyTolerance )
                {
                    Array2D_Item ( indicesNR, nNonRedundant, 0 ) = p ;
                    Array2D_Item ( indicesNR, nNonRedundant, 1 ) = a ;
                    Array1D_Item     ( aDiagonal, nNonRedundant    ) = energyDifference / occupancyDifference ;
                    Array1D_Item     ( qNR      , nNonRedundant    ) = f                / occupancyDifference ;
                    nNonRedundant += 1 ;
                }
                else if ( fabs ( f ) > ZeroTolerance )
                {
                    if ( fabs ( energyDifference ) > OrbitalDegeneracyTolerance )
                    {
                        Array2D_Item ( indicesR, nRedundant, 0 ) = p ;
                        Array2D_Item ( indicesR, nRedundant, 1 ) = a ;
                        Array1D_Item     ( qR      , nRedundant    ) = f / energyDifference ;
                        nRedundant += 1 ;
                    }
                    else nDegenerateRedundant += 1 ;
                }
            }
        }
        /* . Remove redundant terms from qNR by using the redundant A matrix and qR. */
        if ( nRedundant > 0 )
        {
            CICPHF_ApplyCPHFMatrix ( nNonRedundant, indicesNR, nRedundant, indicesR, NULL, qR, orbitals, twoElectronIntegrals, work1, work2, preconditioner ) ;
            RealArray1D_Add ( qNR, -1.0e+00, preconditioner, NULL ) ;
        }
        /* . Determine the preconditioner.*/
        /* . Should there be a square root here? */
        for ( i = 0 ; i < nNonRedundant ; i++ )
        {
            f = Array1D_Item ( aDiagonal, i ) ;
            if ( fabs ( f ) > PreconditionerTolerance ) Array1D_Item ( preconditioner, i ) = 1.0e+00 / sqrt ( fabs ( f ) ) ;
            else                                        Array1D_Item ( preconditioner, i ) = 1.0e+00 / PreconditionerTolerance ;
        }
        /* . Finish up. */
FinishUp:
        (*numberDegenerateRedundant) = nDegenerateRedundant ;
        (*numberNonRedundant       ) = nNonRedundant        ;
        (*numberRedundant          ) = nRedundant           ;
        RealArray2D_Deallocate ( &twoXY ) ;
    }
}
# undef OccupancyTolerance
# undef OrbitalDegeneracyTolerance
# undef PreconditionerTolerance
# undef ZeroTolerance

/*----------------------------------------------------------------------------------------------------------------------------------
! . Create a quantity of the form Kpa = Sum_qrs Gamma_pqrs TEI234_aqrs: a runs over all AOs.
! . This is then transformed to the MO basis and scaled by 2.
!---------------------------------------------------------------------------------------------------------------------------------*/
RealArray2D *CICPHF_CalculateKPA ( const Integer                nActive  ,
                                   const Integer                nBasis   ,
                                   const DoubleSymmetricMatrix *twoPDM   ,
                                   const RealArray2D           *orbitals ,
                                   const RealArrayND           *moTEI234 ,
                                         Status                *status   )
{
    RealArray2D *kpaMO = NULL ;
    if ( ( twoPDM   != NULL ) &&
         ( orbitals != NULL ) &&
         ( moTEI234 != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto RealArray2D *kpa ;
        kpa   = RealArray2D_AllocateWithExtents ( nBasis, nActive, status ) ;
        kpaMO = RealArray2D_AllocateWithExtents ( nBasis, nActive, status ) ;
        if ( ( kpa != NULL ) && ( kpaMO != NULL ) )
        {
            auto Integer      i, p, pq, q, r, rs, s ;
            auto Real         t ;
            RealArray2D_Set ( kpa, 0.0e+00 ) ;
            for ( p = 0 ; p < nActive ; p++ )
            {
                for ( q = 0 ; q < nActive ; q++, pq++ )
                {
                    for ( r = rs = 0 ; r < nActive ; r++ )
                    {
                        for ( s = 0 ; s <= r ; s++, rs++ )
                        {
                            if ( r == s ) t = DoubleSymmetricMatrix_GetItem ( twoPDM, p, q, r, r, NULL ) ;
                            else          t = DoubleSymmetricMatrix_GetItem ( twoPDM, p, q, r, s, NULL ) + DoubleSymmetricMatrix_GetItem ( twoPDM, p, q, s, r, NULL ) ;
                            /* . This is an increment so ultimately should use slices. */
                            for ( i = 0 ; i < nBasis ; i++ ) Array2D_Item ( kpa, i, p ) += t * ArrayND_Item3D ( moTEI234, i, q, rs ) ;
                        }
                    }
                }
            }
            /* . Transform to the full MO basis. */
            RealArray2D_MatrixMultiply ( True, False, 1.0e+00, orbitals, kpa, 0.0e+00, kpaMO, NULL ) ;
            RealArray2D_Scale ( kpaMO, 2.0e+00 ) ;
        }
    }
    return kpaMO ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Convert packed vectors indexed by M.O.s to the A.O. basis.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . Note X1 and X2 are not symmetric although the output vector, Z, will be.
! . No index pair occurs more than once and there are no diagonal pairs.
! . The process of symmetrizing introduces a factor of 2 which should be corrected for.
*/
void CICPHF_Transform ( const Integer          n1       ,
                        const IntegerArray2D  *in1      ,
                        const RealArray1D     *x1       ,
                        const Integer          n2       ,
                        const IntegerArray2D  *in2      ,
                        const RealArray1D     *x2       ,
                        const RealArray2D     *orbitals ,
                        const Boolean          doScale  ,
                              SymmetricMatrix *work     ,
                              SymmetricMatrix *z        )
{
    auto Integer  i, j, n ;
    SymmetricMatrix_Set ( work, 0.0e+00 ) ;
    SymmetricMatrix_Set ( z   , 0.0e+00 ) ;
    /* . Unpack the elements (first index always less than second). */
    if ( n1 > 0 )
    {
        for ( n = 0 ; n < n1 ; n++ )
        {
            i = Array2D_Item ( in1, n, 0 ) ;
            j = Array2D_Item ( in1, n, 1 ) ;
            SymmetricMatrix_Item ( work, j, i ) = Array1D_Item ( x1, n ) ;
        }
    }
    if ( n2 > 0 )
    {
        for ( n = 0 ; n < n2 ; n++ )
        {
            i = Array2D_Item ( in2, n, 0 ) ;
            j = Array2D_Item ( in2, n, 1 ) ;
            SymmetricMatrix_Item ( work, j, i ) = Array1D_Item ( x2, n ) ;
        }
    }
    /* . Transform the Z-matrix to the A.O. basis. */
    SymmetricMatrix_Transform ( work, orbitals, True, z, NULL ) ;
    /* . Scale if necessary. */
    if ( doScale ) SymmetricMatrix_Scale ( z, 0.5e+00 ) ;
}
