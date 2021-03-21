/*==================================================================================================================================
! . Functions to build CI matrices.
!=================================================================================================================================*/

# include <math.h>
# include <stdlib.h>

# include "CIConfigurationContainer.h"
# include "NumericalMacros.h"

/*==================================================================================================================================
! . Local procedures declarations.
!=================================================================================================================================*/
static Real CIMatrix_Diagonal        ( const Integer                nActive ,
                                       const IntegerArray1D        *iAlphas ,
                                       const IntegerArray1D        *iBetas  ,
                                       const SymmetricMatrix       *fCore   ,
                                       const DoubleSymmetricMatrix *moTEIs  ) ;
static Real CIMatrix_OneAlphaOneBeta ( const Integer                nActive ,
                                       const IntegerArray1D        *iAlphas ,
                                       const IntegerArray1D        *iBetas  ,
                                       const IntegerArray1D        *jAlphas ,
                                       const IntegerArray1D        *jBetas  ,
                                       const DoubleSymmetricMatrix *moTEIs  ) ;
static Real CIMatrix_OneOrbital      ( const Integer                nActive ,
                                       const IntegerArray1D        *iAlphas ,
                                       const IntegerArray1D        *iBetas  ,
                                       const IntegerArray1D        *jAlphas ,
                                       const SymmetricMatrix       *fCore   ,
                                       const DoubleSymmetricMatrix *moTEIs  ) ;
static Real CIMatrix_TwoOrbitals     ( const Integer                nActive ,
                                       const IntegerArray1D        *iAlphas ,
                                       const IntegerArray1D        *jAlphas ,
                                       const DoubleSymmetricMatrix *moTEIs  ) ;

/*==================================================================================================================================
! . Public procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the CI matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Check matrix dimensions somewhere? */
void CIConfigurationContainer_MakeCIMatrix ( const CIConfigurationContainer *self           ,
                                             const SymmetricMatrix          *fCoreMO        ,
                                             const DoubleSymmetricMatrix    *moTEIs         ,
                                                   SymmetricMatrix          *ciMatrixFull   ,
                                                   SparseSymmetricMatrix    *ciMatrixSparse ,
                                                   Status                   *status         )
{
    auto Boolean doFull   = ( ciMatrixFull   != NULL ) ,
                 doSparse = ( ciMatrixSparse != NULL ) ;
    if ( ( self   != NULL     ) &&
         ( doFull || doSparse ) &&
         Status_IsOK ( status ) )
    {
        auto Integer          i, j, k, nA, nActive, nAi, nAj, nB ;
        auto IntegerArray1D  *iAlphas, *iBetas, *jAlphas, *jBetas ;
        auto Real             value = 0.0e+00 ;
        /* . Initialization. */
        if ( doFull   ) SymmetricMatrix_Set ( ciMatrixFull, 0.0e+00 ) ;
        if ( doSparse ) SparseSymmetricMatrix_Clear ( ciMatrixSparse ) ;
        /* . Double loop over configurations. */
        nActive = self->nActive ;
        for ( i = 0 ; i < self->nConfigurations ; i++ )
        {
            nAi     = self->configurations[i].nAlphas ;
            iAlphas = self->configurations[i].alphas  ;
            iBetas  = self->configurations[i].betas   ;
            for ( j = 0 ; j < i ; j++ )
            {
                nAj     = self->configurations[j].nAlphas ;
                jAlphas = self->configurations[j].alphas  ;
                jBetas  = self->configurations[j].betas   ;
                /* . Skip if there are different numbers of alpha orbitals in the two configurations. */
                if ( nAi != nAj ) continue ;
                /* . Find the differences in the numbers of alpha and beta orbitals including positional information. */
                for ( k = nA = nB = 0 ; k < nActive ; k++ )
                {
                    nA += abs ( Array1D_Item ( iAlphas, k ) - Array1D_Item ( jAlphas, k ) ) ;
                    nB += abs ( Array1D_Item ( iBetas , k ) - Array1D_Item ( jBetas , k ) ) ;
                }
                /* . Skip if more than two orbitals are different. */
                if ( ( nA + nB ) > 4 ) continue ;
                /* . Two orbitals different. */
                if ( ( nA + nB ) == 4 )
                {
                    /* . Two beta orbitals. */
                    if      ( nA == 0 ) value = CIMatrix_TwoOrbitals ( nActive, iBetas, jBetas, moTEIs ) ;
                    /* . One alpha and one beta orbital. */
                    else if ( nA == 2 ) value = CIMatrix_OneAlphaOneBeta ( nActive, iAlphas, iBetas, jAlphas, jBetas, moTEIs ) ;
                    /* . Two alpha orbitals. */
                    else                value = CIMatrix_TwoOrbitals ( nActive, iAlphas, jAlphas, moTEIs ) ;
                }
                /* . One alpha orbital different. */
                else if ( nA == 2 )     value = CIMatrix_OneOrbital ( nActive, iAlphas, iBetas , jAlphas, fCoreMO, moTEIs ) ;
                /* . One beta orbital different. */
                else if ( nB == 2 )     value = CIMatrix_OneOrbital ( nActive, iBetas , iAlphas, jBetas , fCoreMO, moTEIs ) ;
                /* . Save the value. */
                if ( doFull   ) SymmetricMatrix_Item ( ciMatrixFull, i, j ) = value ;
                if ( doSparse ) SparseSymmetricMatrix_AppendItem ( ciMatrixSparse, i, j, value, NULL ) ;
            }
            /* . Diagonal elements. */
            value = CIMatrix_Diagonal ( nActive, iAlphas, iBetas, fCoreMO, moTEIs ) ;
            if ( doFull   ) SymmetricMatrix_Item ( ciMatrixFull, i, i ) = value ;
            if ( doSparse ) SparseSymmetricMatrix_AppendItem ( ciMatrixSparse, i, i, value, NULL ) ;
        }
        /* . Finalization. */
        if ( doSparse ) SparseSymmetricMatrix_Canonicalize ( ciMatrixSparse, NULL ) ;
    }
}

/*==================================================================================================================================
! . CI matrix element procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Diagonal elements.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real CIMatrix_Diagonal ( const Integer                nActive ,
                                const IntegerArray1D        *iAlphas ,
                                const IntegerArray1D        *iBetas  ,
                                const SymmetricMatrix       *fCore   ,
                                const DoubleSymmetricMatrix *moTEIs  )
{
    Integer  i, j ;
    Real     hij = 0.0e+00 ;

    /* . Loop over active orbitals. */
    for ( i = 0 ; i < nActive ; i++ )
    {
        /* . Alpha orbital. */
        if ( Array1D_Item ( iAlphas, i ) != 0 )
        {
            hij += SymmetricMatrix_Item ( fCore, i, i ) ;

            /* . Alpha/alpha terms. */
            for ( j = 0 ; j < i ; j++ )
            {
                if ( Array1D_Item ( iAlphas, j ) != 0 ) hij += ( DoubleSymmetricMatrix_GetItem ( moTEIs, i, i, j, j, NULL ) - DoubleSymmetricMatrix_GetItem ( moTEIs, i, j, i, j, NULL ) ) ;
            }

            /* . Alpha/beta terms. */
            for ( j = 0 ; j < nActive ; j++ )
            {
                if ( Array1D_Item ( iBetas, j ) != 0 ) hij += DoubleSymmetricMatrix_GetItem ( moTEIs, i, i, j, j, NULL ) ;
            }
        }

        /* . Beta orbital. */
        if ( Array1D_Item ( iBetas, i ) != 0 )
        {
            hij += SymmetricMatrix_Item ( fCore, i, i ) ;

            /* . Beta/beta terms. */
            for ( j = 0 ; j < i ; j++ )
            {
                if ( Array1D_Item ( iBetas, j ) != 0 ) hij += ( DoubleSymmetricMatrix_GetItem ( moTEIs, i, i, j, j, NULL ) - DoubleSymmetricMatrix_GetItem ( moTEIs, i, j, i, j, NULL ) ) ;
            }
        }
    }
    return hij ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . States differing by one alpha and one beta orbital.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real CIMatrix_OneAlphaOneBeta ( const Integer                nActive ,
                                       const IntegerArray1D        *iAlphas ,
                                       const IntegerArray1D        *iBetas  ,
                                       const IntegerArray1D        *jAlphas ,
                                       const IntegerArray1D        *jBetas  ,
                                       const DoubleSymmetricMatrix *moTEIs  )
{
    Integer  i, j, k, l, n, p ;
    Real     hij = 0.0e+00 ;

    /* . Find the alpha orbitals that differ (i and j with j > i). */
    for ( i = nActive, n = 0   ; n < nActive ; n++ ) { if ( Array1D_Item ( iAlphas, n ) != Array1D_Item ( jAlphas, n ) ) { i = n ; break ; } }
    for ( j = nActive, n = i+1 ; n < nActive ; n++ ) { if ( Array1D_Item ( iAlphas, n ) != Array1D_Item ( jAlphas, n ) ) { j = n ; break ; } }

    /* . Find the beta orbitals that differ (k and l with l > k). */
    for ( k = nActive, n = 0   ; n < nActive ; n++ ) { if ( Array1D_Item ( iBetas , n ) != Array1D_Item ( jBetas , n ) ) { k = n ; break ; } }
    for ( l = nActive, n = k+1 ; n < nActive ; n++ ) { if ( Array1D_Item ( iBetas , n ) != Array1D_Item ( jBetas , n ) ) { l = n ; break ; } }

    /* . Calculate the matrix element. */
    hij = DoubleSymmetricMatrix_GetItem ( moTEIs, i, j, k, l, NULL ) ;

    /* . Check the parity. */
    p = 0 ;
    /* . i in state 2, j in state 1. */
    if ( Array1D_Item ( iAlphas, i ) == 0 )
    {
        for ( n = 0 ; n <= j ; n++ ) p += Array1D_Item ( iAlphas, n ) ;
        for ( n = 0 ; n <= i ; n++ ) p -= Array1D_Item ( jAlphas, n ) ;
    }
    /* . i in state 1, j in state 2. */
    else
    {
        for ( n = 0 ; n <= j ; n++ ) p += Array1D_Item ( jAlphas, n ) ;
        for ( n = 0 ; n <= i ; n++ ) p -= Array1D_Item ( iAlphas, n ) ;
    }
    /* . k in state 2, l in state 1. */
    if ( Array1D_Item ( iBetas, k ) == 0 )
    {
        for ( n = 0 ; n <= l ; n++ ) p += Array1D_Item ( iBetas , n ) ;
        for ( n = 0 ; n <= k ; n++ ) p -= Array1D_Item ( jBetas , n ) ;
    }
    /* . k in state 1, l in state 2. */
    else
    {
        for ( n = 0 ; n <= l ; n++ ) p += Array1D_Item ( jBetas , n ) ;
        for ( n = 0 ; n <= k ; n++ ) p -= Array1D_Item ( iBetas , n ) ;
    }
    if ( IsOdd ( p ) ) hij *= -1.0e+00 ;

    return hij ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . States differing by one alpha or one beta orbital.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real CIMatrix_OneOrbital ( const Integer                nActive ,
                                  const IntegerArray1D        *iAlphas ,
                                  const IntegerArray1D        *iBetas  ,
                                  const IntegerArray1D        *jAlphas ,
                                  const SymmetricMatrix       *fCore   ,
                                  const DoubleSymmetricMatrix *moTEIs  )
{
    Integer  i, j, n, p ;
    Real     hij ;

    /* . Find the orbitals that differ (j > i). */
    for ( i = nActive, n = 0   ; n < nActive ; n++ ) { if ( Array1D_Item ( iAlphas, n ) != Array1D_Item ( jAlphas, n ) ) { i = n ; break ; } }
    for ( j = nActive, n = i+1 ; n < nActive ; n++ ) { if ( Array1D_Item ( iAlphas, n ) != Array1D_Item ( jAlphas, n ) ) { j = n ; break ; } }

    /* . Calculate the matrix element. */
    hij = SymmetricMatrix_Item ( fCore, j, i ) ;
    for ( n = 0 ; n < nActive ; n++ )
    {
        /* . Common alpha. */
        if ( ( Array1D_Item ( iAlphas, n ) != 0 ) && ( Array1D_Item ( jAlphas, n ) != 0 ) ) hij += ( DoubleSymmetricMatrix_GetItem ( moTEIs, i, j, n, n, NULL ) - DoubleSymmetricMatrix_GetItem ( moTEIs, i, n, j, n, NULL ) ) ;
        /* . Common beta. */
        if (   Array1D_Item ( iBetas , n ) != 0                                                  ) hij +=   DoubleSymmetricMatrix_GetItem ( moTEIs, i, j, n, n, NULL ) ;
    }

    /* . Check the parity. */
    p = 0 ;
    /* . i in state 2, j in state 1. */
    if ( Array1D_Item ( iAlphas, i ) == 0 )
    {
        for ( n = 0 ; n <= j ; n++ ) p += Array1D_Item ( iAlphas, n ) ;
        for ( n = 0 ; n <= i ; n++ ) p -= Array1D_Item ( jAlphas, n ) ;
    }
    /* . i in state 1, j in state 2. */
    else
    {
        for ( n = 0 ; n <= j ; n++ ) p += Array1D_Item ( jAlphas, n ) ;
        for ( n = 0 ; n <= i ; n++ ) p -= Array1D_Item ( iAlphas, n ) ;
    }
    if ( IsOdd ( p ) ) hij *= -1.0e+00 ;

    return hij ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . States differing by two alpha or two beta orbitals.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real CIMatrix_TwoOrbitals ( const Integer                nActive ,
                                   const IntegerArray1D        *iAlphas ,
                                   const IntegerArray1D        *jAlphas ,
                                   const DoubleSymmetricMatrix *moTEIs  )
{
    Integer  i, j, k, l, n, p ;
    Real     hij = 0.0e+00 ;

    /* . Find the orbitals that differ (i and j in state 2 with j > i and k and l in state 1 with l > k). */
    for ( i = nActive, n = 0   ; n < nActive ; n++ ) { if ( Array1D_Item ( iAlphas, n ) < Array1D_Item ( jAlphas, n ) ) { i = n ; break ; } }
    for ( j = nActive, n = i+1 ; n < nActive ; n++ ) { if ( Array1D_Item ( iAlphas, n ) < Array1D_Item ( jAlphas, n ) ) { j = n ; break ; } }
    for ( k = nActive, n = 0   ; n < nActive ; n++ ) { if ( Array1D_Item ( iAlphas, n ) > Array1D_Item ( jAlphas, n ) ) { k = n ; break ; } }
    for ( l = nActive, n = k+1 ; n < nActive ; n++ ) { if ( Array1D_Item ( iAlphas, n ) > Array1D_Item ( jAlphas, n ) ) { l = n ; break ; } }

    /* . Calculate the matrix element. */
    hij = DoubleSymmetricMatrix_GetItem ( moTEIs, i, k, j, l, NULL ) - DoubleSymmetricMatrix_GetItem ( moTEIs, i, l, k, j, NULL ) ;

    /* . Check the parity. */
    p = 0 ;
    for ( n = k+1 ; n <= l ; n++ ) p += Array1D_Item ( iAlphas, n ) ;
    for ( n = i+1 ; n <= j ; n++ ) p -= Array1D_Item ( jAlphas, n ) ;
    if ( IsOdd ( p ) ) hij *= -1.0e+00 ;

    return hij ;
}
