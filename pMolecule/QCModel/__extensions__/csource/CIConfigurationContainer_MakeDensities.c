/*==================================================================================================================================
! . Functions to make CI densities.
!=================================================================================================================================*/

# include <math.h>
# include <stdlib.h>

# include "CIConfigurationContainer.h"
# include "NumericalMacros.h"

/*==================================================================================================================================
! . Local procedures declarations.
!=================================================================================================================================*/
static void CIDensity_Diagonal        ( const Integer                nActive ,
                                        const IntegerArray1D        *iAlphas ,
                                        const IntegerArray1D        *iBetas  ,
                                        const Real                   aIaJ    ,
                                              SymmetricMatrix       *onePDMa ,
                                              SymmetricMatrix       *onePDMb ,
                                              DoubleSymmetricMatrix *twoPDM  ) ;
static void CIDensity_OneAlphaOneBeta ( const Integer                nActive ,
                                        const IntegerArray1D        *iAlphas ,
                                        const IntegerArray1D        *iBetas  ,
                                        const IntegerArray1D        *jAlphas ,
                                        const IntegerArray1D        *jBetas  ,
                                        const Real                   aIaJ    ,
                                              DoubleSymmetricMatrix *twoPDM  ) ;
static void CIDensity_OneOrbital      ( const Integer                nActive ,
                                        const IntegerArray1D        *iAlphas ,
                                        const IntegerArray1D        *iBetas  ,
                                        const IntegerArray1D        *jAlphas ,
                                        const Real                   aIaJ    ,
                                              SymmetricMatrix       *onePDM  ,
                                              DoubleSymmetricMatrix *twoPDM  ) ;
static void CIDensity_TwoOrbitals     ( const Integer                nActive ,
                                        const IntegerArray1D        *iAlphas ,
                                        const IntegerArray1D        *jAlphas ,
                                        const Real                   aIaJ    ,
                                              DoubleSymmetricMatrix *twoPDM  ) ;

/*==================================================================================================================================
! . Public procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the total and spin CI densities in the MO basis.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CIConfigurationContainer_MakeDensities ( const CIConfigurationContainer *self      ,
                                              const RealArray1D              *ciVector  ,
                                                    SymmetricMatrix          *onePDMMOa ,
                                                    SymmetricMatrix          *onePDMMOb ,
                                                    DoubleSymmetricMatrix    *twoPDM    ,
                                                    Status                   *status    )
{
    if ( ( self      != NULL ) &&
         ( ciVector  != NULL ) &&
         ( onePDMMOa != NULL ) &&
         ( onePDMMOb != NULL ) &&
         ( twoPDM    != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer          i, j, k, nA, nActive, nAi, nAj, nB ;
        auto IntegerArray1D  *iAlphas, *iBetas, *jAlphas, *jBetas ;
        auto Real             aIaJ ;
        /* . Double loop over configurations. */
        nActive = self->nActive ;
        SymmetricMatrix_Set       ( onePDMMOa, 0.0e+00 ) ;
        SymmetricMatrix_Set       ( onePDMMOb, 0.0e+00 ) ;
        DoubleSymmetricMatrix_Set ( twoPDM   , 0.0e+00 ) ;
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
                    nB += abs ( Array1D_Item ( iBetas,  k ) - Array1D_Item ( jBetas,  k ) ) ;
                }

                /* . Skip if more than two orbitals are different. */
                if ( ( nA + nB ) > 4 ) continue ;

                /* . Get the coefficient factor. */
                aIaJ = 2.0e+00 * Array1D_Item ( ciVector, i ) * Array1D_Item ( ciVector, j ) ;

                /* . Two orbitals different. */
                if ( ( nA + nB ) == 4 )
                {
                    /* . Two beta orbitals. */
                    if      ( nA == 0 ) CIDensity_TwoOrbitals ( nActive, iBetas, jBetas, aIaJ, twoPDM ) ;
                    /* . One alpha and one beta orbital. */
                    else if ( nA == 2 ) CIDensity_OneAlphaOneBeta ( nActive, iAlphas, iBetas, jAlphas, jBetas, aIaJ, twoPDM ) ;
                    /* . Two alpha orbitals. */
                    else                CIDensity_TwoOrbitals ( nActive, iAlphas, jAlphas, aIaJ, twoPDM ) ;
                }
                /* . One alpha orbital different. */
                else if ( nA == 2 )     CIDensity_OneOrbital ( nActive, iAlphas, iBetas,  jAlphas, aIaJ, onePDMMOa, twoPDM ) ;
                /* . One beta orbital different. */
                else if ( nB == 2 )     CIDensity_OneOrbital ( nActive, iBetas,  iAlphas, jBetas,  aIaJ, onePDMMOb, twoPDM ) ;
            }

            /* . Diagonal elements. */
            aIaJ = Array1D_Item ( ciVector, i ) * Array1D_Item ( ciVector, i ) ;
            CIDensity_Diagonal ( nActive, iAlphas, iBetas, aIaJ, onePDMMOa, onePDMMOb, twoPDM ) ;
        }
        /* . Unweight the matrices. */
        SymmetricMatrix_ScaleOffDiagonal ( onePDMMOa, 0.5e+00 ) ;
        SymmetricMatrix_ScaleOffDiagonal ( onePDMMOb, 0.5e+00 ) ;
        DoubleSymmetricMatrix_Unweight   ( twoPDM ) ;
        /* . Convert to total and spin densities. */
        SymmetricMatrix_SumDifference ( onePDMMOa, onePDMMOb, NULL ) ;
    }
}

/*==================================================================================================================================
! . Private procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Diagonal elements.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CIDensity_Diagonal ( const Integer                nActive ,
                                 const IntegerArray1D        *iAlphas ,
                                 const IntegerArray1D        *iBetas  ,
                                 const Real                   aIaJ    ,
                                       SymmetricMatrix       *onePDMa ,
                                       SymmetricMatrix       *onePDMb ,
                                       DoubleSymmetricMatrix *twoPDM  )
{
    Integer  i, j ;

    /* . Loop over active alpha and beta orbitals. */
    for ( i = 0 ; i < nActive ; i++ )
    {
        if ( Array1D_Item ( iAlphas, i ) != 0 )
        {
            SymmetricMatrix_Item ( onePDMa, i, i ) += aIaJ ;

            /* . Alpha/alpha terms. */
            for ( j = 0 ; j < i ; j++ )
            {
                if ( Array1D_Item ( iAlphas, j ) != 0 )
                {
                    DoubleSymmetricMatrix_IncrementItem ( twoPDM, i, i, j, j,  aIaJ, NULL ) ;
                    DoubleSymmetricMatrix_IncrementItem ( twoPDM, i, j, i, j, -aIaJ, NULL ) ;
                }
            }

            /* . Alpha/beta terms. */
            for ( j = 0 ; j < nActive ; j++ )
            {
                if ( Array1D_Item ( iBetas, j ) != 0 ) DoubleSymmetricMatrix_IncrementItem ( twoPDM, i, i, j, j, aIaJ, NULL ) ;
            }
        }
        if ( Array1D_Item ( iBetas, i ) != 0 )
        {
            SymmetricMatrix_Item ( onePDMb, i, i ) += aIaJ ;

            /* . Beta/beta terms. */
            for ( j = 0 ; j < i ; j++ )
            {
                if ( Array1D_Item ( iBetas, j ) != 0 )
                {
                    DoubleSymmetricMatrix_IncrementItem ( twoPDM, i, i, j, j,  aIaJ, NULL ) ;
                    DoubleSymmetricMatrix_IncrementItem ( twoPDM, i, j, i, j, -aIaJ, NULL ) ;
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . States differing by one alpha and one beta orbital.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CIDensity_OneAlphaOneBeta ( const Integer                nActive ,
                                        const IntegerArray1D        *iAlphas ,
                                        const IntegerArray1D        *iBetas  ,
                                        const IntegerArray1D        *jAlphas ,
                                        const IntegerArray1D        *jBetas  ,
                                        const Real                   aIaJ    ,
                                              DoubleSymmetricMatrix *twoPDM  )
{
    Integer  i, j, k, l, n, p ;
    Real    factor ;

    /* . Find the alpha orbitals that differ (i and j with j > i). */
    for ( i = nActive, n = 0   ; n < nActive ; n++ ) { if ( Array1D_Item ( iAlphas, n ) != Array1D_Item ( jAlphas, n ) ) { i = n ; break ; } }
    for ( j = nActive, n = i+1 ; n < nActive ; n++ ) { if ( Array1D_Item ( iAlphas, n ) != Array1D_Item ( jAlphas, n ) ) { j = n ; break ; } }

    /* . Find the beta orbitals that differ (k and l with l > k). */
    for ( k = nActive, n = 0   ; n < nActive ; n++ ) { if ( Array1D_Item ( iBetas , n ) != Array1D_Item ( jBetas , n ) ) { k = n ; break ; } }
    for ( l = nActive, n = k+1 ; n < nActive ; n++ ) { if ( Array1D_Item ( iBetas , n ) != Array1D_Item ( jBetas , n ) ) { l = n ; break ; } }

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
    if ( IsOdd ( p ) ) factor = -aIaJ ;
    else               factor =  aIaJ ;

    /* . Calculate P2. */
    DoubleSymmetricMatrix_IncrementItem ( twoPDM, i, j, k, l, factor, NULL ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . States differing by one alpha or one beta orbital.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CIDensity_OneOrbital ( const Integer                nActive ,
                                   const IntegerArray1D        *iAlphas ,
                                   const IntegerArray1D        *iBetas  ,
                                   const IntegerArray1D        *jAlphas ,
                                   const Real                   aIaJ    ,
                                         SymmetricMatrix       *onePDM  ,
                                         DoubleSymmetricMatrix *twoPDM  )
{
    Integer  i, j, n, p ;
    Real     factor ;

    /* . Find the orbitals that differ (j > i). */
    for ( i = nActive, n = 0   ; n < nActive ; n++ ) { if ( Array1D_Item ( iAlphas, n ) != Array1D_Item ( jAlphas, n ) ) { i = n ; break ; } }
    for ( j = nActive, n = i+1 ; n < nActive ; n++ ) { if ( Array1D_Item ( iAlphas, n ) != Array1D_Item ( jAlphas, n ) ) { j = n ; break ; } }

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
    if ( IsOdd ( p ) ) factor = -aIaJ ;
    else               factor =  aIaJ ;

    /* . Calculate P1 and P2. */
    SymmetricMatrix_Item ( onePDM, j, i ) += factor ;
    for ( n = 0 ; n < nActive ; n++ )
    {
        /* . Common alpha. */
        if ( ( Array1D_Item ( iAlphas, n ) != 0 ) && ( Array1D_Item ( jAlphas, n ) != 0 ) )
        {
            DoubleSymmetricMatrix_IncrementItem ( twoPDM, i, j, n, n,  factor, NULL ) ;
            DoubleSymmetricMatrix_IncrementItem ( twoPDM, i, n, j, n, -factor, NULL ) ;
        }
        /* . Common beta. */
        if ( Array1D_Item ( iBetas , n ) != 0 ) DoubleSymmetricMatrix_IncrementItem ( twoPDM, i, j, n, n,  factor, NULL ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . States differing by two alpha or two beta orbitals.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CIDensity_TwoOrbitals ( const Integer                nActive ,
                                    const IntegerArray1D        *iAlphas ,
                                    const IntegerArray1D        *jAlphas ,
                                    const Real                   aIaJ    ,
                                          DoubleSymmetricMatrix *twoPDM  )
{
    Integer  i, j, k, l, n, p ;
    Real     factor ;

    /* . Find the orbitals that differ (i and j in state 2 with j > i and k and l in state 1 with l > k). */
    for ( i = nActive, n = 0   ; n < nActive ; n++ ) { if ( Array1D_Item ( iAlphas, n ) < Array1D_Item ( jAlphas, n ) ) { i = n ; break ; } }
    for ( j = nActive, n = i+1 ; n < nActive ; n++ ) { if ( Array1D_Item ( iAlphas, n ) < Array1D_Item ( jAlphas, n ) ) { j = n ; break ; } }
    for ( k = nActive, n = 0   ; n < nActive ; n++ ) { if ( Array1D_Item ( iAlphas, n ) > Array1D_Item ( jAlphas, n ) ) { k = n ; break ; } }
    for ( l = nActive, n = k+1 ; n < nActive ; n++ ) { if ( Array1D_Item ( iAlphas, n ) > Array1D_Item ( jAlphas, n ) ) { l = n ; break ; } }

    /* . Check the parity. */
    p = 0 ;
    for ( n = k+1 ; n <= l ; n++ ) p += Array1D_Item ( iAlphas, n ) ;
    for ( n = i+1 ; n <= j ; n++ ) p -= Array1D_Item ( jAlphas, n ) ;
    if ( IsOdd ( p ) ) factor = -aIaJ ;
    else               factor =  aIaJ ;

    /* . Calculate P2. */
    DoubleSymmetricMatrix_IncrementItem ( twoPDM, i, k, j, l,  factor, NULL ) ;
    DoubleSymmetricMatrix_IncrementItem ( twoPDM, i, l, k, j, -factor, NULL ) ;
}
