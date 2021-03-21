/*==================================================================================================================================
! . Functions to make CI transition dipoles.
!=================================================================================================================================*/

# include <math.h>
# include <stdlib.h>

# include "CIConfigurationContainer.h"
# include "NumericalMacros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real CITD_OneOrbital ( const Integer          nActive ,
                              const IntegerArray1D  *iAlphas ,
                              const IntegerArray1D  *jAlphas ,
                              const SymmetricMatrix *tdMOs   ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the TD matrix between configurations in tdMatrix - dense only for the moment.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CIConfigurationContainer_TransitionDipoles ( const CIConfigurationContainer *self     ,
                                                        SymmetricMatrix          *tdMOs    ,
                                                        SymmetricMatrix          *tdMatrix )
{
    if ( ( self     != NULL ) &&
         ( tdMOs    != NULL ) &&
         ( tdMatrix != NULL ) )
    {
        auto Integer          i, j, k, nA, nActive, nAi, nAj, nB ;
        auto IntegerArray1D  *iAlphas, *iBetas, *jAlphas, *jBetas ;
        /* . Initialization. */
        SymmetricMatrix_Set ( tdMatrix, 0.0e+00 ) ;
        /* . Double loop over configurations. */
        nActive = self->nActive ;
        for ( i = 1 ; i < self->nConfigurations ; i++ )
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
                /* . One alpha orbital different. */
                     if ( ( nA == 2 ) && ( nB == 0 ) ) SymmetricMatrix_Item ( tdMatrix, i, j ) = CITD_OneOrbital ( nActive, iAlphas, jAlphas, tdMOs ) ;
                /* . One beta orbital different. */
                else if ( ( nA == 0 ) && ( nB == 2 ) ) SymmetricMatrix_Item ( tdMatrix, i, j ) = CITD_OneOrbital ( nActive, iBetas , jBetas , tdMOs ) ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . TD matrix elements for states differing by one orbital.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real CITD_OneOrbital ( const Integer  nActive, const IntegerArray1D  *iAlphas, const IntegerArray1D  *jAlphas, const SymmetricMatrix *tdMOs )
{
    Integer  i, j, n, p ;
    Real     hIJ ;
    /* . Find the orbitals that differ (j > i). */
    for ( i = -1, n = 0   ; n < nActive ; n++ ) { if ( Array1D_Item ( iAlphas, n ) != Array1D_Item ( jAlphas, n ) ) { i = n ; break ; } }
    for ( j = -1, n = i+1 ; n < nActive ; n++ ) { if ( Array1D_Item ( iAlphas, n ) != Array1D_Item ( jAlphas, n ) ) { j = n ; break ; } }
    /* . Calculate the matrix element. */
    hIJ = SymmetricMatrix_Item ( tdMOs, j, i ) ;
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
    if ( IsOdd ( p ) ) hIJ *= -1.0e+00 ;
    return hIJ ;
}
