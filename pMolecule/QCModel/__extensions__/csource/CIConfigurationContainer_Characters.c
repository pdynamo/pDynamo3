/*==================================================================================================================================
! . Functions to make CI state characters.
!=================================================================================================================================*/

# include <math.h>
# include <stdlib.h>

# include "CIConfigurationContainer.h"
# include "DenseDeterminants.h"
# include "NumericalMacros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real CharacterDeterminant ( const Integer          activeElectrons       ,
                                   const Boolean          includeCoreOrbitals   ,
                                   const Integer          numberCoreOrbitals    ,
                                   const IntegerArray1D  *iActiveIndices        ,
                                   const IntegerArray1D  *jActiveIndices        ,
                                   const RealArray2D     *orbitalTransformation ,
                                         RealArray2D     *work                  ,
                                         Status          *status                ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . CI state characters.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . The configurations are stored as ( core with alpha and beta alternating ) ( alpha active ) ( beta active ).
! . Phases are unnecessary when cores are not included because alphas and betas are already ordered.
! . Likewise when cores are included the phase to get all core betas to the active betas is determined by whether
! . Nc * Na + ( Nc * ( Nc - 1 ) ) / 2 is odd or even. Here Nc is the number of core orbitals and Na is the number
! . of active alpha. However, these phases are always the same for non-zero interactions between states and so
! . can be ignored as they multiply to 1.
*/
void CIConfigurationContainer_Characters ( const CIConfigurationContainer *self                  ,
                                           const Boolean                   includeCoreOrbitals   ,
                                           const Integer                   coreOrbitals          ,
                                           const RealArray2D              *orbitalTransformation ,
                                                 RealArray2D              *stateTransformation   ,
                                                 Status                   *status                )
{
    if ( ( self                  != NULL ) &&
         ( orbitalTransformation != NULL ) &&
         ( stateTransformation   != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer      o ;
        auto RealArray2D *work ;
        /* . Allocate space. */
        o    = View2D_Columns ( orbitalTransformation ) ;
        work = RealArray2D_AllocateWithExtents ( o, o, status ) ;
        if ( work != NULL )
        {
            auto Integer          i, j, nAi, nAj, nBi ;
            auto IntegerArray1D  *iAlphas, *iBetas, *jAlphas, *jBetas ;
             /* . Get the state transformation. */
            RealArray2D_Set ( stateTransformation, 0.0e+00 ) ;
            /* . Double loop over configurations. */
            for ( i = 0 ; i < self->nConfigurations ; i++ )
            {
                nAi     = self->configurations[i].nAlphas ;
                nBi     = self->nElectrons - nAi ;
                iAlphas = self->configurations[i].alphas  ;
                iBetas  = self->configurations[i].betas   ;
                for ( j = 0 ; j < self->nConfigurations ; j++ )
                {
                    nAj     = self->configurations[j].nAlphas ;
                    jAlphas = self->configurations[j].alphas  ;
                    jBetas  = self->configurations[j].betas   ;
                    /* . Skip if there are different numbers of alpha orbitals in the two configurations. */
                    if ( nAi != nAj ) continue ;
                    /* . Set the value. */
                    Array2D_Item ( stateTransformation, i, j ) = CharacterDeterminant ( nAi, includeCoreOrbitals, coreOrbitals, iAlphas, jAlphas, orbitalTransformation, work, status ) *
                                                                     CharacterDeterminant ( nBi, includeCoreOrbitals, coreOrbitals, iBetas , jBetas , orbitalTransformation, work, status ) ;
                }
            }
        }
        /* . Out of space. */
        else Status_Set ( status, Status_OutOfMemory ) ;
        /* . Finish up. */
        RealArray2D_Deallocate ( &work ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Determine the alpha or beta contribution to an element of the CI state transformation matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real CharacterDeterminant ( const Integer          activeElectrons       ,
                                   const Boolean          includeCoreOrbitals   ,
                                   const Integer          numberCoreOrbitals    ,
                                   const IntegerArray1D  *iActiveIndices        ,
                                   const IntegerArray1D  *jActiveIndices        ,
                                   const RealArray2D     *orbitalTransformation ,
                                         RealArray2D     *work                  ,
                                         Status          *status                )
{
    Real determinant = 1.0e+00 ;
    if ( ( activeElectrons > 0 ) && ( work != NULL ) )
    {
        auto Integer     coreIncrement, i, iActive, iIndex, j, jActive, jIndex, numberActive, totalElectrons ;
        auto RealArray2D matrix ;
        /* . Get space. */
        totalElectrons = activeElectrons ;
        if ( includeCoreOrbitals ) totalElectrons += numberCoreOrbitals ;
        RealArray2D_View ( work, 0, 0, totalElectrons, totalElectrons, 1, 1, False, &matrix, status ) ;
        RealArray2D_Set  ( &matrix, 0.0e+00 ) ;
        /* . Cores. */
        if ( includeCoreOrbitals )
        {
            coreIncrement = numberCoreOrbitals ;
            for ( i = 0 ; i < numberCoreOrbitals ; i++ )
            {
                for ( j = 0 ; j < numberCoreOrbitals ; j++ ) Array2D_Item ( &matrix, i, j ) = Array2D_Item ( orbitalTransformation, i, j ) ;
            }
        }
        else coreIncrement = 0 ;
        /* . Active space. */
        numberActive = View1D_Extent ( iActiveIndices ) ;
        for ( i = coreIncrement, iActive = 0 ; iActive < numberActive ; iActive++ )
        {
            iIndex = Array1D_Item ( iActiveIndices, iActive ) ;
            if ( iIndex > 0 )
            {
                for ( j = coreIncrement, jActive = 0 ; jActive < numberActive ; jActive++ )
                {
                    jIndex = Array1D_Item ( jActiveIndices, jActive ) ;
                    if ( jIndex > 0 )
                    {
                        Array2D_Item ( &matrix, i, j ) = Array2D_Item ( orbitalTransformation, iActive + coreIncrement, jActive + coreIncrement ) ;
                        j++ ;
                    }
                }
                i++ ;
            }
        }
        /* . Get the determinant. */
        determinant = SquareMatrix_Determinant ( &matrix, status ) ;
    }
    return determinant ;
}
