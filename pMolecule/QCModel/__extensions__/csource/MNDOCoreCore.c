/*==================================================================================================================================
! . MNDO core-core interactions.
!=================================================================================================================================*/

# include <math.h>

# include "Boolean.h"
# include "Integer.h"
# include "MNDOCoreCore.h"
# include "MNDODefinitions.h"
# include "MNDOParameters.h"
# include "Real.h"
# include "Units.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CoreCoreInteractions ( const MNDOParameters *iData, const MNDOParameters *jData, const Real R, Real *fCore, Real *dCore ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . The core-core energy and gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Real MNDO_CoreCoreEnergy ( const MNDOParametersContainer *parameters   ,
                                  const Coordinates3            *coordinates3 ,
                                        Coordinates3            *gradients3   )
{
    Real energy = 0.0e+00 ;
    if ( ( parameters   != NULL ) &&
         ( coordinates3 != NULL ) )
    {
        auto Boolean         doGradients = ( gradients3 != NULL ) ; 
        auto Integer         i, j ;
        auto MNDOParameters *iData, *jData ;
        auto Real            f, g, *gP = NULL, R, xI, xIJ, yI, yIJ, zI, zIJ ;
        if ( doGradients ) gP = &g ;
        for ( i = 0 ; i < Coordinates3_Rows ( coordinates3 ) ; i++ )
        {
            iData = parameters->entries[i] ;
            xI    = Coordinates3_Item ( coordinates3, i, 0 ) ;
            yI    = Coordinates3_Item ( coordinates3, i, 1 ) ;
            zI    = Coordinates3_Item ( coordinates3, i, 2 ) ;
            for ( j = 0 ; j < i ; j++ )
            {
                jData = parameters->entries[j] ;
                xIJ   = Coordinates3_Item ( coordinates3, j, 0 ) - xI ;
                yIJ   = Coordinates3_Item ( coordinates3, j, 1 ) - yI ;
                zIJ   = Coordinates3_Item ( coordinates3, j, 2 ) - zI ;
                R     = sqrt ( xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ) ;
                CoreCoreInteractions ( iData, jData, R, &f, gP ) ;
                energy += f ;
                if ( doGradients )
                {
                    g   /= ( - R ) ;
                    xIJ *= g ; yIJ *= g ; zIJ *= g ;
                    Coordinates3_IncrementRow ( gradients3, i, xIJ, yIJ, zIJ ) ;
                    Coordinates3_DecrementRow ( gradients3, j, xIJ, yIJ, zIJ ) ;
                }
            }
        }
    }
    return energy ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . The core-core interactions.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CoreCoreInteractions ( const MNDOParameters *iData, const MNDOParameters *jData, const Real R, Real *fCore, Real *dCore )
{
    Integer  i, j ;
    Real     anam1, aij, ax, dd, dedr, dgdr, dscale, enuc, exi, exj, f3, gam, scale, xij, zaf, zbf ;

    /* . Initialization. */
    enuc = 0.0e+00 ;
    dedr = 0.0e+00 ;

    /* . Calculate the core integral and its derivative. */
    gam  = 1.0e+00 / sqrt ( R*R + pow ( ( iData->po[8] + jData->po[8] ), 2 ) ) ;
    dgdr = - R * gam * gam * gam ;

    /* . Diatomic terms. */
    if ( iData->QDIATOMIC && jData->QDIATOMIC )
    {
        /* . Initialization. */
        dscale = 0.0e+00 ;
        scale  = 1.0e+00 ;

        /* . Determine aij and xij parameter-dependent terms. */
        if ( ( jData->atomicNumber < iData->ndiatomic ) && iData->QDIATOMICFLAGS[jData->atomicNumber] )
        {
            aij =           iData->diatomica[jData->atomicNumber] ;
            xij = 2.0e+00 * iData->diatomicx[jData->atomicNumber] ; /* . Factor of 2 - see Stewart's am1/d paper! */

            /* . C-H, N-H and O-H. */
            if ( ( ( iData->atomicNumber == 1 ) && ( ( jData->atomicNumber == 6 ) || ( jData->atomicNumber == 7 ) || ( jData->atomicNumber == 8 ) ) ) ||
                 ( ( jData->atomicNumber == 1 ) && ( ( iData->atomicNumber == 6 ) || ( iData->atomicNumber == 7 ) || ( iData->atomicNumber == 8 ) ) ) )
            {
                f3 = xij * exp ( - aij * R * R * Units_Length_Bohrs_To_Angstroms ) ;
                dscale -= 2.0e+00 * aij * R * Units_Length_Bohrs_To_Angstroms * f3 ;
                scale  += f3 ;
            }
            /* . All others. */
            else
            {
                dd = 0.0003e+00 * pow ( R * Units_Length_Bohrs_To_Angstroms, 5 ) ;
                f3 = xij * exp ( - aij * R * ( 1.0e+00 + dd ) ) ;
                dscale -= aij * ( 1.0e+00 + 6.0e+00 * dd ) * f3 ;
                scale  += f3 ;
            }
        }

        /* . Element-specific extra terms independent of aij and xij. */
        /* . C-C. */
        if ( ( iData->atomicNumber == 6 ) && ( jData->atomicNumber == 6 ) )
        {
            f3 = 9.28e+00 * exp ( - 5.98e+00 * R * Units_Length_Bohrs_To_Angstroms ) ;
            dscale -= 5.98e+00 * Units_Length_Bohrs_To_Angstroms * f3 ;
            scale  += f3 ;
        }
        /* . Si-O. */
        if ( ( ( iData->atomicNumber ==  8 ) && ( jData->atomicNumber == 14 ) ) ||
             ( ( iData->atomicNumber == 14 ) && ( jData->atomicNumber ==  8 ) ) )
        {
            dd = R * Units_Length_Bohrs_To_Angstroms - 2.9e+00 ;
            f3 = - 0.0007e+00 * exp ( - pow ( dd, 2 ) ) ;
            dscale -= 2.0e+00 * dd * Units_Length_Bohrs_To_Angstroms * f3 ;
            scale  += f3 ;
        }

        /* . Initial term. */
        enuc   = iData->zcore * jData->zcore * gam * scale ;
        dedr   = iData->zcore * jData->zcore * ( dgdr * scale + gam * dscale ) ;

        /* . Unpolarizable core. */
        enuc  += PM6_UNPOLARIZABLECORE * pow ( ( ( pow ( iData->zcore, 1.0e+00 / 3.0e+00 ) + pow ( jData->zcore, 1.0e+00 / 3.0e+00 ) ) / ( R * Units_Length_Bohrs_To_Angstroms ) ), 12 ) ;
        dedr  -= 12.0e+00 * PM6_UNPOLARIZABLECORE * pow ( ( ( pow ( iData->zcore, 1.0e+00 / 3.0e+00 ) + pow ( jData->zcore, 1.0e+00 / 3.0e+00 ) ) / ( R * Units_Length_Bohrs_To_Angstroms ) ), 12 ) / R ;
        scale  = 0.0e+00 ;
    }
    /* . Monoatomic terms. */
    else
    {
        exi   = exp ( -iData->alp * R ) ;
        exj   = exp ( -jData->alp * R ) ;
        scale = exi + exj ;
	if ( ( iData->atomicNumber == 1 ) && ( ( jData->atomicNumber == 7 ) || ( jData->atomicNumber == 8 ) ) )
        {
	    f3     = 1.0e+00 + exi + Units_Length_Bohrs_To_Angstroms * R * exj ;
	    dd     = dgdr * f3 - gam * ( iData->alp * exi + Units_Length_Bohrs_To_Angstroms * ( jData->alp * R - 1.0e+00 ) * exj ) ;
            scale += ( Units_Length_Bohrs_To_Angstroms * R - 1.0e+00 ) * exj ;
	}
        else if ( ( ( iData->atomicNumber == 7 ) || ( iData->atomicNumber == 8 ) ) && ( jData->atomicNumber == 1 ) )
        {
	    f3     = 1.0e+00 + exj + Units_Length_Bohrs_To_Angstroms * R * exi ;
	    dd     = dgdr * f3 - gam * ( jData->alp * exj + Units_Length_Bohrs_To_Angstroms * ( iData->alp * R - 1.0e+00 ) * exi ) ;
            scale += ( Units_Length_Bohrs_To_Angstroms * R - 1.0e+00 ) * exi ;
	}
        else
        {
	    f3     = 1.0e+00 + exi + exj ;
            dd     = dgdr * f3 - gam * ( iData->alp * exi + jData->alp * exj ) ;
	}
	dedr  = iData->zcore * jData->zcore * dd ;
        enuc  = iData->zcore * jData->zcore * gam ;
        scale = fabs ( scale * enuc ) ;
    }

    /* . Compute the AM1/PM3-specific terms. */
    anam1 = 0.0e+00 ;
    for ( i = 0 ; i < iData->nam1pm3g ; i++ )
    {
        dd = R - iData->fn3[i] ;
        ax = iData->fn2[i] * dd * dd ;
        if ( ax <= EXPONENT_TOLERANCE )
        {
            anam1 += iData->fn1[i] * ( 1.0e+00 / ( R * R ) + iData->fn2[i] * 2.0e+00 * dd / R ) * exp( -ax ) ;
            scale += iData->gphot * jData->gphot * iData->zcore * jData->zcore / R * iData->fn1[i] * exp( -ax ) ;
        }
    }
    for ( i = 0 ; i < jData->nam1pm3g ; i++ )
    {
        dd = R - jData->fn3[i] ;
        ax = jData->fn2[i] * dd * dd ;
        if ( ax <= EXPONENT_TOLERANCE )
        {
            anam1 += jData->fn1[i] * ( 1.0e+00 / ( R * R ) + jData->fn2[i] * 2.0e+00 * dd / R ) * exp( -ax ) ;
            scale += iData->gphot * jData->gphot * iData->zcore * jData->zcore / R * jData->fn1[i] * exp( -ax ) ;
        }
    }
    dedr -= anam1 * iData->gphot * jData->gphot * iData->zcore * jData->zcore ;

    /* . Compute the PDDG-specific terms. */
    if ( ( iData->npddg > 0 ) && ( jData->npddg > 0 ) )
    {
        anam1 = 0.0e+00 ;
        zaf   = iData->zcore / ( iData->zcore + jData->zcore ) ;
        zbf   = jData->zcore / ( iData->zcore + jData->zcore ) ;
        for ( i = 0 ; i < iData->npddg ; i++ )
        {
            for ( j = 0 ; j < jData->npddg ; j++ )
            {
                dd = R - iData->pddge[i] - jData->pddge[j] ;
                ax = PDDG_EXPONENT * dd * dd ;
                anam1 += ( zaf * iData->pddgc[i] + zbf * jData->pddgc[j] ) * 2.0e+00 * PDDG_EXPONENT * dd * exp ( -ax ) ;
                scale += ( zaf * iData->pddgc[i] + zbf * jData->pddgc[j] ) * exp ( -ax ) ;
            }
        }
        dedr -= anam1 ;
    }

    /* . Add in the correction factor. */
    enuc += scale ;

    /* . Finish up. */
    if ( dCore != NULL ) (*dCore) = dedr ;
    if ( fCore != NULL ) (*fCore) = enuc ;
}

