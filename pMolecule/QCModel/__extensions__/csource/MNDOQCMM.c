/*==================================================================================================================================
! . MNDO QC/MM integrals and their derivatives.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>

# include "Integer.h"
# include "Memory.h"
# include "MNDOIntegralDefinitions.h"
# include "MNDOIntegralsMM.h"
# include "MNDOQCMM.h"
# include "Real.h"
# include "Units.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _BlockSize           1024
# define _DefaultCutOff       1.0e+6
# define _ConversionFactorE   Units_Energy_Hartrees_To_Kilojoules_Per_Mole
# define _ConversionFactorG ( Units_Energy_Hartrees_To_Kilojoules_Per_Mole * Units_Length_Angstroms_To_Bohrs )
# define _NumberOfLFOEIs      10      /* . s, sp, spd: 1,  4, 10 - see NCUNIQUE MNDOIntegralDefinitions. */
# define _NumberOfMFOEIs      45      /* . s, sp, spd: 1, 10, 45 = (n*(n+1))/2. */
# define _UnderFlow           1.0e-12

/*----------------------------------------------------------------------------------------------------------------------------------
! . Gradients in normal units.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDO_QCMMGradients ( const IntegerArray1D  *atomIndices  , /* . The atom to which a basis function belongs. */
                          const SymmetricMatrix *dTotal       ,
                                BlockStorage    *integrals    ,
                                Coordinates3    *qcGradients3 ,
                                Coordinates3    *mmGradients3 ,
                                Status          *status       )
{
    if ( ( atomIndices  != NULL ) &&
         ( dTotal       != NULL ) &&
         ( integrals    != NULL ) &&
         ( mmGradients3 != NULL ) &&
         ( qcGradients3 != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Block      *block ;
        auto Cardinal16 *indices16 ;
        auto Cardinal32 *indices32 ;
        auto Integer     c, c2, c3, m, q, u, v ;
        auto Real       *data, gX, gY, gZ, p ;
        List_Iterate_Initialize ( integrals->blocks ) ;
        while ( ( block = BlockStorage_Iterate ( integrals ) ) != NULL )
        {
            data      = block->data      ;
            indices16 = block->indices16 ;
            indices32 = block->indices32 ;
            for ( c = 0 ; c < block->count ; c++ )
            {
                c2 = 2 * c ;
                c3 = 3 * c ;
                u  = indices16[c2  ] ;
                v  = indices16[c2+1] ;
                p  = SymmetricMatrix_Item ( dTotal, u, v ) ;
                if ( u != v ) p *= 2.0e+00 ; /* . Scale off-diagonal values. */
                gX = p * data [c3  ] ;
                gY = p * data [c3+1] ;
                gZ = p * data [c3+2] ;
                m  = indices32[c   ] ;
                q  = Array1D_Item ( atomIndices, u ) ;
                Coordinates3_IncrementRow ( qcGradients3, q, gX, gY, gZ ) ; /* . Positive as g terms already have -1 factor. */
                Coordinates3_DecrementRow ( mmGradients3, m, gX, gY, gZ ) ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Integrals, derivative integrals, core energy and gradients.
! . Integrals in atomic units, all other quantities in normal units.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real MNDO_QCMMIntegrals ( const MNDOParametersContainer *parameters          ,
                          const IntegerArray1D          *basisIndices        ,
                          const CubicSplineContainer    *splines             ,
                          const Real                     cutOff              ,
                          const Real                     eScale              ,
                          const Coordinates3            *qcCoordinates3      ,
                          const Coordinates3            *mmCoordinates3      ,
                          const RealArray1D             *mmCharges           ,
                                PairList                *pairList            ,
                                SymmetricMatrix         *oneElectronMatrix   ,
                                Coordinates3            *qcGradients3        ,
                                Coordinates3            *mmGradients3        ,
                                BlockStorage           **derivativeIntegrals ,
                                Status                  *status              )
{
    Real eNuclear = 0.0e+00 ;
    if ( ( parameters        != NULL    ) &&
         ( basisIndices      != NULL    ) &&
         ( eScale            != 0.0e+00 ) &&
         ( mmCharges         != NULL    ) &&
         ( mmCoordinates3    != NULL    ) &&
         ( qcCoordinates3    != NULL    ) &&
         ( pairList          != NULL    ) &&
         ( oneElectronMatrix != NULL    ) &&
         Status_IsOK ( status ) )
    {
        auto BlockStorage     *dOEIs     = NULL ;
        auto Boolean           doGradients, useSplines ;
        auto Cardinal16       *indices16 = NULL ;
        auto Cardinal32       *indices32 = NULL ;
        auto CubicSpline      *qSpline ;
        auto Integer           c, c2, c3, i0, m, n, nI, nT, q, u, v, w ;
        auto MNDOParameters   *qData ;
        auto PairListIterator  iterator ;
        auto PairRecord       *record ;
        auto Real              cutOffSquared, eNuc, gNuc, gX, gXt = 0.0e+00, gY, gYt = 0.0e+00, gZ, gZt = 0.0e+00,
                               qM, r, r2, scale, xM, xQ, xQM, yM, yQ, yQM, zM, zQ, zQM ;
        auto Real             *dIntegrals = NULL ;
        auto RealArray1D      *gLocal = NULL , *gMolecularX = NULL , *gMolecularY = NULL , *gMolecularZ = NULL ,
                              *iLocal = NULL , *iMolecular  = NULL , *iTotal      = NULL ;
        /* . Options */
        if ( cutOff > 0.0e+00 ) cutOffSquared = pow (         cutOff, 2 ) ;
        else                    cutOffSquared = pow ( _DefaultCutOff, 2 ) ;
        doGradients = ( ( derivativeIntegrals != NULL ) && ( mmGradients3 != NULL ) && ( qcGradients3 != NULL ) ) ;
        useSplines  = ( splines != NULL ) ;
        /* . Initialization. */
        if ( derivativeIntegrals != NULL ) (*derivativeIntegrals) = NULL ;
        if ( doGradients )
        {
            /* . Allocation. */
            nI         = MNDOParametersContainer_LargestBasis ( parameters ) ;
            nT         = ( nI * ( nI + 1 ) ) / 2 ;
            n          = PairList_MaximumRecordSize  ( pairList ) * nT ;
            dOEIs      = BlockStorage_Allocate       ( NULL ) ;
            dIntegrals = Memory_AllocateArrayOfTypes ( 3 * n, Real       ) ;
            indices16  = Memory_AllocateArrayOfTypes ( 2 * n, Cardinal16 ) ;
            indices32  = Memory_AllocateArrayOfTypes (     n, Cardinal32 ) ;
            if ( ( dOEIs      == NULL ) ||
                 ( dIntegrals == NULL ) ||
                 ( indices16  == NULL ) ||
                 ( indices32  == NULL ) ) { Status_Set ( status, Status_OutOfMemory ) ; goto FinishUp ; }
            /* . Initialization. */
            dOEIs->blockSize      = _BlockSize ;
            dOEIs->checkUnderFlow = True ;
            dOEIs->nIndices16     = 2 ;
            dOEIs->nIndices32     = 1 ;
            dOEIs->nReal          = 3 ;
            dOEIs->underFlow      = _UnderFlow ;
        }
        /* . Loop over QC/MM records. */
        PairListIterator_Initialize ( &iterator, pairList ) ;
        while ( ( record = PairListIterator_Next ( &iterator ) ) != NULL )
        {
            q     = record->index ;
            qData = parameters->entries[q] ;
            if ( useSplines ) qSpline = splines->entries[q] ;
            nI    = qData->norbitals ;
            if ( nI <= 0 ) continue ;
            i0    = Array1D_Item ( basisIndices, q ) ;
            nT    = ( nI * ( nI + 1 ) ) / 2 ;
            Coordinates3_GetRow ( qcCoordinates3, q, xQ, yQ, zQ ) ; /* . In Angstroms. */
            /* . Allocation. */
            iLocal      = RealArray1D_AllocateWithExtent ( nT, status ) ;
            iMolecular  = RealArray1D_AllocateWithExtent ( nT, status ) ;
            iTotal      = RealArray1D_AllocateWithExtent ( nT, status ) ;
            if ( doGradients )
            {
                gLocal      = RealArray1D_AllocateWithExtent ( nT, status ) ;
                gMolecularX = RealArray1D_AllocateWithExtent ( nT, status ) ;
                gMolecularY = RealArray1D_AllocateWithExtent ( nT, status ) ;
                gMolecularZ = RealArray1D_AllocateWithExtent ( nT, status ) ;
            }
            if ( ! Status_IsOK ( status ) ) goto Deallocate ;
            /* . Initialization. */
            RealArray1D_Set ( iTotal, 0.0e+00 ) ;
            if ( doGradients ) { gXt = gYt = gZt = 0.0e+00 ; }
            /* . Loop over interactions. */
       	    for ( c = n = 0 ; n < record->capacity ; n++ )
            {
	        m  = record->indices[n] ;
                Coordinates3_GetRow ( mmCoordinates3, m, xM, yM, zM ) ; /* . In Angstroms. */
	        xQM = xQ - xM ;
                yQM = yQ - yM ;
                zQM = zQ - zM ;
                r2  = ( xQM * xQM + yQM * yQM + zQM * zQM ) ;
                if ( r2 < cutOffSquared )
                {
                    r    = sqrt ( r2 ) * Units_Length_Angstroms_To_Bohrs ;
                    qM   = eScale * Array1D_Item ( mmCharges, m ) ;
                    xQM *=               Units_Length_Angstroms_To_Bohrs ;
                    yQM *=               Units_Length_Angstroms_To_Bohrs ;
                    zQM *=               Units_Length_Angstroms_To_Bohrs ;
                    RealArray1D_Set ( iLocal, 0.0e+00 ) ;
                    RealArray1D_Set ( gLocal, 0.0e+00 ) ;
                    /* . Core term and local integrals in atomic units. */
                    if ( useSplines )
                    {
                        MNDOIntegralsMM_FromSpline ( qData, qSpline, qM, r, &eNuc, &gNuc, iLocal, gLocal ) ;
                    }
                    else
                    {
                        auto Real eNuc0, eNuc1, gNuc0, gNuc1 ;
                        MNDOIntegralsMM_CoreCharge ( qData, qM, r, &eNuc0, &eNuc1, &gNuc0, &gNuc1 ) ; eNuc = eNuc0 + eNuc1 ; gNuc = gNuc0 + gNuc1 ;
                        MNDOIntegralsMM_LocalFrame ( qData,     r, iLocal        , gLocal         ) ;
                    }
                    eNuclear += eNuc ;
                    MNDOIntegralsMM_MolecularFrame ( nI, r, xQM, yQM, zQM, iLocal, gLocal, iMolecular, gMolecularX, gMolecularY, gMolecularZ ) ; /* . I and -dI/dx (for some reason). */
                    RealArray1D_Add ( iTotal, qM, iMolecular, NULL ) ;
/*
printf ( "Core and Integrals (local, molecule, total) for %d/%d\n", q, m ) ; fflush ( stdout ) ;
printf ( "Core, distance, charge = %20.5f %10.3f %10.3f\n", eNuc, r, qM ) ; fflush ( stdout ) ;
RealArray1D_Print ( iLocal     ) ; fflush ( stdout ) ;
RealArray1D_Print ( iMolecular ) ; fflush ( stdout ) ;
RealArray1D_Print ( iTotal     ) ; fflush ( stdout ) ;
*/
                    if ( doGradients )
                    {
                        /* . Core term. */
                        scale = _ConversionFactorG ;
                        gNuc *= ( scale / r ) ;
                        gX    = gNuc * xQM ; gXt += gX ;
                        gY    = gNuc * yQM ; gYt += gY ;
                        gZ    = gNuc * zQM ; gZt += gZ ;
                        Coordinates3_DecrementRow ( mmGradients3, m, gX, gY, gZ ) ;
                        /* . Electron term. */
                        scale *= qM ;
                        for ( u = i0, w = 0 ; u < ( i0+nI ) ; u++ )
                        {
                            for ( v = i0 ; v <= u ; c++, v++, w++ )
                            {
                                c2 = 2 * c ;
                                c3 = 3 * c ;
                                indices16 [c2  ] = u ;
                                indices16 [c2+1] = v ;
                                indices32 [c   ] = m ;
                                dIntegrals[c3  ] = scale * Array1D_Item ( gMolecularX, w ) ;
                                dIntegrals[c3+1] = scale * Array1D_Item ( gMolecularY, w ) ;
                                dIntegrals[c3+2] = scale * Array1D_Item ( gMolecularZ, w ) ;
                            }
                        }
                    }
                }
            }
            /* . Accumulate the terms for q. */
/*
printf ( "Final Integrals for %d\n", q ) ; fflush ( stdout ) ;
RealArray1D_Print ( iTotal ) ; fflush ( stdout ) ;
printf ( "Counters : %d %d %d %d\n", i0, i0 + nI, nI, nT ) ; fflush ( stdout ) ;
*/
            for ( u = i0, w = 0 ; u < ( i0+nI ) ; u++ )
            {
                for ( v = i0 ; v <= u ; v++, w++ ) SymmetricMatrix_Item ( oneElectronMatrix, u, v ) -= Array1D_Item ( iTotal, w ) ; /* . -ve as electrons. */
            }
            if ( doGradients )
            {
                Coordinates3_IncrementRow ( qcGradients3, q, gXt, gYt, gZt ) ;
                BlockStorage_AddData ( dOEIs, c, dIntegrals, indices16, indices32, status ) ;
                if ( ! Status_IsOK ( status ) ) goto FinishUp ;
            }
            /* . Finish up. */
        Deallocate:
            RealArray1D_Deallocate ( &gLocal      ) ;
            RealArray1D_Deallocate ( &gMolecularX ) ;
            RealArray1D_Deallocate ( &gMolecularY ) ;
            RealArray1D_Deallocate ( &gMolecularZ ) ;
            RealArray1D_Deallocate ( &iLocal      ) ;
            RealArray1D_Deallocate ( &iMolecular  ) ;
            RealArray1D_Deallocate ( &iTotal      ) ;
            if ( ! Status_IsOK ( status ) ) break ;
        }
        /* . Finish up. */
    FinishUp:
        eNuclear *= _ConversionFactorE ;
        if ( doGradients )
        {
            if ( Status_IsOK ( status ) ) (*derivativeIntegrals) = dOEIs ;
            else BlockStorage_Deallocate ( &dOEIs ) ;
            Memory_Deallocate ( dIntegrals ) ;
            Memory_Deallocate ( indices16  ) ;
            Memory_Deallocate ( indices32  ) ;
        }
    }
    return eNuclear ;
}

# undef _BlockSize
# undef _ConversionFactorE
# undef _ConversionFactorG
# undef _DefaultCutOff
# undef _NumberOfLFOEIs
# undef _NumberOfMFOEIs
# undef _UnderFlow
