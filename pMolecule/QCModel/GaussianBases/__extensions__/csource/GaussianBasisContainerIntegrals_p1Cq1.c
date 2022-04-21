/*==================================================================================================================================
! . Container integrals - 0 basis, 0 electron, 2 nuclei/points.
!=================================================================================================================================*/

# include <math.h>

# include "GaussianBasis.h" /* . For PI252. */
# include "GaussianNucleus.h"
# include "GaussianBasisContainerIntegrals_p1Cq1.h"
# include "RysQuadrature.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _Selected( selection, i ) ( (selection) == NULL ? True : Block_Item ( selection->flags, i ) )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Nuclear-nuclear energy and, optionally, gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real GaussianBasisContainer_m1Cn1ER1 ( const RealArray1D  *chargesI      ,
                                       const RealArray1D  *chargesJ      ,
                                       const Coordinates3 *coordinates3I ,
                                       const Coordinates3 *coordinates3J ,
                                             Selection    *selectionI    ,
                                             Selection    *selectionJ    ,
                                       const RealArray1D  *widthsEI      ,
                                       const RealArray1D  *widthsEJ      ,
                                       const RealArray1D  *widthsNI      ,
                                       const RealArray1D  *widthsNJ      ,
                                             Coordinates3 *gradients3I   ,
                                             Coordinates3 *gradients3J   )
{
    Real energy = 0.0e+00 ;
    if ( ( chargesI      != NULL ) &&
         ( chargesJ      != NULL ) &&
         ( coordinates3I != NULL ) &&
         ( coordinates3J != NULL ) )
    {
        auto Boolean       doGradients = ( ( gradients3I != NULL ) && ( gradients3J != NULL ) ), iIsJ ; 
        auto Integer       i, j, jUpper ;
        auto Real          dF, eI, eJ, factor, iAndj, ij, nI, nJ, qI, qJ, rho, r2, u2, xI, xIJ, yI, yIJ, zI, zIJ ;
        auto RysQuadrature roots ;
        Selection_MakeFlags ( selectionI, Coordinates3_Rows ( coordinates3I ), NULL ) ;
        Selection_MakeFlags ( selectionJ, Coordinates3_Rows ( coordinates3J ), NULL ) ;
        iIsJ   = ( chargesI == chargesJ ) && ( coordinates3I == coordinates3J ) ;
        jUpper = Coordinates3_Rows ( coordinates3J ) ;
        for ( i = 0 ; i < Coordinates3_Rows ( coordinates3I ) ; i++ )
        {
            if ( _Selected ( selectionI, i ) )
            {
                qI = Array1D_Item      ( chargesI, i ) ;
                eI = _GetWidthE        ( widthsEI, i ) ;
                nI = _GetWidthN        ( widthsNI, i ) ;
                xI = Coordinates3_Item ( coordinates3I, i, 0 ) ;
                yI = Coordinates3_Item ( coordinates3I, i, 1 ) ;
                zI = Coordinates3_Item ( coordinates3I, i, 2 ) ;
                if ( iIsJ ) jUpper = i ;
                for ( j = 0 ; j < jUpper ; j++ )
                {
                    if ( _Selected ( selectionJ, j ) )
                    {
                        qJ      = Array1D_Item ( chargesJ, j ) ;
                        eJ      = _GetWidthE   ( widthsEJ, j ) ;
                        nJ      = _GetWidthN   ( widthsNJ, j ) ;
                        xIJ     = xI - Coordinates3_Item ( coordinates3J, j, 0 ) ;
                        yIJ     = yI - Coordinates3_Item ( coordinates3J, j, 1 ) ;
                        zIJ     = zI - Coordinates3_Item ( coordinates3J, j, 2 ) ;
                        r2      = xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ;
                        iAndj   = eI + eJ ;
                        ij      = eI * eJ ;
                        rho     = ij / iAndj ;
                        RysQuadrature_Roots ( &roots, 1, ( rho * r2 ) ) ;
                        factor  = PI252 * nI * nJ * qI * qJ * roots.weights[0] / sqrt ( iAndj ) ;
                        energy += ( factor / ij ) ;
                        if ( doGradients )
                        {
                            u2    = rho * roots.roots[0] ;
                            dF    = - 2.0e+00 * factor * u2 / ( ij + u2 * iAndj ) ;
                            xIJ  *= dF ; yIJ *= dF ; zIJ *= dF ;
                            Coordinates3_IncrementRow ( gradients3I, i, xIJ, yIJ, zIJ ) ;
                            Coordinates3_DecrementRow ( gradients3J, j, xIJ, yIJ, zIJ ) ;
                        }
                    }
                }
            }
        }
    }
    return energy ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Nuclear-nuclear potentials (at I due to J).
! . Potentials should be initialized before entry.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainer_m1Cp1V ( const RealArray1D  *chargesJ      ,
                                     const Coordinates3 *coordinates3I ,
                                     const Coordinates3 *coordinates3J ,
                                           Selection    *selectionI    ,
                                           Selection    *selectionJ    ,
                                     const RealArray1D  *widthsEI      ,
                                     const RealArray1D  *widthsEJ      ,
                                     const RealArray1D  *widthsNI      ,
                                     const RealArray1D  *widthsNJ      ,
                                           RealArray1D  *potentialsI   )
{
    if ( ( chargesJ      != NULL ) &&
         ( coordinates3I != NULL ) &&
         ( coordinates3J != NULL ) &&
         ( potentialsI   != NULL ) )
    {
        auto Boolean       iIsJ ; 
        auto Integer       i, j, jUpper ;
        auto Real          eI, eJ, iAndj, ij, nI, nJ, pI, qJ, rho, r2, xI, xIJ, yI, yIJ, zI, zIJ ;
        auto RysQuadrature roots ;
        Selection_MakeFlags ( selectionI, Coordinates3_Rows ( coordinates3I ), NULL ) ;
        Selection_MakeFlags ( selectionJ, Coordinates3_Rows ( coordinates3J ), NULL ) ;
        iIsJ   = ( coordinates3I == coordinates3J ) ;
        jUpper = Coordinates3_Rows ( coordinates3J ) ;
        for ( i = 0 ; i < Coordinates3_Rows ( coordinates3I ) ; i++ )
        {
            if ( _Selected ( selectionI, i ) )
            {
                eI = _GetWidthE ( widthsEI, i ) ;
                nI = _GetWidthN ( widthsNI, i ) ;
                pI = 0.0e+00 ;
                xI = Coordinates3_Item ( coordinates3I, i, 0 ) ;
                yI = Coordinates3_Item ( coordinates3I, i, 1 ) ;
                zI = Coordinates3_Item ( coordinates3I, i, 2 ) ;
                if ( iIsJ ) jUpper = i ;
                for ( j = 0 ; j < jUpper ; j++ )
                {
                    if ( _Selected ( selectionJ, j ) )
                    {
                        qJ    = Array1D_Item ( chargesJ, j ) ;
                        eJ    = _GetWidthE   ( widthsEJ, j ) ;
                        nJ    = _GetWidthN   ( widthsNJ, j ) ;
                        xIJ   = xI - Coordinates3_Item ( coordinates3J, j, 0 ) ;
                        yIJ   = yI - Coordinates3_Item ( coordinates3J, j, 1 ) ;
                        zIJ   = zI - Coordinates3_Item ( coordinates3J, j, 2 ) ;
                        r2    = xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ;
                        iAndj = eI + eJ ;
                        ij    = eI * eJ ;
                        rho   = ij / iAndj ;
                        RysQuadrature_Roots ( &roots, 1, ( rho * r2 ) ) ;
                        pI   += ( qJ * PI252 * nI * nJ * roots.weights[0] / ( ij * sqrt ( iAndj ) ) ) ;
                    }
                }
                Array1D_Item ( potentialsI, i ) += pI ;
            }
        }
    }
}

# undef _Selected
