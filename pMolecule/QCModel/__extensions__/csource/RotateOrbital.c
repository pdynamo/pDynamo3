/*==================================================================================================================================
! . A function for rotating orbitals.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "Integer.h"
# include "NumericalMacros.h"
# include "Real.h"
# include "RealArray2D.h"
# include "RotateOrbital.h"


/*----------------------------------------------------------------------------------------------------------------------------------
! . Rotate an orbital by either a proper or an improper rotation.
! . Only valid for MNDO currently (or minimal basis sets with d functions or less).
!---------------------------------------------------------------------------------------------------------------------------------*/
void RotateOrbital (       IntegerArray1D  *orbitalBasisIndices ,
                     const Matrix33        *rotation            ,
                     const IntegerArray1D  *mapping             ,
                     const RealArray1D     *inOrbital           ,
                           RealArray1D     *outOrbital          )
{
    if ( ( orbitalBasisIndices != NULL ) &&
         ( rotation            != NULL ) &&
         ( mapping             != NULL ) &&
         ( inOrbital           != NULL ) &&
         ( outOrbital          != NULL ) )
    {
        auto Integer      iAtom, iFirstOrbital, jAtom, jFirstOrbital, nAtoms, numberOrbitals ;
        auto RealArray1D  inView, outView ;
        auto RealArray2D *dTransformation = NULL, *pTransformation = NULL ;
        /* . Initialization. */
        RealArray1D_Set ( outOrbital, 0.0e+00 ) ;
        /* . Find maximum number of orbitals per atom and calculate the appropriate transformation matrices. */
        nAtoms         = View1D_Extent ( orbitalBasisIndices ) - 1 ;
        numberOrbitals = 0 ;
        for ( iAtom = 0 ; iAtom < nAtoms ; iAtom++ )
        {
            numberOrbitals = Maximum ( numberOrbitals, ( Array1D_Item ( orbitalBasisIndices, iAtom+1 ) - Array1D_Item ( orbitalBasisIndices, iAtom ) ) ) ;
        }
        /* . p and d transformations. */
        if ( numberOrbitals > 1 )
        {
            auto Real r00, r0m, r0p, rm0, rmm, rmp, rp0, rpm, rpp ;
            /* . Get matrix elements. */
            r00 = Matrix33_Item ( rotation, 2, 2 ) ;
            r0p = Matrix33_Item ( rotation, 2, 0 ) ;
            r0m = Matrix33_Item ( rotation, 2, 1 ) ;
            rp0 = Matrix33_Item ( rotation, 0, 2 ) ;
            rpp = Matrix33_Item ( rotation, 0, 0 ) ;
            rpm = Matrix33_Item ( rotation, 0, 1 ) ;
            rm0 = Matrix33_Item ( rotation, 1, 2 ) ;
            rmp = Matrix33_Item ( rotation, 1, 0 ) ;
            rmm = Matrix33_Item ( rotation, 1, 1 ) ;
            /* . p transformation - 10, 11, 1-1 = z, x, y. */
            pTransformation = RealArray2D_AllocateWithExtents ( 3, 3, NULL ) ;
            Array2D_Item ( pTransformation, 0, 0 ) = r00 ;
            Array2D_Item ( pTransformation, 0, 1 ) = r0p ;
            Array2D_Item ( pTransformation, 0, 2 ) = r0m ;
            Array2D_Item ( pTransformation, 1, 0 ) = rp0 ;
            Array2D_Item ( pTransformation, 1, 1 ) = rpp ;
            Array2D_Item ( pTransformation, 1, 2 ) = rpm ;
            Array2D_Item ( pTransformation, 2, 0 ) = rm0 ;
            Array2D_Item ( pTransformation, 2, 1 ) = rmp ;
            Array2D_Item ( pTransformation, 2, 2 ) = rmm ;
            /* . d transformation - 20, 21, 2-1, 22, 2-2. */
            if ( numberOrbitals > 4 )
            {
                auto Real sqrt3 ;
                sqrt3 = sqrt ( 3.0e+00 ) ;
                dTransformation = RealArray2D_AllocateWithExtents ( 5, 5, NULL ) ;
                Array2D_Item ( dTransformation, 0, 0 ) = ( 3.0e+00 * r00 * r00 - 1.0e+00 ) / 2.0e+00 ;
                Array2D_Item ( dTransformation, 0, 1 ) =  sqrt3 * r00 * r0p ;
                Array2D_Item ( dTransformation, 0, 2 ) =  sqrt3 * r00 * r0m ;
                Array2D_Item ( dTransformation, 0, 3 ) =  sqrt3 * ( r0p * r0p - r0m * r0m ) / 2.0e+00 ;
                Array2D_Item ( dTransformation, 0, 4 ) =  sqrt3 * r0p * r0m ;
                Array2D_Item ( dTransformation, 1, 0 ) =  sqrt3 * rp0 * r00 ;
                Array2D_Item ( dTransformation, 1, 1 ) =  rpp * r00 + rp0 * r0p ;
                Array2D_Item ( dTransformation, 1, 2 ) =  rpm * r00 + rp0 * r0m ;
                Array2D_Item ( dTransformation, 1, 3 ) =  rpp * r0p - rpm * r0m ;
                Array2D_Item ( dTransformation, 1, 4 ) =  rpp * r0m + r0p * rpm ;
                Array2D_Item ( dTransformation, 2, 0 ) =  sqrt3 * rm0 * r00 ;
                Array2D_Item ( dTransformation, 2, 1 ) =  rmp * r00 + r0p * rm0 ;
                Array2D_Item ( dTransformation, 2, 2 ) =  rmm * r00 + r0m * rm0 ;
                Array2D_Item ( dTransformation, 2, 3 ) =  rmp * r0p - rmm * r0m ;
                Array2D_Item ( dTransformation, 2, 4 ) =  rmp * r0m + r0p * rmm ;
                Array2D_Item ( dTransformation, 3, 0 ) =  sqrt3 * ( rp0 * rp0 - rm0 * rm0 ) / 2.0e+00 ;
                Array2D_Item ( dTransformation, 3, 1 ) =  rpp * rp0 - rmp * rm0 ;
                Array2D_Item ( dTransformation, 3, 2 ) =  rpm * rp0 - rmm * rm0 ;
                Array2D_Item ( dTransformation, 3, 3 ) =  ( rpp * rpp + rmm * rmm - rmp * rmp - rpm * rpm ) / 2.0e+00 ;
                Array2D_Item ( dTransformation, 3, 4 ) =  rpp * rpm - rmp * rmm ;
                Array2D_Item ( dTransformation, 4, 0 ) =  sqrt3 * rp0 * rm0 ;
                Array2D_Item ( dTransformation, 4, 1 ) =  rpp * rm0 + rmp * rp0 ;
                Array2D_Item ( dTransformation, 4, 2 ) =  rpm * rm0 + rmm * rp0 ;
                Array2D_Item ( dTransformation, 4, 3 ) =  rpp * rmp - rpm * rmm ;
                Array2D_Item ( dTransformation, 4, 4 ) =  rpp * rmm + rmp * rpm ;
            }
        }
        /* . Loop over atoms and rotate each block of orbitals separately. */
        for ( iAtom = 0 ; iAtom < nAtoms ; iAtom++ )
        {
            iFirstOrbital  =   Array1D_Item ( orbitalBasisIndices, iAtom   ) ;
            numberOrbitals = ( Array1D_Item ( orbitalBasisIndices, iAtom+1 ) - iFirstOrbital ) ;
            if ( numberOrbitals > 0 )
            {
                jAtom         = Array1D_Item ( mapping, iAtom ) ;
                jFirstOrbital = Array1D_Item ( orbitalBasisIndices, jAtom ) ;
                /* . s. */
                Array1D_Item ( outOrbital, jFirstOrbital ) = Array1D_Item ( inOrbital, iFirstOrbital ) ;
                /* . p. */
                if ( numberOrbitals > 1 )
                {
                    RealArray1D_View (  inOrbital, iFirstOrbital+1, 3, 1, False,  &inView, NULL ) ;
                    RealArray1D_View ( outOrbital, jFirstOrbital+1, 3, 1, False, &outView, NULL ) ;
                    RealArray2D_VectorMultiply ( False, 1.0e+00, pTransformation, &inView, 0.0e+00, &outView, NULL ) ;
                }
                /* . d. */
                if ( numberOrbitals > 4 )
                {
                    RealArray1D_View (  inOrbital, iFirstOrbital+4, 5, 1, False,  &inView, NULL ) ;
                    RealArray1D_View ( outOrbital, jFirstOrbital+4, 5, 1, False, &outView, NULL ) ;
                    RealArray2D_VectorMultiply ( False, 1.0e+00, dTransformation, &inView, 0.0e+00, &outView, NULL ) ;
                }
            }
        }
        /* . Finish up. */
        RealArray2D_Deallocate ( &dTransformation ) ;
        RealArray2D_Deallocate ( &pTransformation ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the rotations up to l = 2.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RotateOrbital_MakeLRotations ( const Integer      L      ,
                                    const Matrix33    *R      ,
                                          RealArray2D *T      ,
                                          Status      *status )
{
    if ( ( L >= 0    ) &&
         ( R != NULL ) &&
         ( T != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer d ;
        d = ( L + 1 ) * ( L + 1 ) ; /* . Dimension of T. */
        if ( ( View2D_Columns ( T ) == d ) && ( View2D_Rows ( T ) == d ) )
        {
            RealArray2D_Set ( T, 0.0e+00 ) ;
            /* . s transformation - 00. */
            Array2D_Item ( T, 0, 0 ) = 1.0e+00 ;
            if ( L > 0 )
            {
                auto Real r00, r0m, r0p, rm0, rmm, rmp, rp0, rpm, rpp ;
                /* . Get matrix elements. */
                r00 = Matrix33_Item ( R, 2, 2 ) ;
                r0p = Matrix33_Item ( R, 2, 0 ) ;
                r0m = Matrix33_Item ( R, 2, 1 ) ;
                rp0 = Matrix33_Item ( R, 0, 2 ) ;
                rpp = Matrix33_Item ( R, 0, 0 ) ;
                rpm = Matrix33_Item ( R, 0, 1 ) ;
                rm0 = Matrix33_Item ( R, 1, 2 ) ;
                rmp = Matrix33_Item ( R, 1, 0 ) ;
                rmm = Matrix33_Item ( R, 1, 1 ) ;
                /* . p transformation - 10, 11, 1-1 = z, x, y. */
                Array2D_Item ( T, 1, 1 ) = r00 ;
                Array2D_Item ( T, 1, 2 ) = r0p ;
                Array2D_Item ( T, 1, 3 ) = r0m ;
                Array2D_Item ( T, 2, 1 ) = rp0 ;
                Array2D_Item ( T, 2, 2 ) = rpp ;
                Array2D_Item ( T, 2, 3 ) = rpm ;
                Array2D_Item ( T, 3, 1 ) = rm0 ;
                Array2D_Item ( T, 3, 2 ) = rmp ;
                Array2D_Item ( T, 3, 3 ) = rmm ;
                /* . d transformation - 20, 21, 2-1, 22, 2-2. */
                if ( L > 1 )
                {
                    auto Real sqrt3 ;
                    sqrt3 = sqrt ( 3.0e+00 ) ;
                    Array2D_Item ( T, 4, 4 ) = ( 3.0e+00 * r00 * r00 - 1.0e+00 ) / 2.0e+00 ;
                    Array2D_Item ( T, 4, 5 ) =  sqrt3 * r00 * r0p ;
                    Array2D_Item ( T, 4, 6 ) =  sqrt3 * r00 * r0m ;
                    Array2D_Item ( T, 4, 7 ) =  sqrt3 * ( r0p * r0p - r0m * r0m ) / 2.0e+00 ;
                    Array2D_Item ( T, 4, 8 ) =  sqrt3 * r0p * r0m ;
                    Array2D_Item ( T, 5, 4 ) =  sqrt3 * rp0 * r00 ;
                    Array2D_Item ( T, 5, 5 ) =  rpp * r00 + rp0 * r0p ;
                    Array2D_Item ( T, 5, 6 ) =  rpm * r00 + rp0 * r0m ;
                    Array2D_Item ( T, 5, 7 ) =  rpp * r0p - rpm * r0m ;
                    Array2D_Item ( T, 5, 8 ) =  rpp * r0m + r0p * rpm ;
                    Array2D_Item ( T, 6, 4 ) =  sqrt3 * rm0 * r00 ;
                    Array2D_Item ( T, 6, 5 ) =  rmp * r00 + r0p * rm0 ;
                    Array2D_Item ( T, 6, 6 ) =  rmm * r00 + r0m * rm0 ;
                    Array2D_Item ( T, 6, 7 ) =  rmp * r0p - rmm * r0m ;
                    Array2D_Item ( T, 6, 8 ) =  rmp * r0m + r0p * rmm ;
                    Array2D_Item ( T, 7, 4 ) =  sqrt3 * ( rp0 * rp0 - rm0 * rm0 ) / 2.0e+00 ;
                    Array2D_Item ( T, 7, 5 ) =  rpp * rp0 - rmp * rm0 ;
                    Array2D_Item ( T, 7, 6 ) =  rpm * rp0 - rmm * rm0 ;
                    Array2D_Item ( T, 7, 7 ) =  ( rpp * rpp + rmm * rmm - rmp * rmp - rpm * rpm ) / 2.0e+00 ;
                    Array2D_Item ( T, 7, 8 ) =  rpp * rpm - rmp * rmm ;
                    Array2D_Item ( T, 8, 4 ) =  sqrt3 * rp0 * rm0 ;
                    Array2D_Item ( T, 8, 5 ) =  rpp * rm0 + rmp * rp0 ;
                    Array2D_Item ( T, 8, 6 ) =  rpm * rm0 + rmm * rp0 ;
                    Array2D_Item ( T, 8, 7 ) =  rpp * rmp - rpm * rmm ;
                    Array2D_Item ( T, 8, 8 ) =  rpp * rmm + rmp * rpm ;
                }
            }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}
