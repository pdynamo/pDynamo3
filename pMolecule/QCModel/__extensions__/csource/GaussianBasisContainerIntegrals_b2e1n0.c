/*==================================================================================================================================
! . Container integrals - 2 basis, 1 electrons, 0 nuclei/points.
!=================================================================================================================================*/

# include <math.h>
# include <stdlib.h>

# include "GaussianBasisContainerIntegrals_b2e1n0.h"
# include "GaussianBasisIntegrals_b2e1n0.h"
# include "Integer.h"
# include "RealArray2D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Coulomb integrals.
! . Integrals is overwritten by this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainerIntegrals_2Coulomb ( const GaussianBasisContainer *self         ,
                                                const IntegerArray1D         *basisIndices ,
                                                const Coordinates3           *coordinates3 ,
                                                      SymmetricMatrix        *integrals    ,
                                                      Status                 *status       )
{
    SymmetricMatrix_Set ( integrals, 0.0e+00 ) ;
    if ( ( self         != NULL ) &&
         ( basisIndices != NULL ) &&
         ( coordinates3 != NULL ) &&
         ( integrals    != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer      n ;
        auto RealArray2D *block ;
        n     = GaussianBasisContainer_LargestBasis ( self, True ) ;
        block = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        if ( block != NULL )
        {
            auto Integer        i, i0, j, j0, nI, nJ, u, v, vUpper ;
            auto GaussianBasis *iBasis, *jBasis ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                iBasis = self->entries[i] ;
                i0     = Array1D_Item ( basisIndices, i   )      ;
                nI     = Array1D_Item ( basisIndices, i+1 ) - i0 ;
                for ( j = 0 ; j <= i ; j++ )
                {
                    jBasis = self->entries[j] ;
                    j0     = Array1D_Item ( basisIndices, j   )      ;
                    nJ     = Array1D_Item ( basisIndices, j+1 ) - j0 ;
                    GaussianBasisIntegrals_2Coulomb ( iBasis                                      ,
                                                      Coordinates3_RowPointer ( coordinates3, i ) ,
                                                      jBasis                                      ,
                                                      Coordinates3_RowPointer ( coordinates3, j ) ,
                                                      block                                       ) ;
                    vUpper = nJ ;
                    for ( u = 0 ; u < nI ; u++ )
                    {
                        if ( i == j ) vUpper = u + 1 ;
                        for ( v = 0 ; v < vUpper ; v++ ) {
                        SymmetricMatrix_Item ( integrals, u+i0, v+j0 ) = Array2D_Item ( block, u, v ) ; }
                    }
                }
            }
/*
printf ( "\nCoulomb 2D Integrals (%d):\n", self->capacity ) ;
SymmetricMatrix_Print ( integrals ) ;
*/
        }
        RealArray2D_Deallocate ( &block ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Coulomb integral derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainerIntegrals_2CoulombD ( const GaussianBasisContainer *self         ,
                                                 const IntegerArray1D         *basisIndices ,
                                                 const Coordinates3           *coordinates3 ,
                                                 const RealArray1D            *fPotential   ,
                                                 const RealArray1D            *wVector      ,
                                                       Coordinates3           *gradients3   ,
                                                       Status                 *status       )
{
    if ( ( self         != NULL ) &&
         ( basisIndices != NULL ) &&
         ( coordinates3 != NULL ) &&
         ( fPotential   != NULL ) &&
         ( gradients3   != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer      n ;
        auto RealArray2D *blockX, *blockY, *blockZ, *fW ;
        n      = GaussianBasisContainer_LargestBasis ( self, True ) ;
        blockX = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        blockY = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        blockZ = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        fW     = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        if ( ( blockX != NULL ) &&
             ( blockY != NULL ) &&
             ( blockZ != NULL ) &&
             ( fW     != NULL ) )
        {
            auto Integer        i, i0, j, j0, nI, nJ, u, v ;
            auto GaussianBasis *iBasis, *jBasis ;
            auto Real           d, dX, dY, dZ, fU, fV, wU, wV ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                iBasis = self->entries[i] ;
                i0     = Array1D_Item ( basisIndices, i   )      ;
                nI     = Array1D_Item ( basisIndices, i+1 ) - i0 ;
                for ( j = 0 ; j < i ; j++ ) /* . Diagonal terms are zero. */
                {
                    jBasis = self->entries[j] ;
                    j0     = Array1D_Item ( basisIndices, j   )      ;
                    nJ     = Array1D_Item ( basisIndices, j+1 ) - j0 ;
                    GaussianBasisIntegrals_2CoulombD ( iBasis                                      ,
                                                       Coordinates3_RowPointer ( coordinates3, i ) ,
                                                       jBasis                                      ,
                                                       Coordinates3_RowPointer ( coordinates3, j ) ,
                                                       blockX                                      ,
                                                       blockY                                      ,
                                                       blockZ                                      ) ;
                    if ( wVector == NULL )
                    {
                        for ( u = 0 ; u < nI ; u++ )
                        {
                            fU = Array1D_Item ( fPotential, u+i0 ) ;
                            for ( v = 0 ; v < nJ ; v++ )
                            {
                                fV = Array1D_Item ( fPotential, v+j0 ) ;
                                Array2D_Item ( fW, u, v ) = - fU * fV ;
                            }
                        }
                    }
                    else
                    {
                        for ( u = 0 ; u < nI ; u++ )
                        {
                            fU = Array1D_Item ( fPotential, u+i0 ) ;
                            wU = Array1D_Item ( wVector   , u+i0 ) ;
                            for ( v = 0 ; v < nJ ; v++ )
                            {
                                fV = Array1D_Item ( fPotential, v+j0 ) ;
                                wV = Array1D_Item ( wVector   , v+j0 ) ;
                                Array2D_Item ( fW, u, v ) = - ( ( fU + wU ) * ( fV + wV ) - wU * wV ) ;
                            }
                        }
                    }
                    dX = dY = dZ = 0.0e+00 ;
                    for ( u = 0 ; u < nI ; u++ )
                    {
                        for ( v = 0 ; v < nJ ; v++ )
                        {
                            d   = Array2D_Item ( fW, u, v ) ;
                            dX += d * Array2D_Item ( blockX, u, v ) ;
                            dY += d * Array2D_Item ( blockY, u, v ) ;
                            dZ += d * Array2D_Item ( blockZ, u, v ) ;
                        }
                    }
                    Coordinates3_IncrementRow ( gradients3, i, dX, dY, dZ ) ;
                    Coordinates3_DecrementRow ( gradients3, j, dX, dY, dZ ) ;
                }
            }
        }
        RealArray2D_Deallocate ( &blockX ) ;
        RealArray2D_Deallocate ( &blockY ) ;
        RealArray2D_Deallocate ( &blockZ ) ;
        RealArray2D_Deallocate ( &fW     ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Overlap integrals.
! . Overlap is overwritten by this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainerIntegrals_2Overlap ( const GaussianBasisContainer *self         ,
                                                const IntegerArray1D         *basisIndices ,
                                                const Coordinates3           *coordinates3 ,
                                                      SymmetricMatrix        *overlap      ,
                                                      Status                 *status       )
{
    SymmetricMatrix_Set ( overlap, 0.0e+00 ) ;
    if ( ( self         != NULL ) &&
         ( basisIndices != NULL ) &&
         ( coordinates3 != NULL ) &&
         ( overlap      != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer      n ;
        auto RealArray2D *block ;
        n     = GaussianBasisContainer_LargestBasis ( self, True ) ;
        block = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        if ( block != NULL )
        {
            auto Integer        i, i0, j, j0, nI, nJ, u, v, vUpper ;
            auto GaussianBasis *iBasis, *jBasis ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                iBasis = self->entries[i] ;
                i0     = Array1D_Item ( basisIndices, i   )      ;
                nI     = Array1D_Item ( basisIndices, i+1 ) - i0 ;
                for ( j = 0 ; j <= i ; j++ )
                {
                    jBasis = self->entries[j] ;
                    j0     = Array1D_Item ( basisIndices, j   )      ;
                    nJ     = Array1D_Item ( basisIndices, j+1 ) - j0 ;
                    GaussianBasisIntegrals_2Overlap ( iBasis                                      ,
                                                      Coordinates3_RowPointer ( coordinates3, i ) ,
                                                      jBasis                                      ,
                                                      Coordinates3_RowPointer ( coordinates3, j ) ,
                                                      block                                       ) ;
                    vUpper = nJ ;
                    for ( u = 0 ; u < nI ; u++ )
                    {
                        if ( i == j ) vUpper = u + 1 ;
                        for ( v = 0 ; v < vUpper ; v++ ) { SymmetricMatrix_Item ( overlap, u+i0, v+j0 ) = Array2D_Item ( block, u, v ) ; }
                    }
                }
            }
        }
        RealArray2D_Deallocate ( &block ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Dipole integrals.
! . The dipole matrices are overwritten by this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainerIntegrals_Dipole ( const GaussianBasisContainer *self         ,
                                              const IntegerArray1D         *basisIndices ,
                                              const Coordinates3           *coordinates3 ,
                                                    Vector3                *center       ,
                                                    SymmetricMatrix        *dipoleX      ,
                                                    SymmetricMatrix        *dipoleY      ,
                                                    SymmetricMatrix        *dipoleZ      ,
                                                    Status                 *status       )
{
    SymmetricMatrix_Set ( dipoleX, 0.0e+00 ) ;
    SymmetricMatrix_Set ( dipoleY, 0.0e+00 ) ;
    SymmetricMatrix_Set ( dipoleZ, 0.0e+00 ) ;
    if ( ( self         != NULL ) &&
         ( basisIndices != NULL ) &&
         ( coordinates3 != NULL ) &&
         ( dipoleX      != NULL ) &&
         ( dipoleY      != NULL ) &&
         ( dipoleZ      != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer      n ;
        auto RealArray2D *blockX, *blockY, *blockZ ;
        auto Vector3     *origin ;
        n      = GaussianBasisContainer_LargestBasis ( self, True ) ;
        blockX = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        blockY = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        blockZ = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        if ( center == NULL ) { origin = Vector3_Allocate ( ) ; Vector3_Set ( origin, 0.0e+00 ) ; }
        else                    origin = center ;
        if ( ( blockX != NULL ) &&
             ( blockY != NULL ) &&
             ( blockZ != NULL ) &&
             ( origin != NULL ) )
        {
            auto Integer        i, i0, j, j0, nI, nJ, u, v, vUpper ;
            auto GaussianBasis *iBasis, *jBasis ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                iBasis = self->entries[i] ;
                i0     = Array1D_Item ( basisIndices, i   )      ;
                nI     = Array1D_Item ( basisIndices, i+1 ) - i0 ;
                for ( j = 0 ; j <= i ; j++ )
                {
                    jBasis = self->entries[j] ;
                    j0     = Array1D_Item ( basisIndices, j   )      ;
                    nJ     = Array1D_Item ( basisIndices, j+1 ) - j0 ;
                    GaussianBasisIntegrals_Dipole ( iBasis                                      ,
                                                    Coordinates3_RowPointer ( coordinates3, i ) ,
                                                    jBasis                                      ,
                                                    Coordinates3_RowPointer ( coordinates3, j ) ,
                                                    Vector3_Data ( origin )                     ,
                                                    blockX                                      ,
                                                    blockY                                      ,
                                                    blockZ                                      ) ;
                    vUpper = nJ ;
                    for ( u = 0 ; u < nI ; u++ )
                    {
                        if ( i == j ) vUpper = u + 1 ;
                        for ( v = 0 ; v < vUpper ; v++ )
                        {
                            SymmetricMatrix_Item ( dipoleX, u+i0, v+j0 ) = Array2D_Item ( blockX, u, v ) ;
                            SymmetricMatrix_Item ( dipoleY, u+i0, v+j0 ) = Array2D_Item ( blockY, u, v ) ;
                            SymmetricMatrix_Item ( dipoleZ, u+i0, v+j0 ) = Array2D_Item ( blockZ, u, v ) ;
                        }
                    }
                }
            }
        }
        RealArray2D_Deallocate ( &blockX ) ;
        RealArray2D_Deallocate ( &blockY ) ;
        RealArray2D_Deallocate ( &blockZ ) ;
        if ( center == NULL ) Vector3_Deallocate ( &origin ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Kinetic and overlap integrals.
! . The matrices must be initialized on entry to this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainerIntegrals_Kinetic2Overlap ( const GaussianBasisContainer *self         ,
                                                       const IntegerArray1D         *basisIndices ,
                                                       const Coordinates3           *coordinates3 ,
                                                             SymmetricMatrix        *kinetic      ,
                                                             SymmetricMatrix        *overlap      ,
                                                             Status                 *status       )
{
    if ( ( self         != NULL ) &&
         ( basisIndices != NULL ) &&
         ( coordinates3 != NULL ) &&
         ( kinetic      != NULL ) &&
         ( overlap      != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer      n ;
        auto RealArray2D *blockS, *blockT ;
        n      = GaussianBasisContainer_LargestBasis ( self, True ) ;
        blockS = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        blockT = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        if ( ( blockS != NULL ) &&
             ( blockT != NULL ) )
        {
            auto Integer        i, i0, j, j0, nI, nJ, u, v, vUpper ;
            auto GaussianBasis *iBasis, *jBasis ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                iBasis = self->entries[i] ;
                i0     = Array1D_Item ( basisIndices, i   )      ;
                nI     = Array1D_Item ( basisIndices, i+1 ) - i0 ;
                for ( j = 0 ; j <= i ; j++ )
                {
                    jBasis = self->entries[j] ;
                    j0     = Array1D_Item ( basisIndices, j   )      ;
                    nJ     = Array1D_Item ( basisIndices, j+1 ) - j0 ;
                    GaussianBasisIntegrals_Kinetic2Overlap ( iBasis                                      ,
                                                             Coordinates3_RowPointer ( coordinates3, i ) ,
                                                             jBasis                                      ,
                                                             Coordinates3_RowPointer ( coordinates3, j ) ,
                                                             blockS                                      ,
                                                             blockT                                      ) ;
                    vUpper = nJ ;
                    for ( u = 0 ; u < nI ; u++ )
                    {
                        if ( i == j ) vUpper = u + 1 ;
                        for ( v = 0 ; v < vUpper ; v++ )
                        {
                            SymmetricMatrix_Item ( kinetic, u+i0, v+j0 ) += Array2D_Item ( blockT, u, v ) ;
                            SymmetricMatrix_Item ( overlap, u+i0, v+j0 ) += Array2D_Item ( blockS, u, v ) ;
                        }
                    }
                }
            }
        }
        RealArray2D_Deallocate ( &blockS ) ;
        RealArray2D_Deallocate ( &blockT ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Kinetic and overlap derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainerIntegrals_Kinetic2OverlapD ( const GaussianBasisContainer *self         ,
                                                        const IntegerArray1D         *basisIndices ,
                                                        const Coordinates3           *coordinates3 ,
                                                        const SymmetricMatrix        *kDensity     ,
                                                        const SymmetricMatrix        *oDensity     ,
                                                              Coordinates3           *gradients3   ,
                                                              Status                 *status       )
{
    if ( ( self         != NULL ) &&
         ( basisIndices != NULL ) &&
         ( coordinates3 != NULL ) &&
         ( kDensity     != NULL ) &&
         ( oDensity     != NULL ) &&
         ( gradients3   != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer      n ;
        auto RealArray2D *blockSX, *blockSY, *blockSZ ,
                         *blockTX, *blockTY, *blockTZ ;
        n      = GaussianBasisContainer_LargestBasis ( self, True ) ;
        blockSX = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        blockSY = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        blockSZ = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        blockTX = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        blockTY = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        blockTZ = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        if ( ( blockSX != NULL ) &&
             ( blockSY != NULL ) &&
             ( blockSZ != NULL ) &&
             ( blockTX != NULL ) &&
             ( blockTY != NULL ) &&
             ( blockTZ != NULL ) )
        {
            auto Integer        i, i0, j, j0, nI, nJ, u, v ;
            auto GaussianBasis *iBasis, *jBasis ;
            auto Real           dS, dT, dX, dY, dZ ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                iBasis = self->entries[i] ;
                i0     = Array1D_Item ( basisIndices, i   )      ;
                nI     = Array1D_Item ( basisIndices, i+1 ) - i0 ;
                for ( j = 0 ; j < i ; j++ ) /* . Diagonal terms are zero. */
                {
                    jBasis = self->entries[j] ;
                    j0     = Array1D_Item ( basisIndices, j   )      ;
                    nJ     = Array1D_Item ( basisIndices, j+1 ) - j0 ;
                    GaussianBasisIntegrals_Kinetic2OverlapD ( iBasis                                      ,
                                                              Coordinates3_RowPointer ( coordinates3, i ) ,
                                                              jBasis                                      ,
                                                              Coordinates3_RowPointer ( coordinates3, j ) ,
                                                              blockSX                                     ,
                                                              blockSY                                     ,
                                                              blockSZ                                     ,
                                                              blockTX                                     ,
                                                              blockTY                                     ,
                                                              blockTZ                                     ) ;
                    dX = dY = dZ = 0.0e+00 ;
                    for ( u = 0 ; u < nI ; u++ )
                    {
                        for ( v = 0 ; v < nJ ; v++ )
                        {
                            dS  =           SymmetricMatrix_Item ( oDensity, u+i0, v+j0 ) ;
                            dT  = 2.0e+00 * SymmetricMatrix_Item ( kDensity, u+i0, v+j0 ) ;
                            dX += ( dS * Array2D_Item ( blockSX, u, v ) + dT * Array2D_Item ( blockTX, u, v ) ) ;
                            dY += ( dS * Array2D_Item ( blockSY, u, v ) + dT * Array2D_Item ( blockTY, u, v ) ) ;
                            dZ += ( dS * Array2D_Item ( blockSZ, u, v ) + dT * Array2D_Item ( blockTZ, u, v ) ) ;
                        }
                    }
                    Coordinates3_IncrementRow ( gradients3, i, dX, dY, dZ ) ;
                    Coordinates3_DecrementRow ( gradients3, j, dX, dY, dZ ) ;
                }
            }
        }
        RealArray2D_Deallocate ( &blockSX ) ;
        RealArray2D_Deallocate ( &blockSY ) ;
        RealArray2D_Deallocate ( &blockSZ ) ;
        RealArray2D_Deallocate ( &blockTX ) ;
        RealArray2D_Deallocate ( &blockTY ) ;
        RealArray2D_Deallocate ( &blockTZ ) ;
    }
}
