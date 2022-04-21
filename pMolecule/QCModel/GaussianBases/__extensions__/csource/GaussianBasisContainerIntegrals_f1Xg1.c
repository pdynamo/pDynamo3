/*==================================================================================================================================
! . Container integrals - 2 basis, 1 electrons, 0 nuclei/points.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "GaussianBasisContainerIntegrals_f1Xg1.h"
# include "GaussianBasisIntegrals_f1Xg1.h"
# include "Integer.h"
# include "IntegerUtilities.h"
# include "RealArray2D.h"
# include "RealUtilities.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Anti-Coulomb integrals.
! . Integrals is overwritten by this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainerIntegrals_f1Af1i ( const GaussianBasisContainer *self         ,
                                              const Coordinates3           *coordinates3 ,
                                                    SymmetricMatrix        *integrals    ,
                                                    Status                 *status       )
{
    SymmetricMatrix_Set ( integrals, 0.0e+00 ) ;
    if ( ( self         != NULL ) &&
         ( coordinates3 != NULL ) &&
         ( integrals    != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer     *iWork , n, s2 ;
        auto Real        *rWork ;
        auto RealArray2D *block ;
        n     = GaussianBasisContainer_LargestBasis ( self, False ) ;
        block = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        n     = GaussianBasisContainer_LargestShell ( self, True ) ;
        s2    = n*n ;
        iWork = Integer_Allocate ( 6*s2, status ) ;
        rWork = Real_Allocate    ( 3*s2, status ) ;
        if ( Status_IsOK ( status ) )
        {
            auto Integer        i, i0, j, j0, nI, nJ, u, v, vUpper ;
            auto GaussianBasis *iBasis, *jBasis ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                iBasis = self->entries[i] ;
                i0     = Array1D_Item ( self->centerFunctionPointers, i   )      ;
                nI     = Array1D_Item ( self->centerFunctionPointers, i+1 ) - i0 ;
                for ( j = 0 ; j <= i ; j++ )
                {
                    jBasis = self->entries[j] ;
                    j0     = Array1D_Item ( self->centerFunctionPointers, j   )      ;
                    nJ     = Array1D_Item ( self->centerFunctionPointers, j+1 ) - j0 ;
                    GaussianBasisIntegrals_f1Ag1i ( iBasis                                      ,
                                                    Coordinates3_RowPointer ( coordinates3, i ) ,
                                                    jBasis                                      ,
                                                    Coordinates3_RowPointer ( coordinates3, j ) ,
                                                    s2                                          ,
                                                    iWork                                       ,
                                                    rWork                                       ,
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
printf ( "\nAnti-Coulomb 2D Integrals (%d):\n", self->capacity ) ;
SymmetricMatrix_Print ( integrals ) ;
*/
        }
        Integer_Deallocate     ( &iWork ) ;
        Real_Deallocate        ( &rWork ) ;
        RealArray2D_Deallocate ( &block ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Coulomb integrals.
! . Integrals is overwritten by this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainerIntegrals_f1Cf1i ( const GaussianBasisContainer *self         ,
                                              const Coordinates3           *coordinates3 ,
                                                    SymmetricMatrix        *integrals    ,
                                                    Status                 *status       )
{
    SymmetricMatrix_Set ( integrals, 0.0e+00 ) ;
    if ( ( self         != NULL ) &&
         ( coordinates3 != NULL ) &&
         ( integrals    != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer     *iWork , n, s2 ;
        auto Real        *rWork ;
        auto RealArray2D *block ;
        n     = GaussianBasisContainer_LargestBasis ( self, False ) ;
        block = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        n     = GaussianBasisContainer_LargestShell ( self, True ) ;
        s2    = n*n ;
        iWork = Integer_Allocate ( 3*s2, status ) ;
        rWork = Real_Allocate    ( 3*s2, status ) ;
        if ( Status_IsOK ( status ) )
        {
            auto Integer        i, i0, j, j0, nI, nJ, u, v, vUpper ;
            auto GaussianBasis *iBasis, *jBasis ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                iBasis = self->entries[i] ;
                i0     = Array1D_Item ( self->centerFunctionPointers, i   )      ;
                nI     = Array1D_Item ( self->centerFunctionPointers, i+1 ) - i0 ;
                for ( j = 0 ; j <= i ; j++ )
                {
                    jBasis = self->entries[j] ;
                    j0     = Array1D_Item ( self->centerFunctionPointers, j   )      ;
                    nJ     = Array1D_Item ( self->centerFunctionPointers, j+1 ) - j0 ;
                    GaussianBasisIntegrals_f1Cg1i ( iBasis                                      ,
                                                    Coordinates3_RowPointer ( coordinates3, i ) ,
                                                    jBasis                                      ,
                                                    Coordinates3_RowPointer ( coordinates3, j ) ,
                                                    s2                                          ,
                                                    iWork                                       ,
                                                    rWork                                       ,
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
        Integer_Deallocate     ( &iWork ) ;
        Real_Deallocate        ( &rWork ) ;
        RealArray2D_Deallocate ( &block ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Dipole integrals.
! . The dipole matrices are overwritten by this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainerIntegrals_f1Df1i ( const GaussianBasisContainer *self         ,
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
         ( coordinates3 != NULL ) &&
         ( dipoleX      != NULL ) &&
         ( dipoleY      != NULL ) &&
         ( dipoleZ      != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer      n, s2 ;
        auto Real        *rWork ;
        auto RealArray2D *blockX, *blockY, *blockZ ;
        auto Vector3     *origin ;
        n      = GaussianBasisContainer_LargestBasis ( self, False ) ;
        blockX = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        blockY = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        blockZ = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        n      = GaussianBasisContainer_LargestShell ( self, True ) ;
        s2     = n*n ;
        rWork  = Real_Allocate ( 4*s2, status ) ;
        if ( center == NULL ) { origin = Vector3_Allocate ( ) ; Vector3_Set ( origin, 0.0e+00 ) ; }
        else                    origin = center ;
        if ( ( origin != NULL ) && Status_IsOK ( status ) )
        {
            auto Integer        i, i0, j, j0, nI, nJ, u, v, vUpper ;
            auto GaussianBasis *iBasis, *jBasis ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                iBasis = self->entries[i] ;
                i0     = Array1D_Item ( self->centerFunctionPointers, i   )      ;
                nI     = Array1D_Item ( self->centerFunctionPointers, i+1 ) - i0 ;
                for ( j = 0 ; j <= i ; j++ )
                {
                    jBasis = self->entries[j] ;
                    j0     = Array1D_Item ( self->centerFunctionPointers, j   )      ;
                    nJ     = Array1D_Item ( self->centerFunctionPointers, j+1 ) - j0 ;
                    GaussianBasisIntegrals_f1Df1i ( iBasis                                      ,
                                                    Coordinates3_RowPointer ( coordinates3, i ) ,
                                                    jBasis                                      ,
                                                    Coordinates3_RowPointer ( coordinates3, j ) ,
                                                    Vector3_Data ( origin )                     ,
                                                    s2                                          ,
                                                    rWork                                       ,
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
        Real_Deallocate        ( &rWork  ) ;
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
void GaussianBasisContainerIntegrals_f1KOf1i ( const GaussianBasisContainer *self         ,
                                               const Coordinates3           *coordinates3 ,
                                                     SymmetricMatrix        *kinetic      ,
                                                     SymmetricMatrix        *overlap      ,
                                                     Status                 *status       )
{
    if ( ( self         != NULL ) &&
         ( coordinates3 != NULL ) &&
         ( kinetic      != NULL ) &&
         ( overlap      != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer      n, s2 ;
        auto Real        *rWork ;
        auto RealArray2D *blockS, *blockT ;
        n      = GaussianBasisContainer_LargestBasis ( self, False ) ;
        blockS = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        blockT = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        n      = GaussianBasisContainer_LargestShell ( self, True ) ;
        s2     = n*n ;
        rWork  = Real_Allocate ( 3*s2, status ) ;
        if ( Status_IsOK ( status ) )
        {
            auto Integer        i, i0, j, j0, nI, nJ, u, v, vUpper ;
            auto GaussianBasis *iBasis, *jBasis ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                iBasis = self->entries[i] ;
                i0     = Array1D_Item ( self->centerFunctionPointers, i   )      ;
                nI     = Array1D_Item ( self->centerFunctionPointers, i+1 ) - i0 ;
                for ( j = 0 ; j <= i ; j++ )
                {
                    jBasis = self->entries[j] ;
                    j0     = Array1D_Item ( self->centerFunctionPointers, j   )      ;
                    nJ     = Array1D_Item ( self->centerFunctionPointers, j+1 ) - j0 ;
                    GaussianBasisIntegrals_f1KOg1i ( iBasis                                      ,
                                                     Coordinates3_RowPointer ( coordinates3, i ) ,
                                                     jBasis                                      ,
                                                     Coordinates3_RowPointer ( coordinates3, j ) ,
                                                     s2                                          ,
                                                     rWork                                       ,
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
        Real_Deallocate        ( &rWork  ) ;
        RealArray2D_Deallocate ( &blockS ) ;
        RealArray2D_Deallocate ( &blockT ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Kinetic and overlap derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainerIntegrals_f1KOf1R1 ( const GaussianBasisContainer *self         ,
                                                const Coordinates3           *coordinates3 ,
                                                const SymmetricMatrix        *kDensity     ,
                                                const SymmetricMatrix        *oDensity     ,
                                                      Coordinates3           *gradients3   ,
                                                      Status                 *status       )
{
    if ( ( self         != NULL ) &&
         ( coordinates3 != NULL ) &&
         ( kDensity     != NULL ) &&
         ( oDensity     != NULL ) &&
         ( gradients3   != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer      n, s2 ;
        auto Real        *rWork ;
        auto RealArray2D *blockSX, *blockSY, *blockSZ ,
                         *blockTX, *blockTY, *blockTZ ;
        n      = GaussianBasisContainer_LargestBasis ( self, False ) ;
        blockSX = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        blockSY = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        blockSZ = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        blockTX = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        blockTY = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        blockTZ = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        n       = GaussianBasisContainer_LargestShell ( self, True ) ;
        s2      = n*n ;
        rWork   = Real_Allocate ( 7*s2, status ) ;
        if ( Status_IsOK ( status ) )
        {
            auto Integer        i, i0, j, j0, nI, nJ, u, v ;
            auto GaussianBasis *iBasis, *jBasis ;
            auto Real           dS, dT, dX, dY, dZ ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                iBasis = self->entries[i] ;
                i0     = Array1D_Item ( self->centerFunctionPointers, i   )      ;
                nI     = Array1D_Item ( self->centerFunctionPointers, i+1 ) - i0 ;
                for ( j = 0 ; j < i ; j++ ) /* . Diagonal terms are zero. */
                {
                    jBasis = self->entries[j] ;
                    j0     = Array1D_Item ( self->centerFunctionPointers, j   )      ;
                    nJ     = Array1D_Item ( self->centerFunctionPointers, j+1 ) - j0 ;
                    GaussianBasisIntegrals_f1KOg1r1 ( iBasis                                      ,
                                                      Coordinates3_RowPointer ( coordinates3, i ) ,
                                                      jBasis                                      ,
                                                      Coordinates3_RowPointer ( coordinates3, j ) ,
                                                      s2                                          ,
                                                      rWork                                       ,
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
        Real_Deallocate        ( &rWork   ) ;
        RealArray2D_Deallocate ( &blockSX ) ;
        RealArray2D_Deallocate ( &blockSY ) ;
        RealArray2D_Deallocate ( &blockSZ ) ;
        RealArray2D_Deallocate ( &blockTX ) ;
        RealArray2D_Deallocate ( &blockTY ) ;
        RealArray2D_Deallocate ( &blockTZ ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Overlap integrals.
! . Overlap is overwritten by this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainerIntegrals_f1Of1i ( const GaussianBasisContainer *self         ,
                                              const Coordinates3           *coordinates3 ,
                                                    SymmetricMatrix        *overlap      ,
                                                    Status                 *status       )
{
    SymmetricMatrix_Set ( overlap, 0.0e+00 ) ;
    if ( ( self         != NULL ) &&
         ( coordinates3 != NULL ) &&
         ( overlap      != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer      n, s2 ;
        auto Real        *rWork ;
        auto RealArray2D *block ;
        n     = GaussianBasisContainer_LargestBasis ( self, False ) ;
        block = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        n     = GaussianBasisContainer_LargestShell ( self, True ) ;
        s2    = n*n ;
        rWork = Real_Allocate ( 2*s2, status ) ;
        if ( Status_IsOK ( status ) )
        {
            auto Integer        i, i0, j, j0, nI, nJ, u, v, vUpper ;
            auto GaussianBasis *iBasis, *jBasis ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                iBasis = self->entries[i] ;
                i0     = Array1D_Item ( self->centerFunctionPointers, i   )      ;
                nI     = Array1D_Item ( self->centerFunctionPointers, i+1 ) - i0 ;
                for ( j = 0 ; j <= i ; j++ )
                {
                    jBasis = self->entries[j] ;
                    j0     = Array1D_Item ( self->centerFunctionPointers, j   )      ;
                    nJ     = Array1D_Item ( self->centerFunctionPointers, j+1 ) - j0 ;
                    GaussianBasisIntegrals_f1Og1i ( iBasis                                      ,
                                                    Coordinates3_RowPointer ( coordinates3, i ) ,
                                                    jBasis                                      ,
                                                    Coordinates3_RowPointer ( coordinates3, j ) ,
                                                    s2                                          ,
                                                    rWork                                       ,
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
        Real_Deallocate        ( &rWork ) ;
        RealArray2D_Deallocate ( &block ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Quadrupole integrals.
! . The quadrupole matrices are overwritten by this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainerIntegrals_f1Qf1i ( const GaussianBasisContainer *self         ,
                                              const Coordinates3           *coordinates3 ,
                                                    Vector3                *center       ,
                                                    SymmetricMatrix        *qXX          ,
                                                    SymmetricMatrix        *qYY          ,
                                                    SymmetricMatrix        *qZZ          ,
                                                    SymmetricMatrix        *qXY          ,
                                                    SymmetricMatrix        *qXZ          ,
                                                    SymmetricMatrix        *qYZ          ,
                                                    Status                 *status       )
{
    SymmetricMatrix_Set ( qXX, 0.0e+00 ) ;
    SymmetricMatrix_Set ( qYY, 0.0e+00 ) ;
    SymmetricMatrix_Set ( qZZ, 0.0e+00 ) ;
    SymmetricMatrix_Set ( qXY, 0.0e+00 ) ;
    SymmetricMatrix_Set ( qXZ, 0.0e+00 ) ;
    SymmetricMatrix_Set ( qYZ, 0.0e+00 ) ;
    if ( ( self         != NULL ) &&
         ( coordinates3 != NULL ) &&
         ( qXX          != NULL ) &&
         ( qYY          != NULL ) &&
         ( qZZ          != NULL ) &&
         ( qXY          != NULL ) &&
         ( qXZ          != NULL ) &&
         ( qYZ          != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer      n, s2 ;
        auto Real        *rWork ;
        auto RealArray2D *blockXX, *blockXY, *blockXZ, *blockYY, *blockYZ, *blockZZ ;
        auto Vector3     *origin ;
        n       = GaussianBasisContainer_LargestBasis ( self, False ) ;
        blockXX = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        blockYY = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        blockZZ = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        blockXY = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        blockXZ = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        blockYZ = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        n       = GaussianBasisContainer_LargestShell ( self, True ) ;
        s2      = n*n ;
        rWork   = Real_Allocate ( 7*s2, status ) ;
        if ( center == NULL ) { origin = Vector3_Allocate ( ) ; Vector3_Set ( origin, 0.0e+00 ) ; }
        else                    origin = center ;
        if ( ( origin != NULL ) && Status_IsOK ( status ) )
        {
            auto Integer        i, i0, j, j0, nI, nJ, u, v, vUpper ;
            auto GaussianBasis *iBasis, *jBasis ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                iBasis = self->entries[i] ;
                i0     = Array1D_Item ( self->centerFunctionPointers, i   )      ;
                nI     = Array1D_Item ( self->centerFunctionPointers, i+1 ) - i0 ;
                for ( j = 0 ; j <= i ; j++ )
                {
                    jBasis = self->entries[j] ;
                    j0     = Array1D_Item ( self->centerFunctionPointers, j   )      ;
                    nJ     = Array1D_Item ( self->centerFunctionPointers, j+1 ) - j0 ;
                    GaussianBasisIntegrals_f1Qf1i ( iBasis                                      ,
                                                    Coordinates3_RowPointer ( coordinates3, i ) ,
                                                    jBasis                                      ,
                                                    Coordinates3_RowPointer ( coordinates3, j ) ,
                                                    Vector3_Data ( origin )                     ,
                                                    s2                                          ,
                                                    rWork                                       ,
                                                    blockXX                                     ,
                                                    blockYY                                     ,
                                                    blockZZ                                     ,
                                                    blockXY                                     ,
                                                    blockXZ                                     ,
                                                    blockYZ                                     ) ;
                    vUpper = nJ ;
                    for ( u = 0 ; u < nI ; u++ )
                    {
                        if ( i == j ) vUpper = u + 1 ;
                        for ( v = 0 ; v < vUpper ; v++ )
                        {
                            SymmetricMatrix_Item ( qXX, u+i0, v+j0 ) = Array2D_Item ( blockXX, u, v ) ;
                            SymmetricMatrix_Item ( qYY, u+i0, v+j0 ) = Array2D_Item ( blockYY, u, v ) ;
                            SymmetricMatrix_Item ( qZZ, u+i0, v+j0 ) = Array2D_Item ( blockZZ, u, v ) ;
                            SymmetricMatrix_Item ( qXY, u+i0, v+j0 ) = Array2D_Item ( blockXY, u, v ) ;
                            SymmetricMatrix_Item ( qXZ, u+i0, v+j0 ) = Array2D_Item ( blockXZ, u, v ) ;
                            SymmetricMatrix_Item ( qYZ, u+i0, v+j0 ) = Array2D_Item ( blockYZ, u, v ) ;
                        }
                    }
                }
            }
        }
        Real_Deallocate        ( &rWork   ) ;
        RealArray2D_Deallocate ( &blockXX ) ;
        RealArray2D_Deallocate ( &blockYY ) ;
        RealArray2D_Deallocate ( &blockZZ ) ;
        RealArray2D_Deallocate ( &blockXY ) ;
        RealArray2D_Deallocate ( &blockXZ ) ;
        RealArray2D_Deallocate ( &blockYZ ) ;
        if ( center == NULL ) Vector3_Deallocate ( &origin ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Integral derivatives for density fitting.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainerIntegrals_f1Xf1R1 ( const GaussianBasisContainer *self         ,
                                               const Coordinates3           *coordinates3 ,
                                               const RealArray1D            *aVector      ,
                                               const RealArray1D            *xVector      ,
                                               const GaussianBasisOperator   operator     ,
                                                     Coordinates3           *gradients3   ,
                                                     Status                 *status       )
{
    if ( ( self         != NULL ) &&
         ( coordinates3 != NULL ) &&
         ( aVector      != NULL ) &&
         ( xVector      != NULL ) &&
         ( gradients3   != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer     *iWork , n, nI = 0, nR = 0, s2 ;
        auto Real        *rWork ;
        auto RealArray2D *aX, *blockX, *blockY, *blockZ ;
        n      = GaussianBasisContainer_LargestBasis ( self, False ) ;
        aX     = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        blockX = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        blockY = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        blockZ = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        n      = GaussianBasisContainer_LargestShell ( self, True ) ;
        s2     = n*n ;
             if ( operator == GaussianBasisOperator_AntiCoulomb ) { nI = 9*s2 ; nR = 5*s2 ; }
        else if ( operator == GaussianBasisOperator_Coulomb     ) { nI = 3*s2 ; nR = 5*s2 ; }
        else if ( operator == GaussianBasisOperator_Overlap     ) { nI = 0    ; nR = 4*s2 ; }
        if ( nI > 0 ) iWork = Integer_Allocate ( nI, status ) ;
        else          iWork = NULL ;
        rWork = Real_Allocate ( nR, status ) ;
        if ( Status_IsOK ( status ) )
        {
            auto Integer        i, i0, j, j0, nI, nJ, u, v ;
            auto GaussianBasis *iBasis, *jBasis ;
            auto Real           aU, aV, d, dX, dY, dZ, xU, xV ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                iBasis = self->entries[i] ;
                i0     = Array1D_Item ( self->centerFunctionPointers, i   )      ;
                nI     = Array1D_Item ( self->centerFunctionPointers, i+1 ) - i0 ;
                for ( j = 0 ; j < i ; j++ ) /* . Diagonal terms are zero. */
                {
                    jBasis = self->entries[j] ;
                    j0     = Array1D_Item ( self->centerFunctionPointers, j   )      ;
                    nJ     = Array1D_Item ( self->centerFunctionPointers, j+1 ) - j0 ;
                    if ( operator == GaussianBasisOperator_AntiCoulomb )
                    {
                        GaussianBasisIntegrals_f1Ag1r1 ( iBasis                                      ,
                                                         Coordinates3_RowPointer ( coordinates3, i ) ,
                                                         jBasis                                      ,
                                                         Coordinates3_RowPointer ( coordinates3, j ) ,
                                                         s2                                          ,
                                                         iWork                                       ,
                                                         rWork                                       ,
                                                         blockX                                      ,
                                                         blockY                                      ,
                                                         blockZ                                      ) ;
                    }
                    else if ( operator == GaussianBasisOperator_Coulomb )
                    {
                        GaussianBasisIntegrals_f1Cg1r1 ( iBasis                                      ,
                                                         Coordinates3_RowPointer ( coordinates3, i ) ,
                                                         jBasis                                      ,
                                                         Coordinates3_RowPointer ( coordinates3, j ) ,
                                                         s2                                          ,
                                                         iWork                                       ,
                                                         rWork                                       ,
                                                         blockX                                      ,
                                                         blockY                                      ,
                                                         blockZ                                      ) ;
                    }
                    else if ( operator == GaussianBasisOperator_Overlap )
                    {
                        GaussianBasisIntegrals_f1Og1r1 ( iBasis                                      ,
                                                         Coordinates3_RowPointer ( coordinates3, i ) ,
                                                         jBasis                                      ,
                                                         Coordinates3_RowPointer ( coordinates3, j ) ,
                                                         s2                                          ,
                                                         rWork                                       ,
                                                         blockX                                      ,
                                                         blockY                                      ,
                                                         blockZ                                      ) ;
                    }
                    for ( u = 0 ; u < nI ; u++ )
                    {
                        aU = Array1D_Item ( aVector, u+i0 ) ;
                        xU = Array1D_Item ( xVector, u+i0 ) ;
                        for ( v = 0 ; v < nJ ; v++ )
                        {
                            aV = Array1D_Item ( aVector, v+j0 ) ;
                            xV = Array1D_Item ( xVector, v+j0 ) ;
                            Array2D_Item ( aX, u, v ) = - 0.5e+00 * ( aU * xV + aV * xU ) ;
                        }
                    }
                    dX = dY = dZ = 0.0e+00 ;
                    for ( u = 0 ; u < nI ; u++ )
                    {
                        for ( v = 0 ; v < nJ ; v++ )
                        {
                            d   =     Array2D_Item ( aX    , u, v ) ;
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
        Integer_Deallocate     ( &iWork  ) ;
        Real_Deallocate        ( &rWork  ) ;
        RealArray2D_Deallocate ( &aX     ) ;
        RealArray2D_Deallocate ( &blockX ) ;
        RealArray2D_Deallocate ( &blockY ) ;
        RealArray2D_Deallocate ( &blockZ ) ;
    }
}
