/*==================================================================================================================================
! . Container integrals - 2 basis, 2 electrons, 1 nucleus/point.
!=================================================================================================================================*/

# include <math.h>
# include <stdlib.h>

# include "GaussianBasisContainerIntegrals_f2Cp1.h"
# include "GaussianBasisIntegrals_f2Cp1.h"
# include "Integer.h"
# include "IntegerUtilities.h"
# include "Real.h"
# include "RealArray2D.h"
# include "RealUtilities.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void _GetDensityFactors ( const Integer          i0      ,
                                 const Integer          nI      ,
                                 const Integer          j0      ,
                                 const Integer          nJ      ,
                                 const Boolean          iIsJ    , 
                                 const SymmetricMatrix *density ,
                                       RealArray2D     *dOneIJ  ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Electron-nuclear/point derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainerIntegrals_f2Cm1R1 ( const GaussianBasisContainer *self          ,
                                               const RealArray1D            *charges       ,
                                               const RealArray1D            *widthsE       ,
                                               const RealArray1D            *widthsN       ,
                                               const Coordinates3           *coordinates3  ,
                                               const Coordinates3           *coordinates3G ,
                                                     Selection              *selectionG    ,
                                               const SymmetricMatrix        *density       ,
                                                     Coordinates3           *gradients3    ,
                                                     Coordinates3           *gradients3G   ,
                                                     Status                 *status        )
{
    if ( ( self          != NULL ) &&
         ( charges       != NULL ) &&
         ( coordinates3  != NULL ) &&
         ( coordinates3G != NULL ) &&
         ( density       != NULL ) &&
         ( gradients3    != NULL ) &&
         ( gradients3G   != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer     *iWork, n, s2 ;
        auto Real         dRi[3], dRj[3], *rWork ;
        auto RealArray2D *block ;
        Selection_MakeFlags ( selectionG, Coordinates3_Rows ( coordinates3G ), status ) ;
        n     = GaussianBasisContainer_LargestBasis ( self, False ) ;
        block = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        n     = GaussianBasisContainer_LargestShell ( self, True ) ;
        s2    = n*n ;
        iWork = Integer_Allocate ( 6*s2, status ) ;
        rWork = Real_Allocate    ( 8*s2, status ) ;
        if ( Status_IsOK ( status ) )
        {
            auto Integer  i, i0, j, j0, nI, nJ ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                i0 = Array1D_Item ( self->centerFunctionPointers, i   )      ;
                nI = Array1D_Item ( self->centerFunctionPointers, i+1 ) - i0 ;
                for ( j = 0 ; j <= i ; j++ )
                {
                    j0 = Array1D_Item ( self->centerFunctionPointers, j   )      ;
                    nJ = Array1D_Item ( self->centerFunctionPointers, j+1 ) - j0 ;
                    _GetDensityFactors ( i0, nI, j0, nJ, ( i == j ), density, block ) ;
                    GaussianBasisIntegrals_f2Cm1R1 ( self->entries[i]                            ,
                                                     Coordinates3_RowPointer ( coordinates3, i ) ,
                                                     self->entries[j]                            ,
                                                     Coordinates3_RowPointer ( coordinates3, j ) ,
                                                     charges                                     ,
                                                     widthsE                                     ,
                                                     widthsN                                     ,
                                                     coordinates3G                               ,
                                                     selectionG                                  ,
                                                     block                                       ,
                                                     s2                                          ,
                                                     iWork                                       ,
                                                     rWork                                       ,
                                                     dRi                                         ,
                                                     dRj                                         ,
                                                     gradients3G                                 ) ;
                    Coordinates3_IncrementRow ( gradients3, i, dRi[0], dRi[1], dRi[2] ) ;
                    Coordinates3_IncrementRow ( gradients3, j, dRj[0], dRj[1], dRj[2] ) ;
                }
            }
        }
        Integer_Deallocate     ( &iWork ) ;
        Real_Deallocate        ( &rWork ) ;
        RealArray2D_Deallocate ( &block ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Electron-nuclear/point integrals.
! . Integrals should be appropriately initialized before entry to this function (often to the kinetic energy).
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainerIntegrals_f2Cm1V ( const GaussianBasisContainer *self              ,
                                              const RealArray1D            *charges           ,
                                              const RealArray1D            *widthsE           ,
                                              const RealArray1D            *widthsN           ,
                                              const Coordinates3           *coordinates3      ,
                                              const Coordinates3           *coordinates3G     ,
                                                    Selection              *selectionG        ,
                                                    SymmetricMatrix        *oneElectronMatrix ,
                                                    Status                 *status            )
{
    if ( ( self              != NULL ) &&
         ( charges           != NULL ) &&
         ( coordinates3      != NULL ) &&
         ( coordinates3G     != NULL ) &&
         ( oneElectronMatrix != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer     *iWork , n, s2 ;
        auto Real        *rWork ;
        auto RealArray2D *block ;
        Selection_MakeFlags ( selectionG, Coordinates3_Rows ( coordinates3G ), status ) ;
        n     = GaussianBasisContainer_LargestBasis ( self, False ) ;
        block = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        n     = GaussianBasisContainer_LargestShell ( self, True ) ;
        s2    = n*n ;
        iWork = Integer_Allocate ( 3*s2, status ) ;
        rWork = Real_Allocate    ( 4*s2, status ) ;
        if ( Status_IsOK ( status ) )
        {
            auto Integer  i, i0, j, j0, nI, nJ, u, v, vUpper ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                i0 = Array1D_Item ( self->centerFunctionPointers, i   )      ;
                nI = Array1D_Item ( self->centerFunctionPointers, i+1 ) - i0 ;
                for ( j = 0 ; j <= i ; j++ )
                {
                    j0 = Array1D_Item ( self->centerFunctionPointers, j   )      ;
                    nJ = Array1D_Item ( self->centerFunctionPointers, j+1 ) - j0 ;
                    GaussianBasisIntegrals_f2Cm1V ( self->entries[i]                            ,
                                                    Coordinates3_RowPointer ( coordinates3, i ) ,
                                                    self->entries[j]                            ,
                                                    Coordinates3_RowPointer ( coordinates3, j ) ,
                                                    charges                                     ,
                                                    widthsE                                     ,
                                                    widthsN                                     ,
                                                    coordinates3G                               ,
                                                    selectionG                                  ,
                                                    s2                                          ,
                                                    iWork                                       ,
                                                    rWork                                       ,
                                                    block                                       ) ;
                    vUpper = nJ ;
                    for ( u = 0 ; u < nI ; u++ )
                    {
                        if ( i == j ) vUpper = u + 1 ;
                        for ( v = 0 ; v < vUpper ; v++ ) { SymmetricMatrix_Item ( oneElectronMatrix, u+i0, v+j0 ) += Array2D_Item ( block, u, v ) ; }
                    }
                }
            }
        }
        Integer_Deallocate     ( &iWork ) ;
        Real_Deallocate        ( &rWork ) ;
        RealArray2D_Deallocate ( &block ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Electron-nuclear/point potentials.
! . Potentials should be appropriately initialized before entry to this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . CHECK SEPARATELY AT END! */
void GaussianBasisContainerIntegrals_f2Cp1V ( const GaussianBasisContainer *self          ,
                                              const RealArray1D            *widthsE       ,
                                              const RealArray1D            *widthsN       ,
                                              const Coordinates3           *coordinates3  ,
                                              const Coordinates3           *coordinates3G ,
                                                    Selection              *selectionG    ,
                                              const SymmetricMatrix        *density       ,
                                                    RealArray1D            *potentials    ,
                                                    Status                 *status        )
{
    if ( ( self          != NULL ) &&
         ( coordinates3  != NULL ) &&
         ( coordinates3G != NULL ) &&
         ( density       != NULL ) &&
         ( potentials    != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer     *iWork , n, s2 ;
        auto Real        *rWork ;
        auto RealArray2D *block ;
        Selection_MakeFlags ( selectionG, Coordinates3_Rows ( coordinates3G ), status ) ;
        n     = GaussianBasisContainer_LargestBasis ( self, False ) ;
        block = RealArray2D_AllocateWithExtents ( n, n, status ) ;
        n     = GaussianBasisContainer_LargestShell ( self, True ) ;
        s2    = n*n ;
        iWork = Integer_Allocate ( 3*s2, status ) ;
        rWork = Real_Allocate    ( 3*s2, status ) ;
        if ( Status_IsOK ( status ) )
        {
            auto Integer  i, i0, j, j0, nI, nJ ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                i0 = Array1D_Item ( self->centerFunctionPointers, i   )      ;
                nI = Array1D_Item ( self->centerFunctionPointers, i+1 ) - i0 ;
                for ( j = 0 ; j <= i ; j++ )
                {
                    j0 = Array1D_Item ( self->centerFunctionPointers, j   )      ;
                    nJ = Array1D_Item ( self->centerFunctionPointers, j+1 ) - j0 ;
                    _GetDensityFactors ( i0, nI, j0, nJ, ( i == j ), density, block ) ;
                    GaussianBasisIntegrals_f2Cp1V ( self->entries[i]                            ,
                                                    Coordinates3_RowPointer ( coordinates3, i ) ,
                                                    self->entries[j]                            ,
                                                    Coordinates3_RowPointer ( coordinates3, j ) ,
                                                    widthsE                                     ,
                                                    widthsN                                     ,
                                                    coordinates3G                               ,
                                                    selectionG                                  ,
                                                    block                                       ,
                                                    s2                                          ,
                                                    iWork                                       ,
                                                    rWork                                       ,
                                                    potentials                                  ) ;
                }
            }
        }
        Integer_Deallocate     ( &iWork ) ;
        Real_Deallocate        ( &rWork ) ;
        RealArray2D_Deallocate ( &block ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get density factors.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void _GetDensityFactors ( const Integer          i0      ,
                                 const Integer          nI      ,
                                 const Integer          j0      ,
                                 const Integer          nJ      ,
                                 const Boolean          iIsJ    , 
                                 const SymmetricMatrix *density ,
                                       RealArray2D     *dOneIJ  )
{
    auto Integer  u, v ;
    if ( iIsJ )
    {
        for ( u = 0 ; u < nI ; u++ )
        {
            for ( v = 0   ; v <= u ; v++ ) { Array2D_Item ( dOneIJ, u, v ) = SymmetricMatrix_Item ( density, u+i0, v+j0 ) ; }
            for ( v = u+1 ; v < nJ ; v++ ) { Array2D_Item ( dOneIJ, u, v ) = SymmetricMatrix_Item ( density, v+j0, u+i0 ) ; }
        }
    }
    else
    {
        for ( u = 0 ; u < nI ; u++ )
        {
            for ( v = 0 ; v < nJ ; v++ ) { Array2D_Item ( dOneIJ, u, v ) = SymmetricMatrix_Item ( density, u+i0, v+j0 ) ; }
        }
    }
}
