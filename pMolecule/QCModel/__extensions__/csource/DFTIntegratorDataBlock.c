/*==================================================================================================================================
! . This module defines a block data structure that is needed for DFT integration.
!=================================================================================================================================*/

# include "DFTIntegratorDataBlock.h"
# include "Memory.h"
# include "NumericalMacros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Accumulation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DFTIntegratorDataBlock_Accumulate ( DFTIntegratorDataBlock *self )
{
    if ( ( self != NULL ) && ( self->hasLocalData ) )
    {
        RealArray1D_Add ( self->eXC           , 1.0e+00 , self->localEXC           , NULL ) ;
        RealArray2D_Add ( self->vLaplacianRho , 1.0e+00 , self->localVLaplacianRho , NULL ) ;
        RealArray2D_Add ( self->vRho          , 1.0e+00 , self->localVRho          , NULL ) ;
        RealArray2D_Add ( self->vSigma        , 1.0e+00 , self->localVSigma        , NULL ) ;
        RealArray2D_Add ( self->vTau          , 1.0e+00 , self->localVTau          , NULL ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
DFTIntegratorDataBlock *DFTIntegratorDataBlock_Allocate ( const Integer numberOfFunctionals ,
                                                          const Integer numberOfPoints      ,
                                                          const Boolean hasSigma            ,
                                                          const Boolean hasLaplacian        ,
                                                          const Boolean hasTau              ,
                                                          const Boolean isSpinRestricted    ,
                                                                Status  *status             )
{
    DFTIntegratorDataBlock *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( DFTIntegratorDataBlock ) ;
        if ( self != NULL )
        {
            auto Integer n ;
            n = Maximum ( numberOfPoints, 0 ) ;
            self->hasLocalData   = ( numberOfFunctionals > 1 ) ;
            self->numberOfPoints = n ;
            /* . Basic initialization. */
            self->eXC                = NULL ;
            self->localEXC           = NULL ;
            self->dRhoX              = NULL ;
            self->dRhoY              = NULL ;
            self->dRhoZ              = NULL ;
            self->localVLaplacianRho = NULL ;
            self->localVRho          = NULL ;
            self->localVSigma        = NULL ;
            self->localVTau          = NULL ;
            self->laplacianRho       = NULL ;
            self->rho                = NULL ;
            self->sigma              = NULL ;
            self->tau                = NULL ;
            self->vLaplacianRho      = NULL ;
            self->vRho               = NULL ;
            self->vSigma             = NULL ;
            self->vTau               = NULL ;
            /* . Allocation. */ 
            if ( n > 0 )
            {
                auto Boolean isOK ;
                auto Integer c, d ;
                /* . Determine array sizes. */
                if ( isSpinRestricted ) { c = 1 ; d = 1 ; }
                else                    { c = 2 ; d = 3 ; }
                /* . Allocation. */
                self->eXC    = RealArray1D_AllocateWithExtent  ( n,    status ) ;
                self->rho    = RealArray2D_AllocateWithExtents ( n, c, status ) ;
                self->vRho   = RealArray2D_AllocateWithExtents ( n, c, status ) ;
                if ( self->hasLocalData )
                {
                    self->localEXC  = RealArray1D_AllocateWithExtent  ( n,    status ) ;
                    self->localVRho = RealArray2D_AllocateWithExtents ( n, c, status ) ;
                }
                else { self->localEXC = self->eXC ; self->localVRho = self->vRho ; }
                isOK = ( ( self->eXC       != NULL ) &&
                         ( self->localEXC  != NULL ) &&
                         ( self->localVRho != NULL ) &&
                         ( self->rho       != NULL ) &&
                         ( self->vRho      != NULL ) ) ;
                if ( hasSigma )
                {
                    self->dRhoX       = RealArray2D_AllocateWithExtents ( n, c, status ) ;
                    self->dRhoY       = RealArray2D_AllocateWithExtents ( n, c, status ) ;
                    self->dRhoZ       = RealArray2D_AllocateWithExtents ( n, c, status ) ;
                    self->sigma       = RealArray2D_AllocateWithExtents ( n, d, status ) ;
                    self->vSigma      = RealArray2D_AllocateWithExtents ( n, d, status ) ;
                    if ( self->hasLocalData ) self->localVSigma = RealArray2D_AllocateWithExtents ( n, d, status ) ;
                    else                      self->localVSigma = self->vSigma ;
                    isOK = isOK && ( ( self->dRhoX       != NULL ) &&
                                     ( self->dRhoY       != NULL ) &&
                                     ( self->dRhoZ       != NULL ) &&
                                     ( self->localVSigma != NULL ) &&
                                     ( self->sigma       != NULL ) &&
                                     ( self->vSigma      != NULL ) ) ;
                }
                if ( hasLaplacian )
                {
                    self->laplacianRho  = RealArray2D_AllocateWithExtents ( n, c, status ) ;
                    self->vLaplacianRho = RealArray2D_AllocateWithExtents ( n, c, status ) ;
                    if ( self->hasLocalData ) self->localVLaplacianRho = RealArray2D_AllocateWithExtents ( n, c, status ) ;
                    else                      self->localVLaplacianRho = self->vLaplacianRho ;
                    isOK = isOK && ( ( self->localVLaplacianRho != NULL ) &&
                                     ( self->laplacianRho       != NULL ) &&
                                     ( self->vLaplacianRho      != NULL ) ) ;
                }
                if ( hasTau )
                {
                    self->tau  = RealArray2D_AllocateWithExtents ( n, c, status ) ;
                    self->vTau = RealArray2D_AllocateWithExtents ( n, c, status ) ;
                    if ( self->hasLocalData ) self->localVTau = RealArray2D_AllocateWithExtents ( n, c, status ) ;
                    else                      self->localVTau = self->vTau ;
                    isOK = isOK && ( ( self->localVTau != NULL ) &&
                                     ( self->tau       != NULL ) &&
                                     ( self->vTau      != NULL ) ) ;
                }
                if ( isOK )
                {   
                    DFTIntegratorDataBlock_InitializeView ( self, 0, &(self->viewP) ) ;
                    if ( ! isSpinRestricted )
                    {
                        DFTIntegratorDataBlock_InitializeView ( self, 1, &(self->viewQ) ) ;
                        RealArray2D_ColumnView ( self->sigma , 1, False, &(self->sigmaPQ ), NULL ) ;
                        RealArray2D_ColumnView ( self->vSigma, 1, False, &(self->vSigmaPQ), NULL ) ;
                    }
                }
                else DFTIntegratorDataBlock_Deallocate ( &self ) ;
            }
        }
    }
    if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DFTIntegratorDataBlock_Deallocate ( DFTIntegratorDataBlock **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        RealArray1D_Deallocate ( &((*self)->eXC          ) ) ;
        RealArray2D_Deallocate ( &((*self)->dRhoX        ) ) ;
        RealArray2D_Deallocate ( &((*self)->dRhoY        ) ) ;
        RealArray2D_Deallocate ( &((*self)->dRhoZ        ) ) ;
        RealArray2D_Deallocate ( &((*self)->laplacianRho ) ) ;
        RealArray2D_Deallocate ( &((*self)->rho          ) ) ;
        RealArray2D_Deallocate ( &((*self)->sigma        ) ) ;
        RealArray2D_Deallocate ( &((*self)->tau          ) ) ;
        RealArray2D_Deallocate ( &((*self)->vLaplacianRho) ) ;
        RealArray2D_Deallocate ( &((*self)->vRho         ) ) ;
        RealArray2D_Deallocate ( &((*self)->vSigma       ) ) ;
        RealArray2D_Deallocate ( &((*self)->vTau         ) ) ;
        if ( (*self)->hasLocalData )
        {
            RealArray1D_Deallocate ( &((*self)->localEXC          ) ) ;
            RealArray2D_Deallocate ( &((*self)->localVLaplacianRho) ) ;
            RealArray2D_Deallocate ( &((*self)->localVRho         ) ) ;
            RealArray2D_Deallocate ( &((*self)->localVSigma       ) ) ;
            RealArray2D_Deallocate ( &((*self)->localVTau         ) ) ;
        }
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DFTIntegratorDataBlock_Initialize ( DFTIntegratorDataBlock *self )
{
    if ( ( self != NULL ) && ( self->hasLocalData ) )
    {
        RealArray1D_Set ( self->eXC          , 0.0e+00 ) ;
        RealArray2D_Set ( self->vLaplacianRho, 0.0e+00 ) ;
        RealArray2D_Set ( self->vRho         , 0.0e+00 ) ;
        RealArray2D_Set ( self->vSigma       , 0.0e+00 ) ;
        RealArray2D_Set ( self->vTau         , 0.0e+00 ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . View initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DFTIntegratorDataBlock_InitializeView ( DFTIntegratorDataBlock *self, const Integer c, DFTIntegratorDataBlockView *view )
{
    if ( ( self != NULL ) && ( view != NULL ) )
    {
        auto Integer d ;
        d = ( c == 0 ? 0 : 2 ) ;
        RealArray2D_ColumnView ( self->dRhoX        , c, False, &(view->dRhoX        ), NULL ) ;
        RealArray2D_ColumnView ( self->dRhoY        , c, False, &(view->dRhoY        ), NULL ) ;
        RealArray2D_ColumnView ( self->dRhoZ        , c, False, &(view->dRhoZ        ), NULL ) ;
        RealArray2D_ColumnView ( self->laplacianRho , c, False, &(view->laplacianRho ), NULL ) ;
        RealArray2D_ColumnView ( self->rho          , c, False, &(view->rho          ), NULL ) ;
        RealArray2D_ColumnView ( self->sigma        , d, False, &(view->sigma        ), NULL ) ;
        RealArray2D_ColumnView ( self->tau          , c, False, &(view->tau          ), NULL ) ;
        RealArray2D_ColumnView ( self->vLaplacianRho, c, False, &(view->vLaplacianRho), NULL ) ;
        RealArray2D_ColumnView ( self->vRho         , c, False, &(view->vRho         ), NULL ) ;
        RealArray2D_ColumnView ( self->vSigma       , d, False, &(view->vSigma       ), NULL ) ;
        RealArray2D_ColumnView ( self->vTau         , c, False, &(view->vTau         ), NULL ) ;
    }
}
