/*==================================================================================================================================
! . This module defines the DFT functional model. It is an interface to the libxc library.
!=================================================================================================================================*/

# include "Array_Macros.h"
# include "DFTFunctionalModel.h"
# include "Memory.h"
# include "NumericalMacros.h"


/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
DFTFunctionalModel *DFTFunctionalModel_Allocate ( const Integer numberOfFunctionals, Status *status )
{
    DFTFunctionalModel *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( DFTFunctionalModel ) ;
        if ( self != NULL )
        {
            auto Integer n = Maximum ( numberOfFunctionals, 0 ) ;
            /* . Basic data. */
            self->functionals         = NULL  ;
            self->hasLaplacian        = False ;
            self->hasSigma            = False ;
            self->hasTau              = False ;
            self->isSpinRestricted    = True  ;
            self->numberOfFunctionals = n     ;
            self->order               = -1    ;
            /* . Functionals array. */
            if ( n > 0 )
            {
                self->functionals = Memory_AllocateArrayOfTypes ( n, xc_func_type ) ;
                if ( self->functionals == NULL ) DFTFunctionalModel_Deallocate ( &self ) ;
            }
        }
        if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
DFTFunctionalModel *DFTFunctionalModel_Clone ( const DFTFunctionalModel *self, Status *status )
{
    DFTFunctionalModel *clone = NULL ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto IntegerArray1D *ids ;
        ids = IntegerArray1D_AllocateWithExtent ( self->numberOfFunctionals, status ) ;
        if ( ids != NULL )
        {
            auto Integer f ;
            for ( f = 0 ; f < self->numberOfFunctionals ; f++ ) Array1D_Item ( ids, f ) = self->functionals[f].info->number ;
            clone = DFTFunctionalModel_MakeFromIDs ( ids, self->isSpinRestricted, status ) ;
            IntegerArray1D_Deallocate ( &ids ) ;
        }
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DFTFunctionalModel_Deallocate ( DFTFunctionalModel **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        auto Integer f ;
        for ( f = 0 ; f < (*self)->numberOfFunctionals ; f++ ) xc_func_end ( &((*self)->functionals[f]) ) ;
        Memory_Deallocate ( (*self)->functionals ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . The exchange scaling.
! . Currently the maximum value is taken - perhaps an error should be flagged if there is more than one?
!---------------------------------------------------------------------------------------------------------------------------------*/
Real DFTFunctionalModel_ExchangeScaling ( const DFTFunctionalModel *self )
{
    Real scaling = 0.0e+00 ;
    if ( self != NULL )
    {
        auto Integer       f ;
        auto xc_func_type *functional ;
        for ( f = 0 ; f < self->numberOfFunctionals ; f++ )
        {
            functional = &(self->functionals[f]) ;
            switch ( functional->info->family )
            {
                case XC_FAMILY_HYB_GGA:
                    scaling = Maximum ( scaling, xc_hyb_gga_exx_coef  ( functional->gga  ) ) ;
                    break ;
                case XC_FAMILY_HYB_MGGA:
                    scaling = Maximum ( scaling, xc_hyb_mgga_exx_coef ( functional->mgga ) ) ;
            }
        }
    }
    return scaling ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Evaluation of the energy density and its first derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DFTFunctionalModel_Evaluate ( const DFTFunctionalModel *self, DFTIntegratorDataBlock *data )
{
    if ( ( self != NULL ) && ( data != NULL ) )
    {
        auto Integer       f ;
        auto xc_func_type *functional ;
        DFTIntegratorDataBlock_Initialize ( data ) ;
        for ( f = 0 ; f < self->numberOfFunctionals ; f++ )
        {
            functional = &(self->functionals[f]) ;
            switch ( functional->info->family )
            {
                case XC_FAMILY_LDA:
                    xc_lda_exc_vxc ( functional                           ,
                                     data->numberOfPoints                 ,
                                     Array_DataPointer ( data->rho       ) ,
                                     Array_DataPointer ( data->localEXC  ) ,
                                     Array_DataPointer ( data->localVRho ) ) ;
                    break ;
                case XC_FAMILY_GGA:
                case XC_FAMILY_HYB_GGA:
                    xc_gga_exc_vxc ( functional                             ,
                                     data->numberOfPoints                   ,
                                     Array_DataPointer ( data->rho         ) ,
                                     Array_DataPointer ( data->sigma       ) ,
                                     Array_DataPointer ( data->localEXC    ) ,
                                     Array_DataPointer ( data->localVRho   ) ,
                                     Array_DataPointer ( data->localVSigma ) ) ;
                    break ;
                case XC_FAMILY_MGGA:
                case XC_FAMILY_HYB_MGGA:
                    xc_mgga_exc_vxc ( functional                                    ,
                                      data->numberOfPoints                          ,
                                      Array_DataPointer ( data->rho                ) ,
                                      Array_DataPointer ( data->sigma              ) ,
                                      Array_DataPointer ( data->laplacianRho       ) ,
                                      Array_DataPointer ( data->tau                ) ,
                                      Array_DataPointer ( data->localEXC           ) ,
                                      Array_DataPointer ( data->localVRho          ) ,
                                      Array_DataPointer ( data->localVSigma        ) ,
                                      Array_DataPointer ( data->localVLaplacianRho ) ,
                                      Array_DataPointer ( data->localVTau          ) ) ;
            }
            DFTIntegratorDataBlock_Accumulate ( data ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Constructor given an array of functional IDs.
!---------------------------------------------------------------------------------------------------------------------------------*/
DFTFunctionalModel *DFTFunctionalModel_MakeFromIDs ( const IntegerArray1D *ids, const Boolean isSpinRestricted, Status *status )
{
    DFTFunctionalModel *self = NULL ;
    if ( ( ids != NULL ) && ( View1D_Extent ( ids ) > 0 ) && Status_IsOK ( status ) )
    {
        self = DFTFunctionalModel_Allocate ( View1D_Extent ( ids ), status ) ;
        if ( self != NULL )
        {
            auto Boolean isOK ;
            auto Integer f, failures, id, spin ;
            self->isSpinRestricted = isSpinRestricted ;
            if ( isSpinRestricted ) spin = XC_UNPOLARIZED ;
            else                    spin = XC_POLARIZED   ;
            for ( f = failures = 0 ; f < self->numberOfFunctionals ; f++ )
            {
                id   = Array1D_Item ( ids, f ) ;
                isOK = xc_func_init ( &(self->functionals[f]), id, spin ) >= 0 ;
                if ( isOK )
                {
                    switch ( self->functionals[f].info->family )
                    {
                        case XC_FAMILY_LDA:
                            self->order        = Maximum ( self->order, 0 ) ;
                            break ;
                        case XC_FAMILY_GGA:
                        case XC_FAMILY_HYB_GGA:
                            self->hasSigma     = True ;
                            self->order        = Maximum ( self->order, 1 ) ;
                            break ;
                        case XC_FAMILY_MGGA:
                        case XC_FAMILY_HYB_MGGA:
                            self->hasLaplacian = True ;
                            self->hasSigma     = True ;
                            self->hasTau       = True ;
                            self->order        = Maximum ( self->order, 2 ) ;
                    }
                }
                else failures += 1 ;
            }
            if ( failures > 0 )
            {
                DFTFunctionalModel_Deallocate ( &self ) ;
                Status_Set ( status, Status_InvalidArgument ) ;
            }
        }
        else Status_Set ( status, Status_OutOfMemory ) ;
    }
    return self ;
}
