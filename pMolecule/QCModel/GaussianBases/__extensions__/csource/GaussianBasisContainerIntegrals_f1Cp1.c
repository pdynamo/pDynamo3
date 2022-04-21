/*==================================================================================================================================
! . Container integrals - 1 basis, 1 electron, 1 nucleus/point.
!=================================================================================================================================*/

# include <math.h>
# include <stdlib.h>

# include "GaussianBasisContainerIntegrals_f1Cp1.h"
# include "GaussianBasisIntegrals_f1Cp1.h"
# include "Integer.h"
# include "IntegerUtilities.h"
# include "Real.h"
# include "RealUtilities.h"

/* . All integral methods increment their results so initialization is required first! */

/*----------------------------------------------------------------------------------------------------------------------------------
! . Fit-nuclear/point derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainerIntegrals_f1Cm1R1 ( const GaussianBasisContainer *self            , 
                                               const RealArray1D            *charges         , 
                                               const RealArray1D            *widthsE         , 
                                               const RealArray1D            *widthsN         , 
                                               const Coordinates3           *coordinates3    , 
                                               const Coordinates3           *coordinates3G   , 
                                                     Selection              *selectionG      , 
                                               const RealArray1D            *fitCoefficients ,
                                                     Coordinates3           *gradients3      ,
                                                     Coordinates3           *gradients3G     ,
                                                     Status                 *status          )
{
    if ( ( self            != NULL ) &&
         ( charges         != NULL ) &&
         ( coordinates3    != NULL ) &&
         ( coordinates3G   != NULL ) &&
         ( fitCoefficients != NULL ) &&
         ( gradients3      != NULL ) &&
         ( gradients3G     != NULL ) )
    {
        auto Integer *iWork, s1 ;
        auto Real    *rWork ;
        Selection_MakeFlags ( selectionG, Coordinates3_Rows ( coordinates3G ), status ) ;
        s1    = GaussianBasisContainer_LargestShell ( self, True ) ;
        iWork = Integer_Allocate ( 3*s1, status ) ;
        rWork = Real_Allocate    ( 5*s1, status ) ;
        if ( Status_IsOK ( status ) )
        {
            auto Integer     i, start, stop ;
            auto Real        dRi[3] ;
            auto RealArray1D view ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                start = Array1D_Item ( self->centerFunctionPointers, i   ) ;
                stop  = Array1D_Item ( self->centerFunctionPointers, i+1 ) ;
                RealArray1D_View ( fitCoefficients, start, ( stop - start ), 1, False, &view, NULL ) ;
                GaussianBasisIntegrals_f1Cm1R1 ( self->entries[i]                            ,
                                                 Coordinates3_RowPointer ( coordinates3, i ) ,
                                                 charges                                     ,
                                                 widthsE                                     ,
                                                 widthsN                                     ,
                                                 coordinates3G                               ,
                                                 selectionG                                  ,
                                                 &view                                       ,
                                                 s1                                          ,
                                                 iWork                                       ,
                                                 rWork                                       ,
                                                 dRi                                         ,
                                                 gradients3G                                 ) ;
                Coordinates3_IncrementRow ( gradients3, i, dRi[0], dRi[1], dRi[2] ) ;
            }
        }
        Integer_Deallocate ( &iWork ) ;
        Real_Deallocate    ( &rWork ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Fit-nuclear/point integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainerIntegrals_f1Cm1V ( const GaussianBasisContainer *self          ,
                                              const RealArray1D            *charges       ,
                                              const RealArray1D            *widthsE       ,
                                              const RealArray1D            *widthsN       ,
                                              const Coordinates3           *coordinates3  ,
                                              const Coordinates3           *coordinates3G ,
                                                    Selection              *selectionG    ,
                                                    RealArray1D            *integrals     ,
                                                    Status                 *status        )
{
    if ( ( self          != NULL ) &&
         ( charges       != NULL ) &&
         ( coordinates3  != NULL ) &&
         ( coordinates3G != NULL ) &&
         ( integrals     != NULL ) )
    {
        auto Integer *iWork, s1 ;
        auto Real    *rWork ;
        Selection_MakeFlags ( selectionG, Coordinates3_Rows ( coordinates3G ), status ) ;
        s1    = GaussianBasisContainer_LargestShell ( self, True ) ;
        iWork = Integer_Allocate ( 3*s1, status ) ;
        rWork = Real_Allocate    ( 3*s1, status ) ;
        if ( Status_IsOK ( status ) )
        {
            auto Integer     i, start, stop ;
            auto RealArray1D view ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                start = Array1D_Item ( self->centerFunctionPointers, i   ) ;
                stop  = Array1D_Item ( self->centerFunctionPointers, i+1 ) ;
                RealArray1D_View ( integrals, start, ( stop - start ), 1, False, &view, NULL ) ;
                GaussianBasisIntegrals_f1Cm1V ( self->entries[i]                            ,
                                                Coordinates3_RowPointer ( coordinates3, i ) ,
                                                charges                                     ,
                                                widthsE                                     ,
                                                widthsN                                     ,
                                                coordinates3G                               ,
                                                selectionG                                  ,
                                                s1                                          ,
                                                iWork                                       ,
                                                rWork                                       ,
                                                &view                                       ) ;
            }
        }
        Integer_Deallocate ( &iWork ) ;
        Real_Deallocate    ( &rWork ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Fit-nuclear/point potentials.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainerIntegrals_f1Cp1V ( const GaussianBasisContainer *self            ,
                                              const RealArray1D            *widthsE         ,
                                              const RealArray1D            *widthsN         ,
                                              const Coordinates3           *coordinates3    ,
                                              const Coordinates3           *coordinates3G   ,
                                                    Selection              *selectionG      ,
                                              const RealArray1D            *fitCoefficients ,
                                                    RealArray1D            *potentials      ,
                                                    Status                 *status          )
{
    if ( ( self            != NULL ) &&
         ( coordinates3    != NULL ) &&
         ( coordinates3G   != NULL ) &&
         ( fitCoefficients != NULL ) &&
         ( potentials      != NULL ) )
    {
        auto Integer *iWork, s1 ;
        auto Real    *rWork ;
        Selection_MakeFlags ( selectionG, Coordinates3_Rows ( coordinates3G ), status ) ;
        s1    = GaussianBasisContainer_LargestShell ( self, True ) ;
        iWork = Integer_Allocate ( 3*s1, status ) ;
        rWork = Real_Allocate    ( 3*s1, status ) ;
        if ( Status_IsOK ( status ) )
        {
            auto Integer     i, start, stop ;
            auto RealArray1D view ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                start = Array1D_Item ( self->centerFunctionPointers, i   ) ;
                stop  = Array1D_Item ( self->centerFunctionPointers, i+1 ) ;
                RealArray1D_View ( fitCoefficients, start, ( stop - start ), 1, False, &view, NULL ) ;
                GaussianBasisIntegrals_f1Cp1V ( self->entries[i]                            ,
                                                Coordinates3_RowPointer ( coordinates3, i ) ,
                                                widthsE                                     ,
                                                widthsN                                     ,
                                                coordinates3G                               ,
                                                selectionG                                  ,
                                                &view                                       ,
                                                s1                                          ,
                                                iWork                                       ,
                                                rWork                                       ,
                                                potentials                                  ) ; /* . No view as indexed by G. */
            }
        }
        Integer_Deallocate ( &iWork ) ;
        Real_Deallocate    ( &rWork ) ;
    }
}
