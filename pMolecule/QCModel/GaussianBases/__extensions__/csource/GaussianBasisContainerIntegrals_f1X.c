/*==================================================================================================================================
! . Container integrals - 1 basis, 1 electron.
!=================================================================================================================================*/

# include <math.h>
# include <stdlib.h>

# include "GaussianBasisContainerIntegrals_f1X.h"
# include "GaussianBasisIntegrals_f1X.h"
# include "Integer.h"
# include "RealUtilities.h"

/* . All output integral arrays are overwritten by these functions. */

/* . Will need to modify output integrals if transforming from a Cartesian to a spherical representation
     by having an intermediate Cartesian integral array. */

/*----------------------------------------------------------------------------------------------------------------------------------
! . Dipole integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainerIntegrals_f1Di ( const GaussianBasisContainer *self         ,
                                            const Coordinates3           *coordinates3 ,
                                                  Vector3                *center       ,
                                                  RealArray1D            *dipoleX      ,
                                                  RealArray1D            *dipoleY      ,
                                                  RealArray1D            *dipoleZ      ,
                                                  Status                 *status       )
{
    RealArray1D_Set ( dipoleX, 0.0e+00 ) ;
    RealArray1D_Set ( dipoleY, 0.0e+00 ) ;
    RealArray1D_Set ( dipoleZ, 0.0e+00 ) ;
    if ( ( self         != NULL ) &&
         ( coordinates3 != NULL ) &&
         ( dipoleX      != NULL ) &&
         ( dipoleY      != NULL ) &&
         ( dipoleZ      != NULL ) &&
         Status_IsOK ( status ) )
    {
        if ( ( View1D_Extent ( dipoleX ) == Array1D_Item ( self->centerFunctionPointers, self->capacity ) ) ||
             ( View1D_Extent ( dipoleY ) == Array1D_Item ( self->centerFunctionPointers, self->capacity ) ) ||
             ( View1D_Extent ( dipoleZ ) == Array1D_Item ( self->centerFunctionPointers, self->capacity ) ) )
        {
            auto Integer  s1 ;
            auto Real    *rWork  ;
            auto Vector3 *origin ;
            s1    = GaussianBasisContainer_LargestShell ( self, True ) ;
            rWork = Real_Allocate ( 4*s1, status ) ;
            if ( center == NULL ) { origin = Vector3_Allocate ( ) ; Vector3_Set ( origin, 0.0e+00 ) ; }
            else                    origin = center ;
            if ( ( origin != NULL ) && Status_IsOK ( status ) )
            {
                auto Integer        i, start, stop ;
                auto GaussianBasis *iBasis ;
                auto RealArray1D    viewX, viewY, viewZ ;
                for ( i = 0 ; i < self->capacity ; i++ )
                {
                    iBasis = self->entries[i] ;
                    start  = Array1D_Item ( self->centerFunctionPointers, i   ) ;
                    stop   = Array1D_Item ( self->centerFunctionPointers, i+1 ) ;
                    RealArray1D_View ( dipoleX, start, ( stop - start ), 1, False, &viewX, NULL ) ;
                    RealArray1D_View ( dipoleY, start, ( stop - start ), 1, False, &viewY, NULL ) ;
                    RealArray1D_View ( dipoleZ, start, ( stop - start ), 1, False, &viewZ, NULL ) ;
                    GaussianBasisIntegrals_f1Di ( iBasis                                      ,
                                                  Coordinates3_RowPointer ( coordinates3, i ) ,
                                                  Vector3_Data ( origin )                     ,
                                                  s1                                          ,
                                                  rWork                                       ,
                                                  &viewX                                      ,
                                                  &viewY                                      ,
                                                  &viewZ                                      ) ;
                }
            }
            if ( center == NULL ) Vector3_Deallocate ( &origin ) ;
            Real_Deallocate ( &rWork ) ;
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Overlap integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainerIntegrals_f1Oi ( const GaussianBasisContainer *self         ,
                                                  RealArray1D            *overlap      ,
                                                  Status                 *status       )
{
    RealArray1D_Set ( overlap, 0.0e+00 ) ;
    if ( ( self         != NULL ) &&
         ( overlap      != NULL ) &&
         Status_IsOK ( status ) )
    {
        if ( View1D_Extent ( overlap ) == Array1D_Item ( self->centerFunctionPointers, self->capacity ) )
        {
            auto Integer  s1 ;
            auto Real    *rWork  ;
            s1    = GaussianBasisContainer_LargestShell ( self, True ) ;
            rWork = Real_Allocate ( 2*s1, status ) ;
            if ( Status_IsOK ( status ) )
            {
                auto Integer     i, start, stop ;
                auto RealArray1D view ;
                for ( i = 0 ; i < self->capacity ; i++ )
                {
                    start = Array1D_Item ( self->centerFunctionPointers, i   ) ;
                    stop  = Array1D_Item ( self->centerFunctionPointers, i+1 ) ;
                    RealArray1D_View ( overlap, start, ( stop - start ), 1, False, &view, NULL ) ;
                    GaussianBasisIntegrals_f1Oi ( self->entries[i], s1, rWork, &view ) ;
                }
            }
            Real_Deallocate ( &rWork ) ;
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Quadrupole integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainerIntegrals_f1Qi ( const GaussianBasisContainer *self         ,
                                            const Coordinates3           *coordinates3 ,
                                                  Vector3                *center       ,
                                                  RealArray1D            *qXX          ,
                                                  RealArray1D            *qYY          ,
                                                  RealArray1D            *qZZ          ,
                                                  RealArray1D            *qXY          ,
                                                  RealArray1D            *qXZ          ,
                                                  RealArray1D            *qYZ          ,
                                                  Status                 *status       )
{
    RealArray1D_Set ( qXX, 0.0e+00 ) ;
    RealArray1D_Set ( qYY, 0.0e+00 ) ;
    RealArray1D_Set ( qZZ, 0.0e+00 ) ;
    RealArray1D_Set ( qXY, 0.0e+00 ) ;
    RealArray1D_Set ( qXZ, 0.0e+00 ) ;
    RealArray1D_Set ( qYZ, 0.0e+00 ) ;
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
        if ( ( View1D_Extent ( qXX ) == Array1D_Item ( self->centerFunctionPointers, self->capacity ) ) ||
             ( View1D_Extent ( qYY ) == Array1D_Item ( self->centerFunctionPointers, self->capacity ) ) ||
             ( View1D_Extent ( qZZ ) == Array1D_Item ( self->centerFunctionPointers, self->capacity ) ) ||
             ( View1D_Extent ( qXY ) == Array1D_Item ( self->centerFunctionPointers, self->capacity ) ) ||
             ( View1D_Extent ( qXZ ) == Array1D_Item ( self->centerFunctionPointers, self->capacity ) ) ||
             ( View1D_Extent ( qYZ ) == Array1D_Item ( self->centerFunctionPointers, self->capacity ) ) )
        {
            auto Integer  s1 ;
            auto Real    *rWork  ;
            auto Vector3 *origin ;
            s1    = GaussianBasisContainer_LargestShell ( self, True ) ;
            rWork = Real_Allocate ( 7*s1, status ) ;
            if ( center == NULL ) { origin = Vector3_Allocate ( ) ; Vector3_Set ( origin, 0.0e+00 ) ; }
            else                    origin = center ;
            if ( ( origin != NULL ) && Status_IsOK ( status ) )
            {
                auto Integer        i, start, stop ;
                auto GaussianBasis *iBasis ;
                auto RealArray1D    viewXX, viewXY, viewXZ, viewYY, viewYZ, viewZZ ;
                for ( i = 0 ; i < self->capacity ; i++ )
                {
                    iBasis = self->entries[i] ;
                    start  = Array1D_Item ( self->centerFunctionPointers, i   ) ;
                    stop   = Array1D_Item ( self->centerFunctionPointers, i+1 ) ;
                    RealArray1D_View ( qXX, start, ( stop - start ), 1, False, &viewXX, NULL ) ;
                    RealArray1D_View ( qYY, start, ( stop - start ), 1, False, &viewYY, NULL ) ;
                    RealArray1D_View ( qZZ, start, ( stop - start ), 1, False, &viewZZ, NULL ) ;
                    RealArray1D_View ( qXY, start, ( stop - start ), 1, False, &viewXY, NULL ) ;
                    RealArray1D_View ( qXZ, start, ( stop - start ), 1, False, &viewXZ, NULL ) ;
                    RealArray1D_View ( qYZ, start, ( stop - start ), 1, False, &viewYZ, NULL ) ;
                    GaussianBasisIntegrals_f1Qi ( iBasis                                      ,
                                                  Coordinates3_RowPointer ( coordinates3, i ) ,
                                                  Vector3_Data ( origin )                     ,
                                                  s1                                          ,
                                                  rWork                                       ,
                                                  &viewXX                                     ,
                                                  &viewYY                                     ,
                                                  &viewZZ                                     ,
                                                  &viewXY                                     ,
                                                  &viewXZ                                     ,
                                                  &viewYZ                                     ) ;
                }
            }
            if ( center == NULL ) Vector3_Deallocate ( &origin ) ;
            Real_Deallocate ( &rWork ) ;
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}
