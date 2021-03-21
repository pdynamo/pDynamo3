/*==================================================================================================================================
! . Symmetry parameter gradient functions.
!=================================================================================================================================*/

/* . Notes:

   Both r/H and f/H representations are catered for but care should be taken to ensure that they are used consistently.

*/

# include <math.h>

# include "Memory.h"
# include "SymmetryParameterGradients.h"
# include "Units.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
SymmetryParameterGradients *SymmetryParameterGradients_Allocate ( void )
{
    SymmetryParameterGradients *self = Memory_AllocateType ( SymmetryParameterGradients ) ;
    self->dEdH = Matrix33_Allocate ( ) ;
    SymmetryParameterGradients_Initialize ( self ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation given dEdH.
!---------------------------------------------------------------------------------------------------------------------------------*/
SymmetryParameterGradients *SymmetryParameterGradients_AllocateWithMatrix ( const Matrix33 *dEdH   ,
                                                                                  Status   *status )
{
    SymmetryParameterGradients *self = NULL ;
    if ( ( dEdH != NULL ) && Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( SymmetryParameterGradients ) ;
        if ( self != NULL )
        {
            self->dEdH = Matrix33_CloneShallow ( dEdH, status ) ;
            SymmetryParameterGradients_Initialize ( self ) ;
        }
        if ( ! Status_IsOK ( status ) ) SymmetryParameterGradients_Deallocate ( &self ) ;
        if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameterGradients_Deallocate ( SymmetryParameterGradients **self )
{
    if ( (*self) != NULL )
    {
        Matrix33_Deallocate ( &((*self)->dEdH) ) ;
        Memory_Deallocate   (   (*self)        ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Convert dEdH to dEda, dEdb, dEdc, dEdalpha, dEdbeta, dEdgamma.
! . This depends upon how H is defined in symmetryParameters.
! . The procedure works for both the r/H and f/H representations.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameterGradients_CrystalDerivatives ( SymmetryParameterGradients *self, const SymmetryParameters *symmetryParameters )
{
    if ( ( self != NULL ) && ( symmetryParameters != NULL ) )
    {
        auto Real alpha, beta, cosAlpha, cosBeta, cosGamma, fact12, fact22, gamma, sinGamma ;
        /* . Some factors. */
        alpha    = symmetryParameters->alpha * Units_Angle_Degrees_To_Radians ;
        beta     = symmetryParameters->beta  * Units_Angle_Degrees_To_Radians ;
        gamma    = symmetryParameters->gamma * Units_Angle_Degrees_To_Radians ;
        cosAlpha = cos ( alpha ) ;
        cosBeta  = cos ( beta  ) ;
        cosGamma = cos ( gamma ) ;
        sinGamma = sin ( gamma ) ;
        fact12   = ( cosAlpha - cosBeta * cosGamma ) ;
        fact22   = sqrt ( 1.0e+00 - cosAlpha * cosAlpha - cosBeta * cosBeta - cosGamma * cosGamma + 2.0e+00 * cosAlpha * cosBeta * cosGamma ) ;
        /* . The derivatives - a, b, c. */
        self->dEda =            Matrix33_Item ( self->dEdH, 0, 0 ) ;
        self->dEdb = cosGamma * Matrix33_Item ( self->dEdH, 0, 1 ) +
                     sinGamma * Matrix33_Item ( self->dEdH, 1, 1 ) ;
        self->dEdc = cosBeta  * Matrix33_Item ( self->dEdH, 0, 2 ) +
                     fact12   * Matrix33_Item ( self->dEdH, 1, 2 ) / sinGamma +
                     fact22   * Matrix33_Item ( self->dEdH, 2, 2 ) / sinGamma ;
        /* . The derivatives - alpha, beta, gamma. */
        self->dEdalpha = - (          Matrix33_Item ( self->dEdH, 1, 2 ) -
                             fact12 * Matrix33_Item ( self->dEdH, 2, 2 ) / fact22 ) *
                           ( symmetryParameters->c * sin ( alpha ) * Units_Angle_Degrees_To_Radians ) / sinGamma ;
        self->dEdbeta  = - (   sinGamma                        * Matrix33_Item ( self->dEdH, 0, 2 ) -
                               cosGamma                        * Matrix33_Item ( self->dEdH, 1, 2 ) +
                             ( cosAlpha * cosGamma - cosBeta ) * Matrix33_Item ( self->dEdH, 2, 2 ) / fact22 ) *
                           ( symmetryParameters->c * sin ( beta ) * Units_Angle_Degrees_To_Radians ) / sinGamma ;

        self->dEdgamma =   ( symmetryParameters->b *
                                ( - sinGamma * Matrix33_Item ( self->dEdH, 0, 1 ) +
                                    cosGamma * Matrix33_Item ( self->dEdH, 1, 1 ) ) +
                             symmetryParameters->c *
                                ( ( cosBeta - fact12 * cosGamma / ( sinGamma * sinGamma ) ) * Matrix33_Item ( self->dEdH, 1, 2 ) +
                                  ( ( cosGamma - cosAlpha * cosBeta ) / fact22 - fact22 * cosGamma / ( sinGamma * sinGamma ) ) *
                                               Matrix33_Item ( self->dEdH, 2, 2 ) ) ) * Units_Angle_Degrees_To_Radians ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Convert r/H to f/H derivatives (in-place).
! . Will need selection.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameterGradients_FractionalDerivatives (       SymmetryParameterGradients *self               ,
                                                        const SymmetryParameters         *symmetryParameters ,
                                                        const Coordinates3               *coordinates3       ,
                                                              Coordinates3               *gradients3         )
{
    if ( ( self               != NULL ) &&
         ( symmetryParameters != NULL ) &&
         ( coordinates3       != NULL ) &&
         ( gradients3         != NULL ) &&
         ( Coordinates3_Rows ( coordinates3 ) == Coordinates3_Rows ( gradients3 ) ) )
    {
        auto Real fx, fy, fz, gx, gy, gz, im00, im01, im02, im10, im11, im12, im20, im21, im22, m00, m01, m02, m10, m11, m12, m20, m21, m22, x, y, z ;
        auto Integer    i ;
        /* . Get local variables for H and inverseH. */
        m00  = Matrix33_Item ( symmetryParameters->H, 0, 0 ) ;
        m01  = Matrix33_Item ( symmetryParameters->H, 0, 1 ) ;
        m02  = Matrix33_Item ( symmetryParameters->H, 0, 2 ) ;
        m10  = Matrix33_Item ( symmetryParameters->H, 1, 0 ) ;
        m11  = Matrix33_Item ( symmetryParameters->H, 1, 1 ) ;
        m12  = Matrix33_Item ( symmetryParameters->H, 1, 2 ) ;
        m20  = Matrix33_Item ( symmetryParameters->H, 2, 0 ) ;
        m21  = Matrix33_Item ( symmetryParameters->H, 2, 1 ) ;
        m22  = Matrix33_Item ( symmetryParameters->H, 2, 2 ) ;
        im00 = Matrix33_Item ( symmetryParameters->inverseH, 0, 0 ) ;
        im01 = Matrix33_Item ( symmetryParameters->inverseH, 0, 1 ) ;
        im02 = Matrix33_Item ( symmetryParameters->inverseH, 0, 2 ) ;
        im10 = Matrix33_Item ( symmetryParameters->inverseH, 1, 0 ) ;
        im11 = Matrix33_Item ( symmetryParameters->inverseH, 1, 1 ) ;
        im12 = Matrix33_Item ( symmetryParameters->inverseH, 1, 2 ) ;
        im20 = Matrix33_Item ( symmetryParameters->inverseH, 2, 0 ) ;
        im21 = Matrix33_Item ( symmetryParameters->inverseH, 2, 1 ) ;
        im22 = Matrix33_Item ( symmetryParameters->inverseH, 2, 2 ) ;
        /* . Loop over the coordinates/gradients. */
        for ( i = 0 ; i < Coordinates3_Rows ( coordinates3 ) ; i++ )
        {
            Coordinates3_GetRow ( coordinates3, i,  x,  y,  z ) ;
            Coordinates3_GetRow ( gradients3,   i, gx, gy, gz ) ;
            /* . Get f. */
            fx = im00 * x + im01 * y + im02 * z ;
            fy = im10 * x + im11 * y + im12 * z ;
            fz = im20 * x + im21 * y + im22 * z ;
            /* . dEdH. */
            Matrix33_Item ( self->dEdH, 0, 0 ) += fx * gx ;
            Matrix33_Item ( self->dEdH, 0, 1 ) += fy * gx ;
            Matrix33_Item ( self->dEdH, 0, 2 ) += fz * gx ;
            Matrix33_Item ( self->dEdH, 1, 0 ) += fx * gy ;
            Matrix33_Item ( self->dEdH, 1, 1 ) += fy * gy ;
            Matrix33_Item ( self->dEdH, 1, 2 ) += fz * gy ;
            Matrix33_Item ( self->dEdH, 2, 0 ) += fx * gz ;
            Matrix33_Item ( self->dEdH, 2, 1 ) += fy * gz ;
            Matrix33_Item ( self->dEdH, 2, 2 ) += fz * gz ;
            /* . dEdf. */
            Coordinates3_Item ( gradients3, i, 0 ) = m00 * gx + m10 * gy + m20 * gz ;
            Coordinates3_Item ( gradients3, i, 1 ) = m01 * gx + m11 * gy + m21 * gz ;
            Coordinates3_Item ( gradients3, i, 2 ) = m02 * gx + m12 * gy + m22 * gz ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the derivatives due to image terms (r/H formalism).
! . |transformation3| is the fractional transformation (without H).
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameterGradients_ImageDerivatives (       SymmetryParameterGradients *self               ,
                                                   const SymmetryParameters         *symmetryParameters ,
                                                   const Transformation3            *transformation3    ,
                                                   const Coordinates3               *coordinates3       ,
                                                   const Coordinates3               *gradients3         )
{
    if ( ( self               != NULL ) &&
         ( symmetryParameters != NULL ) &&
         ( transformation3    != NULL ) &&
         ( coordinates3       != NULL ) &&
         ( gradients3         != NULL ) &&
         ( Coordinates3_Rows ( coordinates3 ) == Coordinates3_Rows ( gradients3 ) ) )
    {
        auto Integer   a, b, i ;
        auto Real      dx, dy, dz, gx, gy, gz, r00, r01, r02, r10, r11, r12, r20, r21, r22, sum, t, x, y, z ;
        auto Matrix33 *di, *ms, *si ;
        /* . Allocate space. */
        di = Matrix33_Allocate ( ) ;
        ms = Matrix33_Allocate ( ) ;
        si = Matrix33_Allocate ( ) ;
        /* . Intermediate quantities (MS and SM^-1). */
        Matrix33_CopyTo ( transformation3->rotation, ms, NULL ) ; Matrix33_PreMultiplyBy  ( ms, symmetryParameters->H        ) ;
        Matrix33_CopyTo ( transformation3->rotation, si, NULL ) ; Matrix33_PostMultiplyBy ( si, symmetryParameters->inverseH ) ;
        /* . Loop over the elements of M. */
        for ( a = 0 ; a < 3 ; a++ )
        {
            for ( b = 0 ; b < 3 ; b++ )
            {
                /* . Rotational contribution. */
                /* . Inverse derivative. */
                Matrix33_InverseDerivative ( symmetryParameters->H, a, b, di ) ;
                Matrix33_PreMultiplyBy ( di, ms ) ;
                /* . Non-inverse derivative. */
                Matrix33_GetRow       ( si, b, dx, dy, dz ) ;
                Matrix33_IncrementRow ( di, a, dx, dy, dz ) ;
                /* . Local variables. */
                r00  = Matrix33_Item ( di, 0, 0 ) ;
                r01  = Matrix33_Item ( di, 0, 1 ) ;
                r02  = Matrix33_Item ( di, 0, 2 ) ;
                r10  = Matrix33_Item ( di, 1, 0 ) ;
                r11  = Matrix33_Item ( di, 1, 1 ) ;
                r12  = Matrix33_Item ( di, 1, 2 ) ;
                r20  = Matrix33_Item ( di, 2, 0 ) ;
                r21  = Matrix33_Item ( di, 2, 1 ) ;
                r22  = Matrix33_Item ( di, 2, 2 ) ;
                /* . Translational contribution. */
                t = Vector3_Item ( transformation3->translation, b ) ;
                /* . Loop over the coordinates/gradients. */
                for ( i = 0, sum = 0.0e+00 ; i < Coordinates3_Rows ( coordinates3 ) ; i++ )
                {
                    Coordinates3_GetRow ( coordinates3, i,  x,  y,  z ) ;
                    Coordinates3_GetRow ( gradients3,   i, gx, gy, gz ) ;
                    /* . Rotational contribution. */
                    dx = r00 * x + r01 * y + r02 * z ;
                    dy = r10 * x + r11 * y + r12 * z ;
                    dz = r20 * x + r21 * y + r22 * z ;
                    /* . Translational contribution. */
                    switch ( a )
                    {
                        case 0: dx += t ; break ;
                        case 1: dy += t ; break ;
                        case 2: dz += t ; break ;
                    }
                    /* . dEdH. */
                    sum += dx * gx + dy * gy + dz * gz ;
                }
                /* . Set the component of dEdH. */
                Matrix33_Item ( self->dEdH, a, b ) += sum ;
            }
        }
        /* . Deallocate space. */
        Matrix33_Deallocate ( &di ) ;
        Matrix33_Deallocate ( &ms ) ;
        Matrix33_Deallocate ( &si ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameterGradients_Initialize ( SymmetryParameterGradients *self )
{
    if ( self != NULL )
    {
        self->dEda     = 0.0e+00 ;
        self->dEdb     = 0.0e+00 ;
        self->dEdc     = 0.0e+00 ;
        self->dEdalpha = 0.0e+00 ;
        self->dEdbeta  = 0.0e+00 ;
        self->dEdgamma = 0.0e+00 ;
        Matrix33_Set ( self->dEdH, 0.0e+00 ) ;
    }
}

