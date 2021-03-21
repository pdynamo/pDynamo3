/*==================================================================================================================================
! . 1-D real arrays.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "NumericalMacros.h"
# include "RealArray1D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Basic definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _ArrayDataFormat          "%20.10f"
# define _ArrayDataPerLine         6
# define _ArrayDataType            Real
# define _ArrayDataTypeInitializer 0.0e+00
# define _UseCBLAS
# define _UseReal
# include "Array1D_Body.i"
# undef _ArrayDataFormat
# undef _ArrayDataPerLine
# undef _ArrayDataType
# undef _ArrayDataTypeInitializer
# undef _UseCBLAS
# undef _UseReal

/*----------------------------------------------------------------------------------------------------------------------------------
! . Specific functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _UseCBLAS
# define _UseReal
# include "Utilities1D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Norm-2.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real RealArray1D_Norm2 ( const RealArray1D  *self )
{
    Real value = 0.0e+00 ;
    if ( self != NULL )
    {
        Utilities1D_Dot ( self->extent, self->data, self->stride, self->extent, self->data, self->stride, value, NULL ) ;
        value = sqrt ( value ) ;
    }
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Normalization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RealArray1D_Normalize ( RealArray1D *self, const Real *nullNormValue, Status *status )
{
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto Real delta, norm2 ;
        /* . Get the 2-norm. */
        norm2 = RealArray1D_Norm2 ( self ) ;
        /* . Get the null-norm value for normalization. */
        if ( nullNormValue == NULL ) delta = Real_SafeMinimum ;
        else                         delta = Maximum ( Real_SafeMinimum, fabs ( (*nullNormValue) ) ) ;
        /* . Normalize or set to zero. */
        if ( norm2 > delta )
        {
            RealArray1D_Scale ( self, 1.0e+00 / norm2 ) ;
        }
        else
        {
            RealArray1D_Set ( self, 0.0e+00 ) ;
            Status_Set ( status, Status_AlgorithmError ) ;
        }
    }
}

# undef _UseCBLAS
# undef _UseReal
