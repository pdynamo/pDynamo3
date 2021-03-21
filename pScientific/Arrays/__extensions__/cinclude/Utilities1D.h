# ifndef _UTILITIES1D_MACROS
# define _UTILITIES1D_MACROS

# include <math.h>

# ifdef _UseCBLAS
# include "cblas.h"
# endif
# include "NumericalMacros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Absolute maximum - result needs to be initialized on entry. */
# ifdef _UseReal
# define Utilities1D_AbsoluteMaximum( extent, data, stride, result ) \
    { \
        auto Integer e, i ; \
        for ( e = i = 0 ; e < extent ; e++ , i += stride ) result = Maximum ( result, fabs ( data[i] ) ) ; \
    }
# else
# define Utilities1D_AbsoluteMaximum( extent, data, stride, result ) \
    { \
        auto Integer e, i ; \
        for ( e = i = 0 ; e < extent ; e++ , i += stride ) result = Maximum ( result,  abs ( data[i] ) ) ; \
    }
# endif

/* . Adding. */
# ifdef _UseCBLAS
# define Utilities1D_Add( extent1, data1, stride1, extent2, data2, stride2, alpha, status ) \
    { \
        if ( extent1 == extent2 ) \
        { \
            cblas_daxpy ( extent1, alpha, data2, stride2, data1, stride1 ) ; \
        } \
        else Status_Set ( status, Status_NonConformableArrays ) ; \
    }
# else
# define Utilities1D_Add( extent1, data1, stride1, extent2, data2, stride2, alpha, status ) \
    { \
        if ( extent1 == extent2 ) \
        { \
            auto Integer e, i1, i2 ; \
            for ( e = i1 = i2 = 0 ; e < extent1 ; e++, i1 += stride1, i2 += stride2 ) data1[i1] += ( alpha * data2[i2] ) ; \
        } \
        else Status_Set ( status, Status_NonConformableArrays ) ; \
    }
# endif

/* . Copying. */
# ifdef _UseCBLAS
# define Utilities1D_CopyTo( extent1, data1, stride1, extent2, data2, stride2, status ) \
    { \
        if ( extent1 == extent2 ) \
        { \
            cblas_dcopy ( extent1, data1, stride1, data2, stride2 ) ; \
        } \
        else Status_Set ( status, Status_NonConformableArrays ) ; \
    }
# else
# define Utilities1D_CopyTo( extent1, data1, stride1, extent2, data2, stride2, status ) \
    { \
        if ( extent1 == extent2 ) \
        { \
            auto Integer e, i1, i2 ; \
            for ( e = i1 = i2 = 0 ; e < extent1 ; e++, i1 += stride1, i2 += stride2 ) data2[i2] = data1[i1] ; \
        } \
        else Status_Set ( status, Status_NonConformableArrays ) ; \
    }
# endif

/* . Dot product - result needs to be initialized on entry. */
# ifdef _UseCBLAS
# define Utilities1D_Dot( extent1, data1, stride1, extent2, data2, stride2, result, status ) \
    { \
        if ( extent1 == extent2 ) \
        { \
            result += cblas_ddot ( extent1, data1, stride1, data2, stride2 ) ; \
        } \
        else Status_Set ( status, Status_NonConformableArrays ) ; \
    }
# else
# define Utilities1D_Dot( extent1, data1, stride1, extent2, data2, stride2, result, status ) \
    { \
        if ( extent1 == extent2 ) \
        { \
            auto Integer e, i1, i2 ; \
            for ( e = i1 = i2 = 0 ; e < extent1 ; e++, i1 += stride1, i2 += stride2 ) result += ( data1[i1] * data2[i2] ) ; \
        } \
        else Status_Set ( status, Status_NonConformableArrays ) ; \
    }
# endif

/* . Incrementing. */
# define Utilities1D_Increment( extent, data, stride, value ) \
    { \
        auto Integer e, i ; \
        for ( e = i = 0 ; e < extent ; e++ , i += stride ) data[i] += value ; \
    }

/* . Maximum. */
# define Utilities1D_Maximum( extent, data, stride, result ) \
    { \
        auto Integer e, i ; \
        result = data[0]  ; \
        for ( e = 1, i = stride ; e < extent ; e++ , i += stride ) result = Maximum ( result,  data[i] ) ; \
    }

/* . Minimum. */
# define Utilities1D_Minimum( extent, data, stride, result ) \
    { \
        auto Integer e, i ; \
        result = data[0]  ; \
        for ( e = 1, i = stride ; e < extent ; e++ , i += stride ) result = Minimum ( result,  data[i] ) ; \
    }

/* . Multiplying. */
# define Utilities1D_Multiply( extent1, data1, stride1, extent2, data2, stride2, status ) \
    { \
        if ( extent1 == extent2 ) \
        { \
            auto Integer e, i1, i2 ; \
            for ( e = i1 = i2 = 0 ; e < extent1 ; e++, i1 += stride1, i2 += stride2 ) data1[i1] *= data2[i2] ; \
        } \
        else Status_Set ( status, Status_NonConformableArrays ) ; \
    }

/* . Scaling. */
# ifdef _UseCBLAS
# define Utilities1D_Scale( extent, data, stride, value ) \
    { \
        cblas_dscal ( extent, value, data, stride ) ; \
    }
# else
# define Utilities1D_Scale( extent, data, stride, value ) \
    { \
        auto Integer e, i ; \
        for ( e = i = 0 ; e < extent ; e++ , i += stride ) data[i] *= value ; \
    }
# endif

/* . Setting. */
# define Utilities1D_Set( extent, data, stride, value ) \
    { \
        auto Integer e, i ; \
        for ( e = i = 0 ; e < extent ; e++ , i += stride ) data[i] = value ; \
    }

# endif

/* . Summing - result needs to be initialized on entry. */
# define Utilities1D_Sum( extent, data, stride, result ) \
    { \
        auto Integer e, i ; \
        for ( e = i = 0 ; e < extent ; e++ , i += stride ) result += data[i] ; \
    }
