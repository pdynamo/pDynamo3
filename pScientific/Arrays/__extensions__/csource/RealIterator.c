/*==================================================================================================================================
! . Real iterator operations.
!=================================================================================================================================*/

# include <errno.h>
# include <math.h>
# include <stdlib.h>

# include "cblas.h"
# include "NumericalMacros.h"
# include "RealIterator.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Unary.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { selfData[i] = fabs ( selfData[i] ) ; }
# define _SelfIsModified
void RealIterator_Absolute ( Iterator *self, Real *selfData )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation                  { result = Maximum ( result, fabs ( selfData[i] ) ) ; }
# define _Reduction
# define _ReductionResultInitializer 0.0e+00
# define _ReductionResultType        Real
Real RealIterator_AbsoluteMaximum ( Iterator *self, const Real *selfData )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _Reduction
# undef _ReductionResultInitializer
# undef _ReductionResultType
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation                  { if ( fabs ( selfData[i] ) <= tolerance ) result += 1 ; }
# define _Reduction
# define _ReductionResultInitializer 0
# define _ReductionResultType        Integer
Integer RealIterator_CountSmall ( Iterator *self, const Real *selfData, const Real tolerance )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _Reduction
# undef _ReductionResultInitializer
# undef _ReductionResultType
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _InnerLoopOperation         result += cblas_ddot ( extent, &selfData[first], stride, &selfData[first], stride ) ;
# define _Operation                  { result += ( selfData[i] * selfData[i] ) ; }
# define _Reduction
# define _ReductionResultInitializer 0.0
# define _ReductionResultType        Real
Real RealIterator_DotSelf ( Iterator *self, const Real *selfData )
{
# include "Iterator_UnaryOperation.i"
}
# undef _InnerLoopOperation
# undef _Operation
# undef _Reduction
# undef _ReductionResultInitializer
# undef _ReductionResultType
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Locals    auto Real t = log ( Real_Largest ) ;
# define _Operation { if ( selfData[i] >= t ) { selfData[i] = Real_Largest ; } else { selfData[i] = exp ( selfData[i] ) ; } }
# define _SelfIsModified
void RealIterator_Exponential ( Iterator *self, Real *selfData, Status *status )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Locals
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { if ( selfData[i] >= tolerance ) selfData[i] = value ; }
# define _SelfIsModified
extern void RealIterator_FilterGreaterThan ( Iterator *self, Real *selfData, const Real tolerance, const Real value, Status *status )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { if ( selfData[i] <= tolerance ) selfData[i] = value ; }
# define _SelfIsModified
extern void RealIterator_FilterLessThan ( Iterator *self, Real *selfData, const Real tolerance, const Real value, Status *status )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { if ( fabs ( selfData[i] ) <= tolerance ) selfData[i] = value ; }
# define _SelfIsModified
extern void RealIterator_FilterSmall ( Iterator *self, Real *selfData, const Real tolerance, const Real value, Status *status )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { selfData[i] += value ; }
# define _SelfIsModified
void RealIterator_Increment ( Iterator *self, Real *selfData, const Real value, Status *status )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation                  { result = Maximum ( result, selfData[i] ) ; }
# define _Reduction
# define _ReductionResultInitializer Real_Smallest
# define _ReductionResultType        Real
Real RealIterator_Maximum ( Iterator *self, const Real *selfData )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _Reduction
# undef _ReductionResultInitializer
# undef _ReductionResultType
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation                  { result = Minimum ( result, selfData[i] ) ; }
# define _Reduction
# define _ReductionResultInitializer Real_Largest
# define _ReductionResultType        Real
Real RealIterator_Minimum ( Iterator *self, const Real *selfData )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _Reduction
# undef _ReductionResultInitializer
# undef _ReductionResultType
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { if ( selfData[i] <= 0.0e+00 ) { selfData[i] = Real_Smallest ; } else { selfData[i] = log ( selfData[i] ) ; } }
# define _SelfIsModified
void RealIterator_NaturalLogarithm ( Iterator *self, Real *selfData, Status *status )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
Real RealIterator_Norm2 ( Iterator *self, const Real *selfData )
{
    return sqrt ( RealIterator_DotSelf ( self, selfData ) ) ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
Boolean RealIterator_Normalize ( Iterator *self, Real *selfData, const Real tolerance, Status *status )
{
    Boolean isOK  = False ;
    Real    norm2 = RealIterator_Norm2 ( self, selfData ) ;
    if ( norm2 >= tolerance ) { RealIterator_Scale ( self, selfData, 1.0 / norm2, status ) ; isOK = True ; }
    return isOK ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { selfData[i] = pow ( selfData[i], power ) ; }
# define _SelfIsModified
void RealIterator_Power ( Iterator *self, Real *selfData, const Real power, Status *status )
{
# include "Iterator_UnaryOperationWithErrorCheck.i"
}
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation                  { result *= selfData[i] ; }
# define _Reduction
# define _ReductionResultInitializer 1.0e+00
# define _ReductionResultType        Real
Real RealIterator_Product ( Iterator *self, const Real *selfData )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _Reduction
# undef _ReductionResultInitializer
# undef _ReductionResultType
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { if ( fabs ( selfData[i] ) >= tolerance ) { selfData[i] = ( 1.0 / selfData[i] ) ; } else { selfData[i] = value ; } }
# define _SelfIsModified
void RealIterator_Reciprocate ( Iterator *self, Real *selfData, const Real tolerance, const Real value, Status *status )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { if ( fabs ( selfData[i] ) >= tolerance ) { selfData[i] = ( 1.0 / pow ( selfData[i], power ) ) ; } else { selfData[i] = value ; } }
# define _SelfIsModified
void RealIterator_ReciprocatePower ( Iterator *self, Real *selfData, const Real power, const Real tolerance, const Real value, Status *status )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
Real RealIterator_RootMeanSquare ( Iterator *self, const Real *selfData )
{
    Real result = 0.0 ;
    if ( self != NULL ) result = sqrt ( RealIterator_DotSelf ( self, selfData ) / ( Real ) self->size ) ;
    return result ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _InnerLoopOperation cblas_dscal ( extent, value, &selfData[first], stride ) ;
# define _Operation          { selfData[i] *= value ; }
# define _SelfIsModified
void RealIterator_Scale ( Iterator *self, Real *selfData, const Real value, Status *status )
{
# include "Iterator_UnaryOperation.i"
}
# undef _InnerLoopOperation
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { selfData[i] = value ; }
# define _SelfIsModified
void RealIterator_Set ( Iterator *self, Real *selfData, const Real value, Status *status )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
Real RealIterator_Sparsity ( Iterator *self, const Real *selfData, const Real tolerance )
{
    Real result = 0.0 ;
    if ( self != NULL ) { result = ( 100.0 * ( Real ) RealIterator_CountSmall ( self, selfData, tolerance ) / ( Real ) self->size ) ; }
    return result ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { selfData[i] = pow ( selfData[i], 2 ) ; }
# define _SelfIsModified
void RealIterator_Square ( Iterator *self, Real *selfData, Status *status )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { selfData[i] = sqrt ( selfData[i] ) ; }
# define _SelfIsModified
void RealIterator_SquareRoot ( Iterator *self, Real *selfData, Status *status )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation                  { result += selfData[i] ; }
# define _Reduction
# define _ReductionResultInitializer 0.0e+00
# define _ReductionResultType        Real
Real RealIterator_Sum ( Iterator *self, const Real *selfData )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _Reduction
# undef _ReductionResultInitializer
# undef _ReductionResultType
/*----------------------------------------------------------------------------------------------------------------------------------
! . Binary.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _InnerLoopOperation cblas_daxpy ( extent, scale, &otherData[first2], stride2, &selfData[first1], stride1 ) ;
# define _Operation          { selfData[i1] += ( scale * otherData[i2] ) ; }
# define _SelfIsModified
void RealIterator_Add ( Iterator *self, Real *selfData, Iterator *other, const Real *otherData, const Real scale, Status *status )
{
# include "Iterator_BinaryOperation.i"
}
# undef _InnerLoopOperation
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _InnerLoopOperation cblas_dcopy ( extent, &selfData[first1], stride1, &otherData[first2], stride2 ) ;
# define _Operation          { otherData[i2] = selfData[i1] ; }
# define _OtherIsModified
void RealIterator_CopyTo ( Iterator *self, const Real *selfData, Iterator *other, Real *otherData, Status *status )
{
# include "Iterator_BinaryOperation.i"
}
# undef _InnerLoopOperation
# undef _Operation
# undef _OtherIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { if ( fabs ( otherData[i2] ) >= tolerance ) { selfData[i1] /= otherData[i2] ; } else { selfData[i1] = value ; } }
# define _SelfIsModified
void RealIterator_Divide ( Iterator *self, Real *selfData, Iterator *other, const Real *otherData, const Real tolerance, const Real value, Status *status )
{
# include "Iterator_BinaryOperation.i"
}
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _InnerLoopOperation         result += cblas_ddot ( extent, &selfData[first1], stride1, &otherData[first2], stride2 ) ;
# define _Operation                  { result += ( selfData[i1] * otherData[i2] ) ; }
# define _Reduction
# define _ReductionResultInitializer 0.0e+00
# define _ReductionResultType        Real
Real RealIterator_Dot ( Iterator *self, const Real *selfData, Iterator *other, const Real *otherData, Status *status )
{
# include "Iterator_BinaryOperation.i"
}
# undef _InnerLoopOperation
# undef _Operation
# undef _Reduction
# undef _ReductionResultInitializer
# undef _ReductionResultType
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { selfData[i1] *= otherData[i2] ; }
# define _SelfIsModified
void RealIterator_Multiply ( Iterator *self, Real *selfData, Iterator *other, const Real *otherData, Status *status )
{
# include "Iterator_BinaryOperation.i"
}
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _InnerLoopOperation cblas_dswap ( extent, &selfData[first1], stride1, &otherData[first2], stride2 ) ;
# define _Locals             auto Real t ;
# define _Operation          { t = otherData[i2] ; otherData[i2] = selfData[i1] ; selfData[i1] = t ; }
# define _OtherIsModified
# define _SelfIsModified
void RealIterator_Swap ( Iterator *self, Real *selfData, Iterator *other, Real *otherData, Status *status )
{
# include "Iterator_BinaryOperation.i"
}
# undef _InnerLoopOperation
# undef _Locals
# undef _Operation
# undef _OtherIsModified
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
