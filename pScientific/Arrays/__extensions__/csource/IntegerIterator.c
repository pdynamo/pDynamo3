/*==================================================================================================================================
! . Integer iterator operations.
!=================================================================================================================================*/

# include <stdlib.h>
# include <stdlib.h>

# include "IntegerIterator.h"
# include "NumericalMacros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Unary.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation                  { result = Maximum ( result, abs ( selfData[i] ) ) ; }
# define _Reduction
# define _ReductionResultInitializer 0
# define _ReductionResultType        Integer
Integer IntegerIterator_AbsoluteMaximum ( Iterator *self, const Integer *selfData )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _Reduction
# undef _ReductionResultInitializer
# undef _ReductionResultType
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation                  { if ( abs ( selfData[i] ) <= tolerance ) result += 1 ; }
# define _Reduction
# define _ReductionResultInitializer 0
# define _ReductionResultType        Integer
Integer IntegerIterator_CountSmall ( Iterator *self, const Integer *selfData, const Integer tolerance )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _Reduction
# undef _ReductionResultInitializer
# undef _ReductionResultType
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation                  { result += ( selfData[i] * selfData[i] ) ; }
# define _Reduction
# define _ReductionResultInitializer 0
# define _ReductionResultType        Integer
Integer IntegerIterator_DotSelf ( Iterator *self, const Integer *selfData )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _Reduction
# undef _ReductionResultInitializer
# undef _ReductionResultType
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { if ( selfData[i] >= tolerance ) selfData[i] = value ; }
# define _SelfIsModified
extern void IntegerIterator_FilterGreaterThan ( Iterator *self, Integer *selfData, const Integer tolerance, const Integer value, Status *status )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { if ( selfData[i] <= tolerance ) selfData[i] = value ; }
# define _SelfIsModified
extern void IntegerIterator_FilterLessThan ( Iterator *self, Integer *selfData, const Integer tolerance, const Integer value, Status *status )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { if ( abs ( selfData[i] ) <= tolerance ) selfData[i] = value ; }
# define _SelfIsModified
extern void IntegerIterator_FilterSmall ( Iterator *self, Integer *selfData, const Integer tolerance, const Integer value, Status *status )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { selfData[i] += value ; }
# define _SelfIsModified
void IntegerIterator_Increment ( Iterator *self, Integer *selfData, const Integer value, Status *status )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation                  { result = Maximum ( result, selfData[i] ) ; }
# define _Reduction
# define _ReductionResultInitializer Integer_Smallest
# define _ReductionResultType        Integer
Integer IntegerIterator_Maximum ( Iterator *self, const Integer *selfData )
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
# define _ReductionResultInitializer Integer_Largest
# define _ReductionResultType        Integer
Integer IntegerIterator_Minimum ( Iterator *self, const Integer *selfData )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _Reduction
# undef _ReductionResultInitializer
# undef _ReductionResultType
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation                  { result *= selfData[i] ; }
# define _Reduction
# define _ReductionResultInitializer 1
# define _ReductionResultType        Integer
Integer IntegerIterator_Product ( Iterator *self, const Integer *selfData )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _Reduction
# undef _ReductionResultInitializer
# undef _ReductionResultType
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { selfData[i] *= value ; }
# define _SelfIsModified
void IntegerIterator_Scale ( Iterator *self, Integer *selfData, const Integer value, Status *status )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { selfData[i] = value ; }
# define _SelfIsModified
void IntegerIterator_Set ( Iterator *self, Integer *selfData, const Integer value, Status *status )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
Real IntegerIterator_Sparsity ( Iterator *self, const Integer *selfData, const Integer tolerance )
{
    Real result = 0.0 ;
    if ( self != NULL ) { result = ( 100.0 * ( Real ) IntegerIterator_CountSmall ( self, selfData, tolerance ) / ( Real ) self->size ) ; }
    return result ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation                  { result += selfData[i] ; }
# define _Reduction
# define _ReductionResultInitializer 0
# define _ReductionResultType        Integer
Integer IntegerIterator_Sum ( Iterator *self, const Integer *selfData )
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
# define _Operation { selfData[i1] += ( scale * otherData[i2] ) ; }
# define _SelfIsModified
void IntegerIterator_Add ( Iterator *self, Integer *selfData, Iterator *other, const Integer *otherData, const Integer scale, Status *status )
{
# include "Iterator_BinaryOperation.i"
}
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { otherData[i2] = selfData[i1] ; }
# define _OtherIsModified
void IntegerIterator_CopyTo ( Iterator *self, const Integer *selfData, Iterator *other, Integer *otherData, Status *status )
{
# include "Iterator_BinaryOperation.i"
}
# undef _Operation
# undef _OtherIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation                  { result += ( selfData[i1] * otherData[i2] ) ; }
# define _Reduction
# define _ReductionResultInitializer 0
# define _ReductionResultType        Integer
Integer IntegerIterator_Dot ( Iterator *self, const Integer *selfData, Iterator *other, const Integer *otherData, Status *status )
{
# include "Iterator_BinaryOperation.i"
}
# undef _Operation
# undef _Reduction
# undef _ReductionResultInitializer
# undef _ReductionResultType
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { selfData[i1] *= otherData[i2] ; }
# define _SelfIsModified
void IntegerIterator_Multiply ( Iterator *self, Integer *selfData, Iterator *other, const Integer *otherData, Status *status )
{
# include "Iterator_BinaryOperation.i"
}
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Locals          auto Integer t ;
# define _Operation       { t = otherData[i2] ; otherData[i2] = selfData[i1] ; selfData[i1] = t ; }
# define _OtherIsModified
# define _SelfIsModified
void IntegerIterator_Swap ( Iterator *self, Integer *selfData, Iterator *other, Integer *otherData, Status *status )
{
# include "Iterator_BinaryOperation.i"
}
# undef _Locals
# undef _Operation
# undef _OtherIsModified
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
