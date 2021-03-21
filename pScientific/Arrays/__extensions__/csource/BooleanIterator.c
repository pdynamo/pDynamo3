/*==================================================================================================================================
! . Boolean iterator operations.
!=================================================================================================================================*/

# include <stdlib.h>

# include "BooleanIterator.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Unary.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation                  { if ( ! selfData[i] ) return False ; }
# define _Reduction
# define _ReductionResultInitializer True
# define _ReductionResultType        Boolean
Boolean BooleanIterator_All ( Iterator *self, const Boolean *selfData )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _Reduction
# undef _ReductionResultInitializer
# undef _ReductionResultType
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation                  { if ( selfData[i] ) return True ; }
# define _Reduction
# define _ReductionResultInitializer False
# define _ReductionResultType        Boolean
Boolean BooleanIterator_Any ( Iterator *self, const Boolean *selfData )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _Reduction
# undef _ReductionResultInitializer
# undef _ReductionResultType
/*--------------------------------------------------------------------------------------------------------------------------------*/
Integer BooleanIterator_CountFalse ( Iterator *self, const Boolean *selfData )
{
    Integer count = 0 ;
    if ( self != NULL ) count = ( self->size - BooleanIterator_CountTrue ( self, selfData ) ) ;
    return count ;
}
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation                  { if ( selfData[i] ) result += 1 ; }
# define _Reduction
# define _ReductionResultInitializer 0
# define _ReductionResultType        Integer
Integer BooleanIterator_CountTrue ( Iterator *self, const Boolean *selfData )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _Reduction
# undef _ReductionResultInitializer
# undef _ReductionResultType
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { selfData[i] = ! selfData[i] ; }
# define _SelfIsModified
void BooleanIterator_Not ( Iterator *self, Boolean *selfData, Status *status )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { selfData[i] = value ; }
# define _SelfIsModified
void BooleanIterator_Set ( Iterator *self, Boolean *selfData, const Boolean value, Status *status )
{
# include "Iterator_UnaryOperation.i"
}
# undef _Operation
# undef _SelfIsModified
/*----------------------------------------------------------------------------------------------------------------------------------
! . Binary.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { selfData[i1] = ( selfData[i1] && otherData[i2] ) ; }
# define _SelfIsModified
void BooleanIterator_And ( Iterator *self, Boolean *selfData, Iterator *other, const Boolean *otherData, Status *status )
{
# include "Iterator_BinaryOperation.i"
}
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { otherData[i2] = selfData[i1] ; }
# define _OtherIsModified
void BooleanIterator_CopyTo ( Iterator *self, const Boolean *selfData, Iterator *other, Boolean *otherData, Status *status )
{
# include "Iterator_BinaryOperation.i"
}
# undef _Operation
# undef _OtherIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { selfData[i1] = ( selfData[i1] || otherData[i2] ) ; }
# define _SelfIsModified
void BooleanIterator_Or ( Iterator *self, Boolean *selfData, Iterator *other, const Boolean *otherData, Status *status )
{
# include "Iterator_BinaryOperation.i"
}
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Operation { selfData[i1] = ( selfData[i1] != otherData[i2] ) ; } /* . OK as booleans either 0 or 1. */
# define _SelfIsModified
void BooleanIterator_Xor ( Iterator *self, Boolean *selfData, Iterator *other, const Boolean *otherData, Status *status )
{
# include "Iterator_BinaryOperation.i"
}
# undef _Operation
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
# define _Locals          auto Boolean t ;
# define _Operation       { t = otherData[i2] ; otherData[i2] = selfData[i1] ; selfData[i1] = t ; }
# define _OtherIsModified
# define _SelfIsModified
void BooleanIterator_Swap ( Iterator *self, Boolean *selfData, Iterator *other, Boolean *otherData, Status *status )
{
# include "Iterator_BinaryOperation.i"
}
# undef _Locals
# undef _Operation
# undef _OtherIsModified
# undef _SelfIsModified
/*--------------------------------------------------------------------------------------------------------------------------------*/
