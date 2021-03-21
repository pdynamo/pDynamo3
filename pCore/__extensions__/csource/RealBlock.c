/*==================================================================================================================================
! . Real memory block.
!=================================================================================================================================*/

# include "RealBlock.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Basic functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _CoreDataType              Real
# define _CoreDataTypeFormat        "%15.10f"
# define _CoreDataTypeInitializer   0.0e+00
# define _UnaryNegationOperator -
# include "TypedMemoryBlock_Body.i"
# undef  _CoreDataType
# undef  _CoreDataTypeFormat
# undef  _CoreDataTypeInitializer
# undef  _UnaryNegationOperator
