/*==================================================================================================================================
! . Boolean memory block.
!=================================================================================================================================*/

# include "BooleanBlock.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Basic functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _CoreDataType              Boolean
# define _CoreDataTypeFormat        "%2u"
# define _CoreDataTypeInitializer   False
# define _UnaryNegationOperator !
# include "TypedMemoryBlock_Body.i"
# undef  _CoreDataType
# undef  _CoreDataTypeFormat
# undef  _CoreDataTypeInitializer
# undef  _UnaryNegationOperator
