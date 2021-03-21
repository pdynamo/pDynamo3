/*==================================================================================================================================
! . Integer memory block.
!=================================================================================================================================*/

# include "IntegerBlock.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Basic functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _CoreDataType            Integer
# define _CoreDataTypeFormat      "%15d"
# define _CoreDataTypeInitializer 0
# define _ExcludeNegate
# include "TypedMemoryBlock_Body.i"
# undef  _CoreDataType
# undef  _CoreDataTypeFormat
# undef  _CoreDataTypeInitializer
# undef  _ExcludeNegate
