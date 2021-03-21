/*==================================================================================================================================
! . Boolean utilities.
!=================================================================================================================================*/

# include "BooleanUtilities.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Array functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _CoreDataType            Boolean
# define _CoreDataTypeFormat      "%2u"
# define _CoreDataTypeInitializer False
# define _ExcludeSort
# include "RawBlock_Body.i"
# undef  _CoreDataType
# undef  _CoreDataTypeFormat
# undef  _CoreDataTypeInitializer
# undef  _ExcludeSort
