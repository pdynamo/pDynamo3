# ifndef _BOOLEANUTILITIES
# define _BOOLEANUTILITIES

# include "Boolean.h"
# include "Integer.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Raw array definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _CoreDataType Boolean
# define _ExcludeSort
# include "RawBlock_Header.i"
# undef  _CoreDataType
# undef  _ExcludeSort

# endif
