# ifndef _INTEGERBLOCK
# define _INTEGERBLOCK

# include "Integer.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Basic definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _CoreDataType Integer
# define _ExcludeNegate
# include "TypedMemoryBlock_Header.i"
# undef  _CoreDataType
# undef  _ExcludeNegate

# endif
