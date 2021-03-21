/*==================================================================================================================================
! . 2-D boolean arrays.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "BooleanArray2D.h"
# include "NumericalMacros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Basic definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _ArrayDataFormat          "%5u"
# define _ArrayDataPerLine         12
# define _ArrayDataType            Boolean
# define _ArrayDataTypeInitializer False
# define _NoNumeric
# include "Array2D_Body.i"
# undef _ArrayDataFormat
# undef _ArrayDataPerLine
# undef _ArrayDataType
# undef _ArrayDataTypeInitializer
# undef _NoNumeric
