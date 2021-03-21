/*==================================================================================================================================
! . N-D integer arrays.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "IntegerArrayND.h"
# include "IntegerIterator.h"
# include "NumericalMacros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Basic definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _ArrayDataFormat          "%20d"
# define _ArrayDataPerLine         6
# define _ArrayDataType            Integer
# define _ArrayDataTypeInitializer 0
# include "ArrayND_Body.i"
# undef _ArrayDataFormat
# undef _ArrayDataPerLine
# undef _ArrayDataType
# undef _ArrayDataTypeInitializer
