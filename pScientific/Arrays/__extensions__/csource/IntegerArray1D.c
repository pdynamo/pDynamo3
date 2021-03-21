/*==================================================================================================================================
! . 1-D integer arrays.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "IntegerArray1D.h"
# include "NumericalMacros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Basic definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _ArrayDataFormat          "%20d"
# define _ArrayDataPerLine         6
# define _ArrayDataType            Integer
# define _ArrayDataTypeInitializer 0
# include "Array1D_Body.i"
# undef _ArrayDataFormat
# undef _ArrayDataPerLine
# undef _ArrayDataType
# undef _ArrayDataTypeInitializer
