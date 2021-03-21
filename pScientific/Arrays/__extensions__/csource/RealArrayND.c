/*==================================================================================================================================
! . N-D real arrays.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "NumericalMacros.h"
# include "RealArrayND.h"
# include "RealIterator.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Basic definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _ArrayDataFormat          "%20.10f"
# define _ArrayDataPerLine         6
# define _ArrayDataType            Real
# define _ArrayDataTypeInitializer 0.0e+00
# define _UseCBLAS
# include "ArrayND_Body.i"
# undef _ArrayDataFormat
# undef _ArrayDataPerLine
# undef _ArrayDataType
# undef _ArrayDataTypeInitializer
# undef _UseCBLAS
