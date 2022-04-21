/* . Definitions for Gaussian nuclei. */
# ifndef _GAUSSIANNUCLEUS
# define _GAUSSIANNUCLEUS

# include <math.h>

# include "RealArray1D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Default width parameters (atomic units). */
# define _DefaultNuclearWidth  1.0e-4
# define _DefaultNuclearWidthE ( 4.0e+00 / ( _DefaultNuclearWidth * _DefaultNuclearWidth * M_PI ) )
# define _DefaultNuclearWidthN sqrt ( pow ( ( _DefaultNuclearWidthE / M_PI ), 3 ) )

/* . Macros to return defaults if necessary. */
# define _GetWidthE( widthsE, i ) ( (widthsE) == NULL ? _DefaultNuclearWidthE : Array1D_Item ( widthsE, i ) )
# define _GetWidthN( widthsN, i ) ( (widthsN) == NULL ? _DefaultNuclearWidthN : Array1D_Item ( widthsN, i ) )

# endif
