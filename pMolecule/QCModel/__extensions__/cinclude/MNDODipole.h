# ifndef _MNDODIPOLE
# define _MNDODIPOLE

# include "Coordinates3.h"
# include "IntegerArray1D.h"
# include "MNDOParametersContainer.h"
# include "Status.h"
# include "SymmetricMatrix.h"
# include "Vector3.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void MNDO_DipoleIntegrals ( const MNDOParametersContainer *parameters   ,
                                   const IntegerArray1D          *basisIndices ,
                                   const Coordinates3            *coordinates3 ,
                                   const Vector3                 *center       ,
                                         SymmetricMatrix         *dX           ,
                                         SymmetricMatrix         *dY           ,
                                         SymmetricMatrix         *dZ           ,
                                         Status                  *status       ) ;

# endif
