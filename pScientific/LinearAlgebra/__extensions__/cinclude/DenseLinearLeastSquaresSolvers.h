# ifndef _DENSELINEARLEASTSQUARESSOLVERS
# define _DENSELINEARLEASTSQUARESSOLVERS

# include "Boolean.h"
# include "Integer.h"
# include "Real.h"
# include "RealArray1D.h"
# include "RealArray2D.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void LinearLeastSquaresSVDSolve ( RealArray2D *self              ,
                                         RealArray1D *rhs               ,
                                         Boolean      preserveInput     ,
                                         Real        *relativeTolerance ,
                                         RealArray1D *solution          ,
                                         Integer     *rank              ,
                                         Real        *condition         ,
                                         Status      *status            ) ;
# endif
