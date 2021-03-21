# ifndef _DENSELINEAREQUATIONSOLVERS
# define _DENSELINEAREQUATIONSOLVERS

# include "Boolean.h"
# include "Integer.h"
# include "Real.h"
# include "RealArray1D.h"
# include "RealArray2D.h"
# include "Status.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void SquareMatrix_LinearEquationsSolve    ( RealArray2D     *self          ,
                                                   RealArray1D     *rhs           ,
                                                   Boolean          preserveInput ,
                                                   RealArray1D     *solution      ,
                                                   Status          *status        ) ;
extern void SymmetricMatrix_LinearEquationsSolve ( SymmetricMatrix *self          ,
                                                   RealArray1D     *rhs           ,
                                                   RealArray1D     *solution      ,
                                                   Status          *status        ) ;
# endif
