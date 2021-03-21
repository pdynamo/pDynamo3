# ifndef _DENSEEIGENVALUESOLVERS
# define _DENSEEIGENVALUESOLVERS

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
extern void SymmetricMatrix_EigenvaluesSolve (       SymmetricMatrix *self          ,
                                               const Boolean          preserveInput ,
                                               const Integer          lower         ,
                                               const Integer          upper         ,
                                                     RealArray1D     *eigenValues   ,
                                                     RealArray2D     *eigenVectors  ,
                                               const Boolean          isColumnMajor ,
                                                     Status          *status        ) ;
# endif
