# ifndef _DENSEMATRIXPOWER
# define _DENSEMATRIXPOWER

# include "Boolean.h"
# include "Integer.h"
# include "Real.h"
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
extern void Matrix_PseudoInverse         (       RealArray2D     *self              ,
                                           const Boolean          preserveInput     ,
                                                 Real            *relativeTolerance ,
                                                 RealArray2D     *inverse           ,
                                                 Integer         *rank              ,
                                                 Real            *condition         ,
                                                 Status          *status            ) ;
extern void SymmetricMatrix_InversePower (       SymmetricMatrix *self              , 
                                           const Boolean          preserveInput     ,
                                           const Real             power             ,
                                           const Real             tolerance         ,
                                                 SymmetricMatrix *result            ,
                                                 Status          *status            ) ;
extern void SymmetricMatrix_Power        (       SymmetricMatrix *self              , 
                                           const Boolean          preserveInput     ,
                                           const Real             power             ,
                                                 SymmetricMatrix *result            ,
                                                 Status          *status            ) ;
# endif
