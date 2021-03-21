# ifndef _DENSEMATRIXPOWER
# define _DENSEMATRIXPOWER

# include "Boolean.h"
# include "Real.h"
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
extern void SymmetricMatrix_InversePower (       SymmetricMatrix *self          , 
                                           const Boolean          preserveInput ,
                                           const Real             power         ,
                                           const Real             tolerance     ,
                                                 SymmetricMatrix *result        ,
                                                 Status          *status        ) ;
extern void SymmetricMatrix_Power        (       SymmetricMatrix *self          , 
                                           const Boolean          preserveInput ,
                                           const Real             power         ,
                                                 SymmetricMatrix *result        ,
                                                 Status          *status        ) ;
# endif
