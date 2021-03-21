# ifndef _CISPARSESOLVER
# define _CISPARSESOLVER

# include "cprimme.h"
# include "Integer.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void CISparseSolver_ApplyMatrix         ( void *xVoid, void *yVoid, Integer *blockSize, primme_params *primme ) ;
extern void CISparseSolver_ApplyPreconditioner ( void *xVoid, void *yVoid, Integer *blockSize, primme_params *primme ) ;

# endif
