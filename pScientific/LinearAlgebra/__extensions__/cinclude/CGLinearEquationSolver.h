# ifndef _CGLINEAREQUATIONSOLVER
# define _CGLINEAREQUATIONSOLVER

# include "Boolean.h"
# include "Integer.h"
# include "Real.h"
# include "RealArray1D.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The report type. */
typedef struct {
    Boolean isConverged     ;
    Integer iterations      ;
    Real    finalResidual   ;
    Real    initialResidual ;
    Real    rhsNorm2        ;
} CGLinearEquationSolverReport ;

/* . The target type. */
typedef struct {
    RealArray1D *rhs      ;
    RealArray1D *solution ;
    void        *object   ;
    void ( * ApplyMatrix         ) ( void *target, RealArray1D *x, RealArray1D *y ) ;
    void ( * ApplyPreconditioner ) ( void *target, RealArray1D *x, RealArray1D *y ) ;
} CGLinearEquationSolverTarget ;

/* . The state type. */
typedef struct {
    RealArray1D *b ;
    RealArray1D *h ;
    RealArray1D *r ;
    RealArray1D *x ;
    CGLinearEquationSolverTarget *target ;
} CGLinearEquationSolverState ;

/* . The solver type. */
typedef struct {
    Integer convergenceMode   ;
    Integer maximumIterations ;
    Real    errorTolerance    ;
} CGLinearEquationSolver ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Report. */
extern void CGLinearEquationSolverReport_Initialize ( CGLinearEquationSolverReport *self ) ;

/* . Solver. */
extern CGLinearEquationSolver      *CGLinearEquationSolver_Allocate   ( void ) ;
extern CGLinearEquationSolver      *CGLinearEquationSolver_Clone      ( const CGLinearEquationSolver  *self ) ;
extern void                         CGLinearEquationSolver_Deallocate (       CGLinearEquationSolver **self ) ;
extern void                         CGLinearEquationSolver_PCGSolver  ( const CGLinearEquationSolver  *self, CGLinearEquationSolverState *state, CGLinearEquationSolverReport *report ) ;

/* . State. */
extern CGLinearEquationSolverState *CGLinearEquationSolverState_Allocate        ( void ) ;
extern void                         CGLinearEquationSolverState_Deallocate      (       CGLinearEquationSolverState **self ) ;
extern CGLinearEquationSolverState *CGLinearEquationSolverState_SetupFromTarget ( CGLinearEquationSolverTarget *target, Status *status ) ;

/* . Target. */
extern void CGLinearEquationSolverTarget_Initialize ( CGLinearEquationSolverTarget *self ) ;

# endif
