# ifndef _JDEIGENVALUESOLVER
# define _JDEIGENVALUESOLVER

# include "cprimme.h"
# include "Boolean.h"
# include "Integer.h"
# include "IntegerArray1D.h"
# include "Real.h"
# include "RealArray1D.h"
# include "RealArray2D.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The report type. */
typedef struct {
    Boolean isConverged        ;
    Boolean solutionChecked    ;
    Integer convergedPairs     ;
    Integer numberMatrixVectorMultiplications  ;
    Integer returnCode         ;
    Real    eigenvalueError    ;
    Real    eigenvectorError   ;
    Real    normalizationError ;
} JDEigenvalueSolverReport ;

/* . The target type. */
typedef struct {
    RealArray1D *eigenvalues  ;
    RealArray2D *eigenvectors ;
    void        *object       ;
    void ( * ApplyMatrix         ) ( void *x, void *y, int *blockSize, struct primme_params *primme ) ;
    void ( * ApplyPreconditioner ) ( void *x, void *y, int *blockSize, struct primme_params *primme ) ;
} JDEigenvalueSolverTarget ;

/* . The state type. */
typedef struct {
    IntegerArray1D           *iWork         ;
    RealArray1D              *eigenvalues   ;
    RealArray1D              *residualNorms ;
    RealArray1D              *rWork         ;
    RealArray2D              *eigenvectors  ;
    JDEigenvalueSolverTarget *target        ;
    struct primme_params      primme        ;
} JDEigenvalueSolverState ;

/* . The solver type. */
typedef struct {
    Boolean usePreconditioning ;
    Integer maximumMatrixVectorMultiplications ;
    Integer printLevel     ;
    Real    errorTolerance ;
} JDEigenvalueSolver ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Report. */
void JDEigenvalueSolverReport_Initialize ( JDEigenvalueSolverReport *self ) ;

/* . Solver. */
extern JDEigenvalueSolver      *JDEigenvalueSolver_Allocate      ( void ) ;
extern void                     JDEigenvalueSolver_CheckSolution ( const JDEigenvalueSolver  *self, JDEigenvalueSolverState *state, RealArray1D *referenceEigenvalues, JDEigenvalueSolverReport *report ) ;
extern JDEigenvalueSolver      *JDEigenvalueSolver_Clone         ( const JDEigenvalueSolver  *self ) ;
extern void                     JDEigenvalueSolver_CopyTo        ( const JDEigenvalueSolver  *self, JDEigenvalueSolver *other ) ;
extern void                     JDEigenvalueSolver_Deallocate    (       JDEigenvalueSolver **self ) ;
extern void                     JDEigenvalueSolver_Initialize    (       JDEigenvalueSolver  *self ) ;
extern void                     JDEigenvalueSolver_Solve         ( const JDEigenvalueSolver  *self, JDEigenvalueSolverState *state, JDEigenvalueSolverReport *report ) ;

/* . State. */
extern JDEigenvalueSolverState *JDEigenvalueSolverState_Allocate        ( void ) ;
extern void                     JDEigenvalueSolverState_Deallocate      ( JDEigenvalueSolverState **self ) ;
extern JDEigenvalueSolverState *JDEigenvalueSolverState_SetupFromTarget ( JDEigenvalueSolverTarget *target, Status *status ) ;

/* . Target. */
void JDEigenvalueSolverTarget_Initialize ( JDEigenvalueSolverTarget *self ) ;

# endif
