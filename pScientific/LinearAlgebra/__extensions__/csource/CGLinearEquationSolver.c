/*==================================================================================================================================
! . This module implements a preconditioned conjugate-gradient linear equation solver.
!=================================================================================================================================*/

# include <stdio.h>
# include <stdlib.h>

# include "CGLinearEquationSolver.h"
# include "Memory.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define Default_ConvergenceMode     1
# define Default_MaximumIterations 500
# define Default_ErrorTolerance    1.0e-10

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Solver. */
static Boolean CGLinearEquationSolver_IsConverged ( const Integer convergenceMode, const Real errorTolerance, const Real rNorm2, const Real hNorm2, const Real bNorm2 ) ;

/*==================================================================================================================================
! . Report.
!=================================================================================================================================*/
void CGLinearEquationSolverReport_Initialize ( CGLinearEquationSolverReport *self )
{
    if ( self != NULL )
    {
        self->finalResidual   = 0.0e+00 ;
        self->initialResidual = 0.0e+00 ;
        self->isConverged     = False   ;
        self->iterations      = 0       ;
        self->rhsNorm2        = 0.0e+00 ;
    }
}

/*==================================================================================================================================
! . Solver.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
CGLinearEquationSolver *CGLinearEquationSolver_Allocate ( void )
{
    CGLinearEquationSolver *self = Memory_AllocateType ( CGLinearEquationSolver ) ;
    if ( self != NULL )
    {
        self->convergenceMode   = Default_ConvergenceMode   ;
        self->maximumIterations = Default_MaximumIterations ;
        self->errorTolerance    = Default_ErrorTolerance    ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
CGLinearEquationSolver *CGLinearEquationSolver_Clone ( const CGLinearEquationSolver *self )
{
    CGLinearEquationSolver *new = NULL ;
    if ( self != NULL )
    {
        new = CGLinearEquationSolver_Allocate ( ) ;
        if ( new != NULL )
        {
            new->convergenceMode   = self->convergenceMode   ;
            new->maximumIterations = self->maximumIterations ;
            new->errorTolerance    = self->errorTolerance    ;
        }
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CGLinearEquationSolver_Deallocate ( CGLinearEquationSolver **self )
{
    if ( (*self) != NULL ) Memory_Deallocate ( (*self) ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for convergence.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Boolean CGLinearEquationSolver_IsConverged ( const Integer convergenceMode, const Real errorTolerance, const Real rNorm2, const Real hNorm2, const Real bNorm2 )
{
    Real stepTest = 0.0e+00 ;
    switch ( convergenceMode )
    {
        case 1: stepTest = rNorm2          ; break ;
        case 2: stepTest = rNorm2 / bNorm2 ; break ;
        case 3: stepTest = hNorm2          ; break ;
        case 4: stepTest = hNorm2 / bNorm2 ; break ;
    }
    return ( stepTest <= errorTolerance ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Solver for SPD matrices with or without preconditioning.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . H can be R if there is no preconditioning. */

void CGLinearEquationSolver_PCGSolver ( const CGLinearEquationSolver *self, CGLinearEquationSolverState *state, CGLinearEquationSolverReport *report )
{
    CGLinearEquationSolverReport_Initialize ( report ) ;
    if ( ( self != NULL ) && ( state != NULL ) )
    {
        auto Boolean      doPreconditioning, isConverged ;
        auto Integer      iterations  ;
        auto Real         alpha, beta, bNorm2, denominator, h0Norm2, hNorm2, oldRH, r0Norm2, rDotH, rNorm2 ;
        auto RealArray1D *b, *h, *r, *x ;
        auto CGLinearEquationSolverTarget *target ;
        auto void                         *object ;

        /* . Aliases. */
        target = state->target ;
        b      = state->b ;
        h      = state->h ;
        r      = state->r ;
        x      = state->x ;
        object = target->object ;

        /* . Initialization. */
        doPreconditioning = ( target->ApplyPreconditioner != NULL ) ;
        bNorm2            = RealArray1D_Norm2 ( b ) ;
        iterations        = 0 ;

        /* . Denominator for convergence checks. */
        denominator = 1.0e+00 ;
             if ( self->convergenceMode == 2 ) denominator = bNorm2 ;
        else if ( self->convergenceMode == 4 )
        {
            if ( doPreconditioning )
            {
                (*target->ApplyPreconditioner) ( object, b, r ) ;
                denominator = RealArray1D_Norm2 ( r ) ;
            }
            else denominator = bNorm2 ;
        }

        /* . Compute initial residual r = b - A*x. */
        (*target->ApplyMatrix)     ( object, x, r ) ;
        RealArray1D_Add ( r, -1.0e+00, b, NULL ) ;
        RealArray1D_Scale          ( r, -1.0e+00  ) ;
        r0Norm2 = RealArray1D_Norm2 ( r ) ;
        rNorm2  = r0Norm2 ;

        /* . Preconditioning h = p * r (in b). */
        if ( doPreconditioning )
        {
            (*target->ApplyPreconditioner) ( object, r, b ) ;
            h0Norm2 = RealArray1D_Norm2 ( b ) ;
        }
        else
        {
            RealArray1D_CopyTo  ( r, b, NULL ) ;
            h0Norm2 = r0Norm2 ;
        }

        /* . Initial convergence check. */
        isConverged = CGLinearEquationSolver_IsConverged ( self->convergenceMode, self->errorTolerance, r0Norm2, h0Norm2, denominator );
        if ( ! isConverged )
        {
            rDotH = RealArray1D_Dot ( r, b, NULL ) ;
            for ( iterations = 1 ; iterations <= self->maximumIterations ; iterations++ )
            {
                /* . New x. */
                (*target->ApplyMatrix) ( object, b, h ) ;
                alpha = rDotH / RealArray1D_Dot ( h, b, NULL ) ;
                RealArray1D_Add ( x, alpha, b, NULL ) ;

                /* . New r. */
                RealArray1D_Add ( r, -alpha, h, NULL ) ;
                rNorm2 = RealArray1D_Norm2 ( r ) ;

                /* . New h. */
                if ( doPreconditioning )
                {
                    (*target->ApplyPreconditioner) ( object, r, h ) ;
                    hNorm2 = RealArray1D_Norm2 ( h ) ;
                }
                else
                {
                    RealArray1D_CopyTo  ( r, h, NULL ) ;
                    hNorm2 = rNorm2 ;
                }

                /* . Check for termination. */
                isConverged = CGLinearEquationSolver_IsConverged ( self->convergenceMode, self->errorTolerance, rNorm2, hNorm2, denominator ) ;
                if ( ( isConverged ) || ( iterations >= self->maximumIterations ) ) break ;

                /* . New p. */
                oldRH = rDotH ;
                rDotH = RealArray1D_Dot ( r, h, NULL ) ;
                beta  = rDotH / oldRH ;
                RealArray1D_Scale ( b, beta ) ;
                RealArray1D_Add ( b, 1.0e+00, h, NULL ) ;
            }
        }

        /* . Finish up. */
        if ( report != NULL )
        {
            report->finalResidual   = rNorm2      ;
            report->initialResidual = r0Norm2     ;
            report->isConverged     = isConverged ;
            report->iterations      = iterations  ;
            report->rhsNorm2        = bNorm2      ;
        }
    }
}

/*==================================================================================================================================
! . State.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
CGLinearEquationSolverState *CGLinearEquationSolverState_Allocate ( void )
{
    CGLinearEquationSolverState *self = Memory_AllocateType ( CGLinearEquationSolverState ) ;
    if ( self != NULL )
    {
        self->b      = NULL ;
        self->h      = NULL ;
        self->r      = NULL ;
        self->x      = NULL ;
        self->target = NULL ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CGLinearEquationSolverState_Deallocate ( CGLinearEquationSolverState **self )
{
    if ( (*self) != NULL )
    {
        RealArray1D_Deallocate ( &((*self)->h) ) ;
        RealArray1D_Deallocate ( &((*self)->r) ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Setup a state given a target.
!---------------------------------------------------------------------------------------------------------------------------------*/
CGLinearEquationSolverState *CGLinearEquationSolverState_SetupFromTarget ( CGLinearEquationSolverTarget *target, Status *status )
{
    CGLinearEquationSolverState *self = NULL ;
    if ( target != NULL )
    {
        auto Integer n = 0 ;
        if ( target->rhs != NULL ) n = View1D_Extent ( target->rhs ) ;
        if ( n > 0 )
        {
            self = CGLinearEquationSolverState_Allocate ( ) ;
            if ( self != NULL )
            {
                /* . Assignment. */
                self->b      = target->rhs      ;
                self->x      = target->solution ;
                self->target = target           ;
                /* . Allocate space. */
                self->h = RealArray1D_AllocateWithExtent ( n, status ) ;
                self->r = RealArray1D_AllocateWithExtent ( n, status ) ;
                if ( ( self->h == NULL ) || ( self->r == NULL ) ) CGLinearEquationSolverState_Deallocate ( &self ) ;
            }
            if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
        }
    }
    return self ;
}

/*==================================================================================================================================
! . Target.
!=================================================================================================================================*/
void CGLinearEquationSolverTarget_Initialize ( CGLinearEquationSolverTarget *self )
{
    if ( self != NULL )
    {
        self->rhs                 = NULL ;
        self->solution            = NULL ;
        self->object              = NULL ;
        self->ApplyMatrix         = NULL ;
        self->ApplyPreconditioner = NULL ;
    }
}
