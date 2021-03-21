/*==================================================================================================================================
! . This module implements a Jacobi-Davidson eigenvalue solver for symmetric matrices.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "JDEigenvalueSolver.h"
# include "Memory.h"
# include "NumericalMacros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define Default_ErrorTolerance                     1.0e-12
# define Default_MaximumMatrixVectorMultiplications 100000
# define Default_PrimmeMethod                       DYNAMIC
# define Default_PrintLevel                         0
# define Default_UsePreconditioning                 True

/* . Procedures to be defined. */
/*
void (*matrixMatvec)         ( void *x, void *y, int *blockSize, struct primme_params *primme ) ;
void (*applyPreconditioner)  ( void *x, void *y, int *blockSize, struct primme_params *primme ) ;
*/

/* . Note that the eigenvectors are stored ROW-WISE! */

/*==================================================================================================================================
! . Report.
!=================================================================================================================================*/
void JDEigenvalueSolverReport_Initialize ( JDEigenvalueSolverReport *self )
{
    if ( self != NULL )
    {
        self->isConverged                       = False ;
        self->solutionChecked                   = False ;
        self->convergedPairs                    = 0 ;
        self->numberMatrixVectorMultiplications = 0 ;
        self->returnCode                        = 0 ;
        self->eigenvalueError                   = 0.0e+00 ;
        self->eigenvectorError                  = 0.0e+00 ;
        self->normalizationError                = 0.0e+00 ;
    }
}

/*==================================================================================================================================
! . Solver.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
JDEigenvalueSolver *JDEigenvalueSolver_Allocate ( void )
{
    JDEigenvalueSolver *self = Memory_AllocateType ( JDEigenvalueSolver ) ;
    JDEigenvalueSolver_Initialize ( self ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check the solution.
!---------------------------------------------------------------------------------------------------------------------------------*/
void JDEigenvalueSolver_CheckSolution ( const JDEigenvalueSolver *self, JDEigenvalueSolverState *state, RealArray1D *referenceEigenvalues, JDEigenvalueSolverReport *report )
{
    if ( ( self != NULL ) && ( state != NULL ) && ( report != NULL ) )
    {
        Integer     dummy = 1, i ;
        Real        deviation1 = 0.0e+00, deviation2 = 0.0e+00, deviation3 = 0.0e+00 ;
        RealArray1D row, *temporary = NULL ;

        /* . Initialization. */
        temporary = RealArray1D_AllocateWithExtent ( state->primme.n, NULL ) ;
        if ( temporary != NULL )
        {
            /* . Loop over eigenvectors. */
            for ( i = 0 ; i < state->primme.numEvals ; i++ )
            {
                /* . Get A * X - mu * X. */
                RealArray2D_RowView ( state->eigenvectors, i, False, &row, NULL ) ;
                state->target->ApplyMatrix ( Array1D_Data ( &row ), Array1D_Data ( temporary ), &dummy, &(state->primme) ) ;
                RealArray1D_Add ( temporary, - Array1D_Item ( state->eigenvalues, i ), &row, NULL ) ;
                deviation1 = Maximum ( deviation1, RealArray1D_AbsoluteMaximum ( temporary ) ) ;
                deviation2 = Maximum ( deviation2, RealArray1D_Norm2 ( &row ) - 1.0e+00 ) ;
            }

            /* . Loop over eigenvalues. */
            if ( referenceEigenvalues != NULL )
            {
                for ( i = 0 ; i < state->primme.numEvals ; i++ )
                {
                    deviation3 = Maximum ( deviation3, fabs ( Array1D_Item ( state->eigenvalues, i ) - Array1D_Item ( referenceEigenvalues, i ) ) ) ;
                }
            }

            /* . Finish up. */
            report->solutionChecked    = True ;
            report->eigenvalueError    = deviation3 ;
            report->eigenvectorError   = deviation1 ;
            report->normalizationError = deviation2 ;
            RealArray1D_Deallocate ( &temporary ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
JDEigenvalueSolver *JDEigenvalueSolver_Clone ( const JDEigenvalueSolver *self )
{
    JDEigenvalueSolver *new = NULL ;
    if ( self != NULL )
    {
        new = JDEigenvalueSolver_Allocate ( ) ;
        JDEigenvalueSolver_CopyTo ( self, new ) ;
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copying.
!---------------------------------------------------------------------------------------------------------------------------------*/
void JDEigenvalueSolver_CopyTo ( const JDEigenvalueSolver *self, JDEigenvalueSolver *other )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        other->usePreconditioning                 = self->usePreconditioning                 ;
        other->errorTolerance                     = self->errorTolerance                     ;
        other->maximumMatrixVectorMultiplications = self->maximumMatrixVectorMultiplications ;
        other->printLevel                         = self->printLevel                         ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void JDEigenvalueSolver_Deallocate ( JDEigenvalueSolver **self )
{
    if ( (*self) != NULL ) Memory_Deallocate ( (*self) ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void JDEigenvalueSolver_Initialize ( JDEigenvalueSolver *self )
{
    if ( self != NULL )
    {
        self->usePreconditioning                 = Default_UsePreconditioning                 ;
        self->errorTolerance                     = Default_ErrorTolerance                     ;
        self->maximumMatrixVectorMultiplications = Default_MaximumMatrixVectorMultiplications ;
        self->printLevel                         = Default_PrintLevel                         ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find eigenvalues and eigenvectors.
!---------------------------------------------------------------------------------------------------------------------------------*/
void JDEigenvalueSolver_Solve ( const JDEigenvalueSolver *self, JDEigenvalueSolverState *state, JDEigenvalueSolverReport *report )
{
    JDEigenvalueSolverReport_Initialize ( report ) ;
    if ( ( self != NULL ) && ( state != NULL ) )
    {
        /* . Set some options. */
        state->primme.eps        = self->errorTolerance                     ;
        state->primme.maxMatvecs = self->maximumMatrixVectorMultiplications ;
        state->primme.printLevel = self->printLevel                         ;
        if ( self->usePreconditioning ) state->primme.correctionParams.precondition = 1 ;
        else                            state->primme.correctionParams.precondition = 0 ;

        /* . Simple case which gives an error in dprimme. */
        if ( state->primme.n == 1 )
        {
            auto Integer dummy ;

            /* . Set eigenvectors. */
            RealArray2D_Set ( state->eigenvectors, 1.0e+00 ) ;

            /* . Get eigenvalues. */
            state->target->ApplyMatrix ( Array2D_Data ( state->eigenvectors ), Array1D_Data ( state->eigenvalues ), &dummy, &(state->primme) ) ;

            /* . Finish up. */
            if ( report != NULL )
            {
                report->convergedPairs                    = 1    ;
                report->isConverged                       = True ;
                report->numberMatrixVectorMultiplications = 1    ;
                report->returnCode                        = 0    ;
            }
        }
        /* . Normal cases. */
        else
        {
            auto Integer returnCode ;

            /* . Get the eigenvalues and eigenvectors. */
            returnCode = dprimme ( Array1D_Data ( state->eigenvalues ), Array2D_Data ( state->eigenvectors ), Array1D_Data ( state->residualNorms ), &state->primme ) ;

            /* . Finish up. */
            if ( report != NULL )
            {
                report->convergedPairs                    = state->primme.initSize ;
                report->isConverged                       = ( returnCode == 0 ) && ( report->convergedPairs == state->primme.numEvals ) ;
                report->numberMatrixVectorMultiplications = state->primme.stats.numMatvecs ;
                report->returnCode                        = returnCode             ;
            }
        }
    }
}

/*==================================================================================================================================
! . State.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
JDEigenvalueSolverState *JDEigenvalueSolverState_Allocate ( void )
{
    JDEigenvalueSolverState *self = Memory_AllocateType ( JDEigenvalueSolverState ) ;
    if ( self != NULL )
    {
        self->eigenvalues   = NULL ;
        self->eigenvectors  = NULL ;
        self->iWork         = NULL ;
        self->residualNorms = NULL ;
        self->rWork         = NULL ;
        self->target        = NULL ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void JDEigenvalueSolverState_Deallocate ( JDEigenvalueSolverState **self )
{
    if ( (*self) != NULL )
    {
        IntegerArray1D_Deallocate ( &((*self)->iWork        ) ) ;
        RealArray1D_Deallocate    ( &((*self)->residualNorms) ) ;
        RealArray1D_Deallocate    ( &((*self)->rWork        ) ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Setup a state given a target.
!---------------------------------------------------------------------------------------------------------------------------------*/
JDEigenvalueSolverState *JDEigenvalueSolverState_SetupFromTarget ( JDEigenvalueSolverTarget *target, Status *status )
{
    JDEigenvalueSolverState *self = NULL ;
    if ( target != NULL )
    {
        auto Integer m = 0, n = 0, v = 0 ;
        m = View1D_Extent  ( target->eigenvalues  ) ; /* . Number of eigenvalues. */
        n = View2D_Columns ( target->eigenvectors ) ; /* . Size of problem. */
        v = View2D_Rows    ( target->eigenvectors ) ; /* . Number of eigenvalues. */
        if ( m != v ) v = 0 ;
        if ( ( n > 0 ) && ( v > 0 ) )
        {
            self = JDEigenvalueSolverState_Allocate ( ) ;
            if ( self != NULL )
            {
                /* . Primme data. */
                /* . Basic initialization. */
                primme_initialize( &(self->primme) ) ;
                self->primme.n                   = n ; /* . Matrix size. */
                self->primme.numEvals            = v ; /* . Number of eigenvectors. */
                self->primme.matrixMatvec        = target->ApplyMatrix         ;
                self->primme.applyPreconditioner = target->ApplyPreconditioner ;
                self->primme.matrix              = ( void * ) self ;
                self->primme.preconditioner      = NULL            ;
                if ( self->primme.applyPreconditioner != NULL ) self->primme.correctionParams.precondition = 1 ;
                primme_set_method ( Default_PrimmeMethod, &(self->primme) ) ;

                /* . Assignment. */
                self->eigenvalues   = target->eigenvalues  ;
                self->eigenvectors  = target->eigenvectors ;
                self->target        = target               ;

                /* . Memory required. */
                dprimme ( NULL, NULL, NULL, &(self->primme) ) ;

                /* . Allocate space. */
                self->residualNorms   = RealArray1D_AllocateWithExtent    ( n, status ) ;
                self->iWork           = IntegerArray1D_AllocateWithExtent ( self->primme.intWorkSize  / sizeof ( Integer ), status ) ;
                self->rWork           = RealArray1D_AllocateWithExtent    ( self->primme.realWorkSize / sizeof ( Real    ), status ) ;
                self->primme.intWork  = Array1D_Data ( self->iWork ) ;
                self->primme.realWork = Array1D_Data ( self->rWork ) ;
                if ( ( self->iWork == NULL ) || ( self->rWork == NULL ) || ( self->residualNorms == NULL ) ) JDEigenvalueSolverState_Deallocate ( &self ) ;
            }
            if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
        }
        else Status_Set ( status, Status_InvalidArgument ) ;
    }
    return self ;
}

/*==================================================================================================================================
! . Target.
!=================================================================================================================================*/
void JDEigenvalueSolverTarget_Initialize ( JDEigenvalueSolverTarget *self )
{
    if ( self != NULL )
    {
        self->eigenvalues         = NULL ;
        self->eigenvectors        = NULL ;
        self->object              = NULL ;
        self->ApplyMatrix         = NULL ;
        self->ApplyPreconditioner = NULL ;
    }
}
