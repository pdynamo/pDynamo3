/*==================================================================================================================================
! . Module for CI sparse matrix diagonalization.
! . Functions for the primme solver.
!=================================================================================================================================*/

# include <stdio.h>
# include <stdlib.h>

# include "CISparseSolver.h"
# include "Real.h"
# include "RealArray1D.h"
# include "SparseSymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Apply the CI matrix to a vector.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CISparseSolver_ApplyMatrix ( void *xVoid, void *yVoid, Integer *blockSize, primme_params *primme )
{
    Real                  *x, *y ;
    RealArray1D            hV, v ;
    SparseSymmetricMatrix *matrix ;
    /* . Casts. */
    matrix = ( SparseSymmetricMatrix * ) primme->matrix ;
    x      = ( Real * ) xVoid ;
    y      = ( Real * ) yVoid ;
    /* . Views. */
    RealArray1D_ViewOfRaw (  &v, 0, primme->n, 1, x ) ;
    RealArray1D_ViewOfRaw ( &hV, 0, primme->n, 1, y ) ;
    /* . Calculate H * v. */
    SparseSymmetricMatrix_VectorMultiply ( matrix, &v, &hV, NULL ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Apply the CI matrix preconditioner to a vector.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CISparseSolver_ApplyPreconditioner ( void *xVoid, void *yVoid, Integer *blockSize, primme_params *primme )
{
    Real        *x, *y ;
    RealArray1D  hV, v, *preconditioner ;
    /* . Casts. */
    preconditioner = ( RealArray1D * ) primme->preconditioner ;
    x              = ( Real * ) xVoid ;
    y              = ( Real * ) yVoid ;
    /* . Views. */
    RealArray1D_ViewOfRaw (  &v, 0, primme->n, 1, x ) ;
    RealArray1D_ViewOfRaw ( &hV, 0, primme->n, 1, y ) ;
    /* . Calculate v / Diagonal ( h ). */
    RealArray1D_CopyTo   ( &v , &hV           , NULL ) ;
    RealArray1D_Multiply ( &hV, preconditioner, NULL ) ;
}
