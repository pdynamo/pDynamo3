/*==================================================================================================================================
! . Dense determinants.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "f2clapack.h"

# include "BooleanArray1D.h"
# include "DenseDeterminants.h"
# include "Integer.h"
# include "Memory.h"
# include "NumericalMacros.h"
# include "Real.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . LUP factorization of a matrix in-place.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . self is overwritten with L in the lower triangle and U in the upper triangle with diagonal.
! . The diagonal elements of L are 1.
*/
void Matrix_LUPFactorizationInPlace ( RealArray2D *self, IntegerArray1D *pivots, Status *status )
{
    if ( ( self != NULL ) && ( pivots != NULL ) )
    {
        auto Integer info, m, n ;

        /* . Array dimensions. */
        info = 0 ;
        m    = View2D_Rows    ( self ) ;
        n    = View2D_Columns ( self ) ;
        if ( ( View1D_Extent ( pivots ) >= Minimum ( m, n ) ) && ( self->stride1 == 1 ) )
        {
            IntegerArray1D_Set ( pivots, -1 ) ;
            /* . Reverse m and n due to C <-> Fortran ordering. */
            dgetrf_( &n, &m, Array2D_Data ( self ), &(self->stride0), Array1D_Data ( pivots ), &info ) ;
        }
        else Status_Set ( status, Status_InvalidArgument ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Determine the determinant of a square matrix using LUP factorization.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . self is destroyed as it is overwritten with the LUP factorization.
*/
Real SquareMatrix_Determinant ( RealArray2D *self, Status *status )
{
    Real determinant = 0.0e+00 ;
    if ( self != NULL )
    {
        /* . The matrix must be square. */
        if ( View2D_Rows ( self ) == View2D_Columns ( self ) )
        {
            auto Integer         current, cycleLength, i, n, p, start, t ;
            auto Real            parity ;
            auto BooleanArray1D *isVisited ;
            auto IntegerArray1D *permutation, *pivots ;

            /* . Allocate space. */
            n = View2D_Rows ( self ) ;
            isVisited   = BooleanArray1D_AllocateWithExtent ( n, status ) ;
            permutation = IntegerArray1D_AllocateWithExtent ( n, status ) ;
            pivots      = IntegerArray1D_AllocateWithExtent ( n, status ) ;

            /* . Do the factorization. */
            Matrix_LUPFactorizationInPlace ( self, pivots, status ) ;

            /* . Find the determinant of U. */
            determinant = 1.0e+00 ;
            for ( i = 0 ; i < n ; i++ ) determinant *= Array2D_Item ( self, i, i ) ;

            /* . Convert pivots to a permutation. */
            for ( i = 0 ; i < n   ; i++ ) Array1D_Item ( permutation, i ) = i ;
            for ( i = 0 ; i < n-1 ; i++ )
            {
                p = Array1D_Item ( pivots,      i ) - 1 ;
                t = Array1D_Item ( permutation, p ) ;
                Array1D_Item ( permutation, p ) = Array1D_Item ( permutation, i ) ;
                Array1D_Item ( permutation, i ) = t ;
            }

            /* . Find the parity of the pivots. */
            parity = 1.0e+00 ;
            BooleanArray1D_Set ( isVisited, False ) ;
            for ( i = 0 ; i < n ; i++ )
            {
                if ( ! Array1D_Item ( isVisited, i ) )
                {
                    current     = Array1D_Item ( permutation, i ) ;
                    cycleLength = 1 ;
                    start       = i ;
                    Array1D_Item ( isVisited, i ) = True ;
                    while ( current != start )
                    {
                        Array1D_Item ( isVisited, current ) = True ;
                        current      = Array1D_Item ( permutation, current ) ;
                        cycleLength += 1 ;
                    }
                    if ( cycleLength % 2 == 0 ) parity *= -1.0e+00 ;
                }
            }
            determinant *= parity ;

            /* . Finish up. */
            BooleanArray1D_Deallocate ( &isVisited   ) ;
            IntegerArray1D_Deallocate ( &permutation ) ;
            IntegerArray1D_Deallocate ( &pivots      ) ;
        }
        else Status_Set ( status, Status_InvalidArgument ) ;
    }
    return determinant ;
}
