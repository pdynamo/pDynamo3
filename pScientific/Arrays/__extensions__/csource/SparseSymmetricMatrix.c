/*==================================================================================================================================
! . Real sparse symmetric matrices.
!=================================================================================================================================*/

/*
! . There are many alternative ways of storing items in the matrix. Currently all diagonal items
! . are stored followed by all off-diagonal items. This is done because:
!
! * The separate storage of diagonal and off-diagonal items facilitates matrix-vector multiplication.
! * The storage of all diagonal items (as opposed to only non-zero ones) is required for some
!   operations (e.g. Cholesky decomposition) and also facilitates some others (e.g. row iteration).
!
! . Equally, at the moment, both diagonal and off-diagonal items employ the same item type although
! . this involves redundant extra storage for diagonal items.
!
! . A matrix is canonical when its off-diagonal items are put in ascending order of (i,j) pairs
! . with i > j.
!
! . These choices may change if matrices with many zero diagonal elements are to be routinely treated.
*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "BooleanArray1D.h"
# include "Memory.h"
# include "NumericalMacros.h"
# include "SparseSymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Integer Value_Compare ( const void *vItem1, const void *vItem2 ) ;

static void SparseSymmetricMatrix_IndexItems              ( SparseSymmetricMatrix *self ) ;
static void SparseSymmetricMatrix_InitializeDiagonalItems ( SparseSymmetricMatrix *self ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Minimum size allocation is extent. */
SparseSymmetricMatrix *SparseSymmetricMatrix_Allocate ( const Integer extent, const Integer size, Status *status )
{
    SparseSymmetricMatrix *self = NULL ;
    if ( extent < 0 ) Status_Set ( status, Status_InvalidArgument ) ;
    else
    {
        self = Memory_AllocateType ( SparseSymmetricMatrix ) ;
        if ( self != NULL )
        {
            /* . Scalars. */
            self->isCanonical            = False  ;
            self->extent                 = extent ;
            self->maximumNonZeroRowItems = 0      ;
            self->numberOfItems          = extent ;
            self->size                   = Maximum ( size, extent ) ;
            /* . Arrays. */
            self->rowIndex               = NULL ;
            self->items                  = NULL ;
            /* . Allocation. */
            self->items    = Memory_AllocateArrayOfTypes ( self->size, SparseSymmetricMatrixItem ) ;
            self->rowIndex = IntegerArray1D_AllocateWithExtent ( extent+1, status ) ;
            if ( ( self->items == NULL ) || ( self->rowIndex == NULL ) ) SparseSymmetricMatrix_Deallocate ( &self ) ;
            else
            {
                IntegerArray1D_Set ( self->rowIndex, -1 ) ;
                SparseSymmetricMatrix_InitializeDiagonalItems ( self ) ;
            }
        }
        if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    }
    return self ;
}


/*----------------------------------------------------------------------------------------------------------------------------------
! . Append an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Use only for off-diagonal items? */

void SparseSymmetricMatrix_AppendItem ( SparseSymmetricMatrix *self, const Integer i, const Integer j, const Real value, Status *status )
{
    if ( self != NULL )
    {
	auto Integer n ;
        n = self->numberOfItems ;
        /* . General checking. */
        if ( ( i < 0 ) || ( i >= self->extent ) || ( j < 0 ) || ( j >= self->extent ) ) Status_Set ( status, Status_IndexOutOfRange ) ;
        /* . Diagonal items. */
	else if ( i == j )
	{
	    if ( self->items[i].value == 0.0e+00 ) self->items[i].value = value ;
	    else Status_Set ( status, Status_InvalidArgument ) ; /* A[i,i] already exists */
	}
        /* . Overflow. */
        else if ( n >= self->size ) Status_Set ( status, Status_NonConformableArrays ) ;
	/* . Off-diagonal items. */
        else
        {
	    self->isCanonical    = False ;
            self->items[n].i     =     i ;
            self->items[n].j     =     j ;
            self->items[n].next  =    -1 ;
            self->items[n].value = value ;
            self->numberOfItems ++ ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Apply an incomplete Cholesky decomposition to a vector.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . The equation M * x = b is solved where M is an incomplete Cholesky decomposition
! . computed by ComputeIncompleteCholeskyDecomposition.
*/

void SparseSymmetricMatrix_ApplyIncompleteCholeskyDecomposition ( const SparseSymmetricMatrix *self, const RealArray1D *b, RealArray1D *x )
{
    if ( self != NULL )
    {
        auto Integer c, i, j, n ;
	auto Real    sum ;

        /* . Extent. */
        n = self->extent ;

        /* . Initialization. */
        RealArray1D_CopyTo ( b, x, NULL ) ;

        /* . Compute x = P^-1 * b = P^T * b. */
/*
        if ( self->permutation == NULL ) RealArray1D_CopyTo ( b, x, NULL ) ;
        else
        {
            for ( i = 0 ; i < n ; i++ ) Array1D_Item ( x, i ) = Array1D_Item ( b, Array1D_Item ( self->permutation, i ) ) ;
        }
*/

        /* . Compute x = L^-1 * x. */
        for ( i = 1 ; i < n ; i++ )
        {
	    sum = 0.0e+00 ;
            for ( j = Array1D_Item ( self->rowIndex, i ) ; j < Array1D_Item ( self->rowIndex, i+1 ) ; j++ )
            {
                c = self->items[j].j ;
                sum += ( self->items[j].value * Array1D_Item ( x, c ) ) ;
            }
	    Array1D_Item ( x, i ) -= sum ;
        }

        /* . x = D^-1 * x. */
        for ( i = 0 ; i < n ; i++ ) Array1D_Item ( x, i ) *= self->items[i].value ;

        /* . Compute x = (L^T)^-1 * x. */
        for ( i = n-1 ; i > 0 ; i-- )
        {
            for ( j = Array1D_Item ( self->rowIndex, i ) ; j < Array1D_Item ( self->rowIndex, i+1 ) ; j++ )
            {
                c = self->items[j].j ;
                Array1D_Item ( x, c ) -= ( self->items[j].value * Array1D_Item ( x, i ) ) ;
            }
        }

        /* . Compute x = (P^T)^-1 * b = P * x. */
/*
        RealArray1D_Permute ( x, self->permutation, NULL ) ;
*/
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Canonicalize the matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SparseSymmetricMatrix_Canonicalize ( SparseSymmetricMatrix *self, Status *status )
{
    if ( ( self != NULL ) && ( ! self->isCanonical ) )
    {
        /* . There are off-diagonal items. */
	if ( self->numberOfItems > self->extent )
	{
            auto Integer i, j, l, n ;

            /* . Ensure i > j for each off-diagonal item. */
            for ( l = self->extent ; l < self->numberOfItems ; l++ )
            {
        	i = self->items[l].i ;
        	j = self->items[l].j ;
        	if ( i < j ) { self->items[l].i = j ; self->items[l].j = i ; }
            }

            /* . Put all off-diagonal items in order. */
	    n = self->numberOfItems - self->extent ;
            qsort ( ( void * ) SparseSymmetricMatrix_ItemPointer ( self, self->extent ), ( size_t ) n, sizeof ( SparseSymmetricMatrixItem ), ( void * ) Value_Compare ) ;

            /* . Remove items with duplicate indices. */
            for ( l = n = self->extent+1 ; l < self->numberOfItems ; l++ )
            {
        	i = self->items[l].i ;
        	j = self->items[l].j ;
        	if ( ( i != self->items[n-1].i ) || ( j != self->items[n-1].j ) )
        	{
                    if ( n != l )
                    {
                	self->items[n].i     = self->items[l].i ;
                	self->items[n].j     = self->items[l].j ;
                	self->items[n].value = self->items[l].value ;
                    }
                    n++ ;
        	}
            }
            if ( n < self->numberOfItems ) Status_Set ( status, Status_AlgorithmError ) ;
            self->numberOfItems = n ;
        }

        /* . Index the items. */
	SparseSymmetricMatrix_IndexItems ( self ) ;

        /* . Finish up. */
        self->isCanonical = True ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Clear all items.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SparseSymmetricMatrix_Clear ( SparseSymmetricMatrix *self )
{
    if ( self != NULL )
    {
        self->isCanonical            = False  ;
        self->maximumNonZeroRowItems =     0  ;
        self->numberOfItems          = self->extent ;
        SparseSymmetricMatrix_InitializeDiagonalItems ( self ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
! . A minimum size array is allocated.
!---------------------------------------------------------------------------------------------------------------------------------*/
SparseSymmetricMatrix *SparseSymmetricMatrix_Clone ( const SparseSymmetricMatrix *self, Status *status )
{
    SparseSymmetricMatrix *clone = NULL ;
    if ( self != NULL )
    {
        clone = SparseSymmetricMatrix_Allocate ( self->extent, self->numberOfItems, status ) ;
        SparseSymmetricMatrix_CopyTo ( self, clone, status ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Compute a modified or unmodified incomplete Cholesky decomposition in-place.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . The decomposition that is computed is A = M + R where M = P*L*D*L^T*P^T.
! . Here L is lower triangular with unit diagonal, D is diagonal and P a permutation.
! . It is stored as (L-I) + D. D can also be inverted D^-1.
! . A non-zero alpha permits the modified ICD to be computed (0 <= alpha <= 1).
!
! . At the moment no fill-in is attempted so P is the identity.
!
! . For fill-in, Markowitz strategy chooses current row which has least number of non-zero elements.
! . When there are ties - choose row of lowest index.
*/

# define _SmallFactor 1.0e-02

void SparseSymmetricMatrix_ComputeIncompleteCholeskyDecomposition ( SparseSymmetricMatrix *self, const Real alpha, Integer *numberOfModifiedPivots, Status *status )
{
    if ( self != NULL )
    {
        auto Boolean                              doFillIn ;
        auto BooleanArray1D                      *visitedIndices ;
	auto Integer                              i, j, jActive, k, kActive, kIndex, l, lActive, lower, n, next, numberActive, pivotRow ;
        auto IntegerArray1D                      *columnIndices, *itemIndices ;
        auto Real                                 f, fillIn, pivot, reciprocalPivot, sum ;
        auto SparseSymmetricMatrixItem           *item ;
	auto SparseSymmetricMatrixRowItemIterator iterator ;

        /* . Initialization. */
	doFillIn = ( alpha != 0.0e+00 ) ;
        n        = self->extent ;
	(*numberOfModifiedPivots) = 0 ;
        SparseSymmetricMatrix_Canonicalize ( self, status ) ;

        /* . Allocate space. */
        columnIndices  = IntegerArray1D_AllocateWithExtent ( self->maximumNonZeroRowItems, status ) ;
        itemIndices    = IntegerArray1D_AllocateWithExtent ( self->maximumNonZeroRowItems, status ) ;
	visitedIndices = BooleanArray1D_AllocateWithExtent ( self->maximumNonZeroRowItems, status ) ;
/*
        iPivot        = IntegerArray1D_AllocateWithExtent ( n, status ) ;
        IntegerArray1D_Set ( iPivot, -1 ) ;
*/

        /* . Perform the decomposition. */
        if ( Status_IsOK ( status ) )
        {
	    /* . Loop over all rows. */
            for ( i = 0 ; i < n ; i++ )
            {

                /* . Select pivot row. */
	        /* . The Markowitz strategy selects the active row with the smallest
	        !  . number of non-zero elements. In the event of a tie, the row of
	        !  . lowest index wins.
	        */

                /* . Use the next row. */
                pivotRow = i ;

                /* . Save iPivot - the row is now inactive. */
/*
	        Array1D_Item ( iPivot, pivotRow ) = i ;
*/

                /* . Find the indices of the non-zero active elements in the row. */
	        /* . The column indices will always be in ascending order. */
	        numberActive = 0 ;
	        sum          = 0.0e+00 ;
                SparseSymmetricMatrixRowItemIterator_Initialize ( &iterator, self, pivotRow, status ) ;
	        while ( ( next = SparseSymmetricMatrixRowItemIterator_Next ( &iterator ) ) >= 0 )
	        {
	            item = &(self->items[next]) ;
		    if ( item->i == pivotRow ) j = item->j ;
		    else                       j = item->i ;
/*
		    if ( Array1D_Item ( iPivot, j ) < 0 )
*/
		    if ( j > i )
		    {
		        Array1D_Item ( columnIndices , numberActive ) = j     ;
		        Array1D_Item ( itemIndices   , numberActive ) = next  ;
		        Array1D_Item ( visitedIndices, numberActive ) = False ;
   	                numberActive += 1 ;
		        sum          += fabs ( item->value ) ;
		    }
	        }

                /* . Find the reciprocal of the diagonal element. */
	        /* . A check is made to ensure positive-definiteness. */
                pivot = self->items[pivotRow].value ;
                if ( pivot <= _SmallFactor * sum )
	        {
                    self->items[pivotRow].value = sum ;
                    if ( sum == 0.0e+00 ) self->items[pivotRow].value = 1.0e+00 ;
                    (*numberOfModifiedPivots) += 1 ;
	        }
                reciprocalPivot = 1.0e+00 / self->items[pivotRow].value ;
/*
printf ( "\nRow Number Active: %d %d %25.15f %25.15f\n", i, numberActive, sum, reciprocalPivot ) ;
IntegerArray1D_Print ( columnIndices ) ;
IntegerArray1D_Print ( itemIndices   ) ;
*/
                /* . Loop over non-zero active rows. */
                for ( jActive = 0 ; jActive < numberActive ; jActive++ )
	        {

                    /* . Initialization - Aij/Aii. */
                    j = Array1D_Item ( columnIndices, jActive ) ;
                    f = reciprocalPivot * self->items[Array1D_Item ( itemIndices, jActive )].value ;
/*
printf ( "\nJ Active>: %d %d %25.15f\n", jActive, j, f ) ;
*/
                    /* . Fill and Markowitz manipulations here. */

                    /* . Loop over non-zero active items of row j - including the diagonal. */
		    /* . k always increases. */
                    lower = 0 ;
                    SparseSymmetricMatrixRowItemIterator_Initialize ( &iterator, self, j, status ) ;
	            while ( ( next = SparseSymmetricMatrixRowItemIterator_Next ( &iterator ) ) >= 0 )
		    {
	                item = &(self->items[next]) ;
		        if ( item->i == j ) k = item->j ;
		        else                k = item->i ;
/*
printf ( "\nk = %d\n", k ) ;
*/
		        /* . Search to see if k also occurs in the active rows of i. */
		        kActive = -1 ;
		        for ( lActive = lower ; lActive < numberActive ; lActive++ )
		        {
			    l = Array1D_Item ( columnIndices, lActive ) ;
		                 if ( l >  k ) break ;
			    else if ( l == k ) { kActive = lActive ; lower = lActive + 1 ; break ; }
			    else lower = lActive ; /* . l < k. */
		        }
		        /*
		        ! . If the column is active and also occurs in the active columns of i,
		        ! . modify Ajk by - Aji*Aki/Aii and flag the item as having been used.
		        */
/*
		        if ( ( Array1D_Item ( iPivot, k ) < 0 ) && ( kActive >= 0 ) )
*/
		        if ( ( k > i ) && ( kActive >= 0 ) )
		        {
		            kIndex      = Array1D_Item ( itemIndices, kActive ) ;
/*printf ( "\nMatch k = %d %d %d %25.15f\n", k, kActive, kIndex, - ( f * self->items[kIndex].value ) ) ;*/
			    item->value -= ( f * self->items[kIndex].value ) ;
                            if ( doFillIn ) Array1D_Item ( visitedIndices, kActive ) = True ;
		        }
		    }

                    /* . Check for fill-in or a modified ICD. */
                    if ( doFillIn )
		    {
                        /* . Reloop over active items. */
	                for ( kActive = 0 ; kActive < numberActive ; kActive++ )
		        {
                	    /* . Reactivate already visited items. */
			    if ( Array1D_Item ( visitedIndices, kActive ) )
			    {
			        Array1D_Item ( visitedIndices, kActive ) = False ;
			    }
			    /* . Item not visited. */
			    /*
			    ! . This means that the pivot row has a non-zero entry but row j does not.
			    ! . This will create fill-in or a correction to the corresponding diagonal
			    ! . items.
			    */
        		    else
			    {
		                /* . Calculate fill-in if appropriate. */
                	        k = Array1D_Item ( columnIndices, kActive ) ;
                                if ( j >= k )
			        {
                        	    fillIn = -f * self->items[Array1D_Item ( itemIndices, kActive )].value ; /* . - Aji*Aki/Aii. */

                        	    /* . Fill-in manipulations - add into Ajk. */

                        	    /* . Add in the correction to the diagonals (modified ICD). */
                        	    fillIn *= alpha ;
                        	    self->items[j].value += fillIn ;
                        	    self->items[k].value += fillIn ;
                                }
        		    }
                        }
		    }
	        }
            }

            /* . Replace diagonal elements by their inverses. */
	    for ( i = 0 ; i < n ; i++ )
	    {
	        self->items[i].value = 1.0e+00 / self->items[i].value ;
	    }

            /* . Scale all lower triangular elements in a column by the column diagonal. */
	    for ( i = n ; i < self->numberOfItems ; i++ )
	    {
	        item = &(self->items[i]) ;
	        item->value *= self->items[item->j].value ;
	    }

            /* . Invert iPivot. */
/*
	    IntegerArray1D_CopyTo ( iPivot, iWork, NULL ) ;
	    for ( i = 0 ; i < n ; i++ ) Array1D_Item ( iPivot, Array1D_Item ( iWork, i ) ) = i ;
*/
        }

        /* . Finish up - recanonicalization unnecesary in principle. */
        IntegerArray1D_Deallocate ( &columnIndices  ) ;
        IntegerArray1D_Deallocate ( &itemIndices    ) ;
        BooleanArray1D_Deallocate ( &visitedIndices ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copying.
! . As much off-diagonal data as possible is copied.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SparseSymmetricMatrix_CopyTo ( const SparseSymmetricMatrix *self, SparseSymmetricMatrix *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        /* . Matrix dimension. */
        if ( self->extent == other->extent )
        {
            auto Integer l, n ;
            if ( self->numberOfItems <= other->size ) n = self->numberOfItems ;
            else
            {
                n = other->size ;
                Status_Set ( status, Status_NonConformableArrays ) ;
            }
            for ( l = 0 ; l < n ; l++ )
            {
                other->items[l].i     = self->items[l].i     ;
                other->items[l].j     = self->items[l].j     ;
                other->items[l].value = self->items[l].value ;
            }
            other->numberOfItems = n ;
            if ( self->isCanonical ) SparseSymmetricMatrix_Canonicalize ( other, status ) ;
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SparseSymmetricMatrix_Deallocate ( SparseSymmetricMatrix **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        IntegerArray1D_Deallocate ( &((*self)->rowIndex) ) ;
        Memory_Deallocate         (   (*self)->items     ) ;
        Memory_Deallocate         (   (*self)            ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Retrieve diagonal elements.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SparseSymmetricMatrix_GetDiagonal ( const SparseSymmetricMatrix *self, RealArray1D *diagonal, Status *status )
{
    RealArray1D_Set ( diagonal, 0.0e+00 ) ;
    if ( ( self != NULL ) && ( diagonal != NULL ) )
    {
        if ( self->extent == diagonal->extent )
        {
            auto Integer l ;
            for ( l = 0 ; l < self->extent ; l++ ) Array1D_Item ( diagonal, l ) = self->items[l].value ;
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Index the items.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void SparseSymmetricMatrix_IndexItems ( SparseSymmetricMatrix *self )
{
    if ( self != NULL )
    {
        auto Integer l, last, n ;

        /* . Find the maximum number of non-zero items in a row. */
        /* . Initialization - use rowIndex as workspace. */
        IntegerArray1D_Set ( self->rowIndex, 1 ) ;
        for ( l = self->extent ; l < self->numberOfItems ; l++ )
        {
            Array1D_Item ( self->rowIndex, self->items[l].i ) += 1 ;
            Array1D_Item ( self->rowIndex, self->items[l].j ) += 1 ;
        }
        self->maximumNonZeroRowItems = IntegerArray1D_Maximum ( self->rowIndex ) ;

        /* . Set up the row index - do in reverse order. */
        /* . Initialization. */
        IntegerArray1D_Set ( self->rowIndex, -1 ) ;

        /* . Upper-triangular items. */
        for ( l = self->numberOfItems - 1 ; l >= self->extent ; l-- )
        {
            last = Array1D_Item ( self->rowIndex, self->items[l].j ) ;
            self->items[l].next = last ;
            Array1D_Item ( self->rowIndex, self->items[l].j ) = l ;
        }

        /* . Diagonal items. */
        for ( l = 0 ; l < self->extent ; l++ )
        {
            last = Array1D_Item ( self->rowIndex, self->items[l].i ) ;
            self->items[l].next = last ;
            Array1D_Item ( self->rowIndex, self->items[l].i ) = l ;
        }

        /* . Row index. */
        IntegerArray1D_Set ( self->rowIndex, 0 ) ;
        for ( l = self->extent ; l < self->numberOfItems ; l++ )
        {
            Array1D_Item ( self->rowIndex, self->items[l].i ) += 1 ;
        }
	last = self->extent ;
        for ( l = 0 ; l < self->extent ; l++ )
        {
	    n = Array1D_Item ( self->rowIndex, l ) ;
            Array1D_Item ( self->rowIndex, l ) = last ;
	    last += n ;
        }
        Array1D_Item ( self->rowIndex, self->extent ) = last ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialize the diagonal items.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void SparseSymmetricMatrix_InitializeDiagonalItems ( SparseSymmetricMatrix *self )
{
    if ( ( self != NULL ) && ( self->items != NULL ) )
    {
        auto Integer l ;
        for ( l = 0 ; l < self->extent ; l++ )
        {
            self->items[l].i     =  l ;
            self->items[l].j     =  l ;
            self->items[l].next  = -1 ;
            self->items[l].value = 0.0e+00 ;
        }
    }
}


/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a diagonal preconditioner from the matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _Small 1.0e-06
void SparseSymmetricMatrix_MakeDiagonalPreconditioner ( const SparseSymmetricMatrix *self, RealArray1D *preconditioner, const Real *tolerance, Status *status )
{
    if ( ( self != NULL ) && ( preconditioner != NULL ) )
    {
        auto Integer i ;
        auto Real    t, v ;
        if ( tolerance == NULL ) t = _Small ;
        else                     t = (*tolerance ) ;
        SparseSymmetricMatrix_GetDiagonal ( self, preconditioner, status ) ;
        for ( i = 0 ; i < preconditioner->extent ; i++ )
        {
            v = Maximum ( fabs ( Array1D_Item ( preconditioner, i ) ), t ) ;
            Array1D_Item ( preconditioner, i ) = 1.0e+00 / v ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Printing.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SparseSymmetricMatrix_Print ( const SparseSymmetricMatrix  *self )
{
    if ( self == NULL ) printf ( "\nNull sparse symmetric matrix.\n" ) ;
    else
    {
        auto Integer l, start ;

        /* . Items. */
        printf ( "\nItems (index, i, j, next-in-column, value):\n" ) ;
        for ( l = 0 ; l < self->numberOfItems ; l++ ) printf ( "%10d %10d %10d %10d %15.10f\n", l, self->items[l].i, self->items[l].j, self->items[l].next, self->items[l].value ) ;

        /* . Row index. */
        if ( self->isCanonical )
        {
            printf ( "\nRow Index (row, start, number of items):\n" ) ;
            for ( l = 0 ; l < self->extent ; l++ )
            {
                start = Array1D_Item ( self->rowIndex, l ) ;
                printf ( "%10d %10d %10d\n", l, start, Array1D_Item ( self->rowIndex, l+1 ) - start ) ;
            }
        }

        /* . Other data. */
        printf ( "\nOther Data (extent, size, numberOfItems, maximumNonZeroRowItems): %10d %10d %10d %10d\n", self->extent, self->size, self->numberOfItems, self->maximumNonZeroRowItems ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Matrix-vector multiplication.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SparseSymmetricMatrix_VectorMultiply ( const SparseSymmetricMatrix *self, const RealArray1D *x, RealArray1D *y, Status *status )
{
    RealArray1D_Set ( y, 0.0e+00 ) ;
    if ( ( self != NULL ) && ( x != NULL ) && ( y != NULL ) )
    {
        if ( ( self->extent == x->extent ) && ( self->extent == y->extent ) )
        {
            auto Integer i, j, l ;
            auto Real    value ;
            for ( l = 0 ; l < self->extent ; l++ ) Array1D_Item ( y, l ) += self->items[l].value * Array1D_Item ( x, l ) ;
            for ( l = self->extent ; l < self->numberOfItems ; l++ )
            {
                i     = self->items[l].i ;
                j     = self->items[l].j ;
                value = self->items[l].value ;
                Array1D_Item ( y, i ) += value * Array1D_Item ( x, j ) ;
                Array1D_Item ( y, j ) += value * Array1D_Item ( x, i ) ;
            }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*==================================================================================================================================
! . Row iterator procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SparseSymmetricMatrixRowItemIterator_Initialize ( SparseSymmetricMatrixRowItemIterator *self, SparseSymmetricMatrix *target, const Integer row, Status *status )
{
    if ( ( self != NULL ) && ( target != NULL ) && ( target->isCanonical ) )
    {
        if ( ( row >= 0 ) && ( row < target->extent ) )
        {
            auto Integer n ;
            /* . Initialization. */
            self->inLT    = True   ;
	    self->current = -1     ;
            self->ltLast  = Array1D_Item ( target->rowIndex, row + 1 ) ;
            self->row     = row    ;
            self->target  = target ;
	    /* . Find the first item in row, either in lower triangle or the diagonal. */
	    n = self->ltLast - Array1D_Item ( target->rowIndex, row ) ;
	    if ( n > 0 ) self->current = Array1D_Item ( target->rowIndex, row ) ;
	    else       { self->current = row ; self->inLT = False ; }
        }
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Next item.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer SparseSymmetricMatrixRowItemIterator_Next ( SparseSymmetricMatrixRowItemIterator *self )
{
    Integer current = -1 ;
    if ( ( self != NULL ) && ( self->current != -1 ) )
    {
        /* . Set current. */
        current = self->current ;
	/* . Find next item. */
        if ( self->inLT )
        {
            /* . Move to diagonal. */
	    if ( current == ( self->ltLast - 1 ) ) { self->current = self->row ; self->inLT = False ; }
            /* . Continue in lower triangle. */
	    else self->current += 1 ;
        }
        /* . Upper triangle. */
        else self->current = self->target->items[current].next ;
/*
if ( ( current >= self->target->numberOfItems ) || ( self->current >= self->target->numberOfItems ) )
{
printf ( "\nInformation: %d %d %d %d %d %d %d\n", current, self->current, self->target->extent, self->inLT, self->ltLast, self->row, self->target->numberOfItems ) ;
}
*/
    }
    return current ;
}

/*==================================================================================================================================
! . Private procedures.
!=================================================================================================================================*/
/* . This is applied to off-diagonal items only. */
/* . For this to work, i >= j. */
static Integer Value_Compare ( const void *vItem1, const void *vItem2 )
{
    Integer i ;
    SparseSymmetricMatrixItem *item1, *item2 ;
    item1 = ( SparseSymmetricMatrixItem * ) vItem1 ;
    item2 = ( SparseSymmetricMatrixItem * ) vItem2 ;
         if ( item1->i < item2->i ) i = -1 ;
    else if ( item1->i > item2->i ) i =  1 ;
    else if ( item1->j < item2->j ) i = -1 ;
    else if ( item1->j > item2->j ) i =  1 ;
    else i = 0 ;
    return i ;
}
