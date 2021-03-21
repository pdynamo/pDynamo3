/*==================================================================================================================================
! . Real antisymmetric matrices.
!=================================================================================================================================*/

# include <stdio.h>
# include <stdlib.h>

# include "AntisymmetricMatrix.h"
# include "Array_Macros.h"
# include "Iterator1D.h"
# include "Memory.h"
# include "NumericalMacros.h"

# define _UseCBLAS
# define _UseReal
# include "Utilities1D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . The abolute maximum value.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real AntisymmetricMatrix_AbsoluteMaximum ( const AntisymmetricMatrix *self )
{
    Real value = 0.0e+00 ;
    if ( self != NULL )
    {
        Utilities1D_AbsoluteMaximum ( self->size, self->data, 1, value ) ;
    }
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Add a scaled matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_Add ( AntisymmetricMatrix *self, const Real alpha, const AntisymmetricMatrix *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL )  && Status_IsOK ( status ) )
    {
        Utilities1D_Add ( self->size, self->data, 1, other->size, other->data, 1, alpha, status ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Raw allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
AntisymmetricMatrix *AntisymmetricMatrix_Allocate ( Status *status )
{
    AntisymmetricMatrix *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( AntisymmetricMatrix ) ;
        if ( self != NULL )
        {
            self->block  = NULL ;
            self->data   = NULL ;
            self->extent = 0 ;
            self->size   = 0 ;
        }
        else Status_Set ( status, Status_OutOfMemory ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation with extent.
!---------------------------------------------------------------------------------------------------------------------------------*/
AntisymmetricMatrix *AntisymmetricMatrix_AllocateWithExtent ( const Integer extent, Status *status )
{
    AntisymmetricMatrix *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        auto RealBlock *block = NULL ;
        block = RealBlock_Allocate ( AntisymmetricMatrix_ViewSize ( extent ), status ) ;
        self  = AntisymmetricMatrix_FromExtentBlock ( extent, block, True   , status ) ;
        if ( self == NULL ) RealBlock_Deallocate ( &block ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Anticommutator of antisymmetric and symmetric matrices A * S + S * A.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_AnticommutatorAS ( const AntisymmetricMatrix *self, const SymmetricMatrix *a, AntisymmetricMatrix *result, Status *status )
{
    if ( ( self != NULL ) && ( a != NULL ) && ( result != NULL ) && Status_IsOK ( status ) )
    {
        if ( ( self->extent == a->extent ) && ( a->extent == result->extent ) )
        {
            auto Integer  i, j, k ;
            auto Real     sum ;
            for ( i = 1 ; i < self->extent ; i++ )
            {
                for ( k = 0 ; k < i ; k++ )
                {
                    sum = 0.0e+00 ;
                    /* . Care needed with signs when swapping indices. */
                    for ( j = 0         ; j < k            ; j++ ) sum += (   AntisymmetricMatrix_Item ( self, i, j ) * SymmetricMatrix_Item ( a, k, j ) - SymmetricMatrix_Item ( a, i, j ) * AntisymmetricMatrix_Item ( self, k, j ) ) ;
                    for ( j = ( k + 1 ) ; j < i            ; j++ ) sum += (   AntisymmetricMatrix_Item ( self, i, j ) * SymmetricMatrix_Item ( a, j, k ) + SymmetricMatrix_Item ( a, i, j ) * AntisymmetricMatrix_Item ( self, j, k ) ) ;
                    for ( j = ( i + 1 ) ; j < self->extent ; j++ ) sum += ( - AntisymmetricMatrix_Item ( self, j, i ) * SymmetricMatrix_Item ( a, j, k ) + SymmetricMatrix_Item ( a, j, i ) * AntisymmetricMatrix_Item ( self, j, k ) ) ;
                    sum += AntisymmetricMatrix_Item ( self, i, k ) * ( SymmetricMatrix_Item ( a, i, i ) + SymmetricMatrix_Item ( a, k, k ) ) ;
                    AntisymmetricMatrix_Item ( result, i, k ) = sum ;
                }
            }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deep cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
AntisymmetricMatrix *AntisymmetricMatrix_CloneDeep ( const AntisymmetricMatrix *self, Status *status )
{
    AntisymmetricMatrix *clone = NULL ;
    if ( self != NULL )
    {
        clone = AntisymmetricMatrix_AllocateWithExtent ( self->extent, status ) ;
        AntisymmetricMatrix_CopyTo ( self, clone, status ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Shallow cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
AntisymmetricMatrix *AntisymmetricMatrix_CloneShallow ( const AntisymmetricMatrix *self, Status *status )
{
    AntisymmetricMatrix *clone = NULL ;
    if ( self != NULL ) clone = AntisymmetricMatrix_FromExtentBlock ( self->extent, self->block, True, status ) ;
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Commutator of antisymmetric and symmetric matrices A * S - S * A.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_CommutatorAS ( const AntisymmetricMatrix *self, const SymmetricMatrix *a, SymmetricMatrix *result, Status *status )
{
    if ( ( self != NULL ) && ( a != NULL ) && ( result != NULL ) && Status_IsOK ( status ) )
    {
        if ( ( self->extent == a->extent ) && ( a->extent == result->extent ) )
        {
            auto Integer  i, j, k ;
            auto Real     sum ;
            for ( i = 0 ; i < self->extent ; i++ )
            {
                for ( k = 0 ; k < i ; k++ )
                {
                    sum = 0.0e+00 ;
                    for ( j = 0         ; j < k            ; j++ ) sum += (   AntisymmetricMatrix_Item ( self, i, j ) * SymmetricMatrix_Item ( a, k, j ) + SymmetricMatrix_Item ( a, i, j ) * AntisymmetricMatrix_Item ( self, k, j ) ) ;
                    for ( j = ( k + 1 ) ; j < i            ; j++ ) sum += (   AntisymmetricMatrix_Item ( self, i, j ) * SymmetricMatrix_Item ( a, j, k ) - SymmetricMatrix_Item ( a, i, j ) * AntisymmetricMatrix_Item ( self, j, k ) ) ;
                    for ( j = ( i + 1 ) ; j < self->extent ; j++ ) sum += ( - AntisymmetricMatrix_Item ( self, j, i ) * SymmetricMatrix_Item ( a, j, k ) - SymmetricMatrix_Item ( a, j, i ) * AntisymmetricMatrix_Item ( self, j, k ) ) ;
                    sum += AntisymmetricMatrix_Item ( self, i, k ) * ( SymmetricMatrix_Item ( a, k, k ) - SymmetricMatrix_Item ( a, i, i ) );
                    SymmetricMatrix_Item ( result, i, k ) = sum ;
                }
                sum = 0.0e+00 ;
                for ( j = 0         ; j < i            ; j++ ) sum += (   AntisymmetricMatrix_Item ( self, i, j ) * SymmetricMatrix_Item ( a, i, j ) + SymmetricMatrix_Item ( a, i, j ) * AntisymmetricMatrix_Item ( self, i, j ) ) ;
                for ( j = ( i + 1 ) ; j < self->extent ; j++ ) sum += ( - AntisymmetricMatrix_Item ( self, j, i ) * SymmetricMatrix_Item ( a, j, i ) - SymmetricMatrix_Item ( a, j, i ) * AntisymmetricMatrix_Item ( self, j, i ) ) ;
                SymmetricMatrix_Item ( result, i, i ) = sum ;
            }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Commutator of two symmetric matrices A * B - B * A. Faster version that requires much more memory.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_CommutatorSS_Fast (       AntisymmetricMatrix *self   ,
                                             const SymmetricMatrix     *a      ,
                                             const SymmetricMatrix     *b      ,
                                                   RealArray2D         *mA     ,
                                                   RealArray2D         *mB     ,
                                                   RealArray2D         *mC     ,
                                                   Status              *status )
{
    if ( ( self != NULL ) && ( a != NULL ) && ( b != NULL ) && ( mA != NULL ) && ( mB != NULL ) && ( mC != NULL ) && Status_IsOK ( status ) )
    {
        if ( ( self->extent == a->extent ) && ( a->extent == b->extent ) )
        {
            SymmetricMatrix_CopyToRealArray2D ( a, mA, status ) ;
            SymmetricMatrix_CopyToRealArray2D ( b, mB, status ) ;
            RealArray2D_MatrixMultiply ( False, False, 1.0e+00, mA, mB, 0.0e+00, mC, status ) ;
            AntisymmetricMatrix_CopyFromRealArray2D ( self, mC, False, status ) ;
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Commutator of two symmetric matrices A * B - B * A. Reference (slow) version.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_CommutatorSS_Reference ( AntisymmetricMatrix *self, const SymmetricMatrix *a, const SymmetricMatrix *b, Status *status )
{
    if ( ( self != NULL ) && ( a != NULL ) && ( b != NULL ) && Status_IsOK ( status ) )
    {
# ifdef USEOPENMP
	# pragma omp parallel /*num_threads(MAXIMUMNUMBEROFTHREADS)*/
# endif
        if ( ( self->extent == a->extent ) && ( a->extent == b->extent ) )
        {
            auto Integer  i, j, k ;
            auto Real     sum ;
# ifdef USEOPENMP
            # pragma omp for schedule(dynamic)
# endif
            for ( i = 0 ; i < self->extent ; i++ )
            {
                for ( k = 0 ; k < i ; k++ )
                {
                    sum = 0.0e+00 ;
                    for ( j = 0         ; j <= k           ; j++ ) sum += ( SymmetricMatrix_Item ( a, i, j ) * SymmetricMatrix_Item ( b, k, j ) - SymmetricMatrix_Item ( b, i, j ) * SymmetricMatrix_Item ( a, k, j ) ) ;
                    for ( j = ( k + 1 ) ; j <= i           ; j++ ) sum += ( SymmetricMatrix_Item ( a, i, j ) * SymmetricMatrix_Item ( b, j, k ) - SymmetricMatrix_Item ( b, i, j ) * SymmetricMatrix_Item ( a, j, k ) ) ;
                    for ( j = ( i + 1 ) ; j < self->extent ; j++ ) sum += ( SymmetricMatrix_Item ( a, j, i ) * SymmetricMatrix_Item ( b, j, k ) - SymmetricMatrix_Item ( b, j, i ) * SymmetricMatrix_Item ( a, j, k ) ) ;
                    AntisymmetricMatrix_Item ( self, i, k ) = sum ;
                }
            }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Commutator of three symmetric matrices A * B * C - C * B * A.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_CommutatorSSS ( AntisymmetricMatrix *self, const SymmetricMatrix *a, const SymmetricMatrix *b, const SymmetricMatrix *c, Status *status )
{
    if ( ( self != NULL ) && ( a != NULL ) && ( b != NULL ) && ( c != NULL ) && Status_IsOK ( status ) )
    {
        if ( ( self->extent == a->extent ) && ( a->extent == b->extent ) && ( b->extent == c->extent ) )
        {
            auto Integer  i, ij, j, jk, k, kl, l ;
            auto Real     sum1, sum2, *tab = NULL, *tcb = NULL ;
            tab = Memory_AllocateArrayOfTypes ( self->extent, Real ) ;
            tcb = Memory_AllocateArrayOfTypes ( self->extent, Real ) ;
            if ( ( tab != NULL ) && ( tcb != NULL ) )
            {
                for ( i = 0 ; i < self->extent ; i++ )
                {
                    for ( k = 0 ; k < self->extent ; k++ )
                    {
                        sum1 = 0.0e+00 ;
                        sum2 = 0.0e+00 ;
                        for ( j = 0 ; j < self->extent ; j++ )
                        {
                            if ( i >= j ) ij = ( i * ( i + 1 ) ) / 2 + j ;
                            else          ij = ( j * ( j + 1 ) ) / 2 + i ;
                            if ( j >= k ) jk = ( j * ( j + 1 ) ) / 2 + k ;
                            else          jk = ( k * ( k + 1 ) ) / 2 + j ;
                            sum1 += a->data[ij] * b->data[jk] ;
                            sum2 += c->data[ij] * b->data[jk] ;
                        }
                        tab[k] = sum1 ;
                        tcb[k] = sum2 ;
                    }
                    for ( l = 0 ; l < i ; l++ )
                    {
                        sum1 = 0.0e+00 ;
                        sum2 = 0.0e+00 ;
                        for ( k = 0 ; k < self->extent ; k++ )
                        {
                            if ( k >= l ) kl = ( k * ( k + 1 ) ) / 2 + l ;
                            else          kl = ( l * ( l + 1 ) ) / 2 + k ;
                            sum1 += tab[k] * c->data[kl] ;
                            sum2 += tcb[k] * a->data[kl] ;
                        }
                        AntisymmetricMatrix_Item ( self, i, l ) = sum1 - sum2 ;
                    }
                }
            }
            else Status_Set ( status, Status_OutOfMemory ) ;
            Memory_Deallocate ( tab ) ;
            Memory_Deallocate ( tcb ) ;
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Commutator of three symmetric matrices and a transformation M^T * ( A * B * C - C * B * A ) * M.
! . The transformation matrix can be transposed, i.e. M * ( A * B * C - C * B * A ) * M^T.
! . Fast version that uses a lot of memory.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_CommutatorTSSST (       AntisymmetricMatrix *self       ,
                                           const SymmetricMatrix     *a          ,
                                           const SymmetricMatrix     *b          ,
                                           const SymmetricMatrix     *c          ,
                                           const RealArray2D         *m          ,
                                           const Boolean              mTranspose ,
                                                 RealArray2D         *u          ,
                                                 RealArray2D         *v          ,
                                                 RealArray2D         *w          ,
                                                 Status              *status     )
{
    if ( ( self != NULL ) && ( a != NULL ) && ( b != NULL ) && ( c != NULL ) && ( m != NULL ) && ( u != NULL ) && ( v != NULL ) && ( w != NULL ) && Status_IsOK ( status ) )
    {
        auto Boolean isOK ;
        isOK = ( a->extent == b->extent ) && ( a->extent == c->extent ) ;
        if ( mTranspose ) isOK = isOK && ( ( c->extent == View2D_Columns ( m ) ) && ( View2D_Rows    ( m ) == self->extent ) ) ;
        else              isOK = isOK && ( ( c->extent == View2D_Rows    ( m ) ) && ( View2D_Columns ( m ) == self->extent ) ) ;
        if ( isOK )
        {
            SymmetricMatrix_CopyToRealArray2D ( a, u, status ) ;
            SymmetricMatrix_CopyToRealArray2D ( b, v, status ) ;
            RealArray2D_MatrixMultiply ( False, False, 1.0e+00, u, v, 0.0e+00, w, status ) ;
            SymmetricMatrix_CopyToRealArray2D ( c, u, status ) ;
            RealArray2D_MatrixMultiply ( False, False, 1.0e+00, w, u, 0.0e+00, v, status ) ;
            if ( View2D_Rows ( m ) == View2D_Columns ( m ) )
            {
                RealArray2D_MatrixMultiply ( False,  mTranspose, 1.0e+00, v, m, 0.0e+00, u, status ) ;
                RealArray2D_MatrixMultiply ( !mTranspose, False, 1.0e+00, m, u, 0.0e+00, w, status ) ;
                AntisymmetricMatrix_CopyFromRealArray2D ( self, w, False, status ) ;
            }
            else
            {
                auto Integer  n ;
                auto RealArray2D uSlice, wSlice ;
                if ( mTranspose ) n = View2D_Rows    ( m ) ;
                else              n = View2D_Columns ( m ) ;
                RealArray2D_View ( u, 0, 0, a->extent, n, 1, 1, False, &uSlice, status ) ;
                RealArray2D_View ( w, 0, 0, n        , n, 1, 1, False, &wSlice, status ) ;
                RealArray2D_MatrixMultiply ( False,  mTranspose, 1.0e+00, v, m      , 0.0e+00, &uSlice, status ) ;
                RealArray2D_MatrixMultiply ( !mTranspose, False, 1.0e+00, m, &uSlice, 0.0e+00, &wSlice, status ) ;
                AntisymmetricMatrix_CopyFromRealArray2D ( self, &wSlice, False, status ) ;
            }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Commutator of two symmetric matrices and two transformations X * A * B * Y - Y^T * B * A * X^T.
! . The transformations can be transposed, i.e X or X^T, Y or Y^T in the first product with the reverse in the second product.
! . Fast version that uses a lot of memory.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_CommutatorXSSY (       AntisymmetricMatrix *self       ,
                                          const SymmetricMatrix     *a          ,
                                          const SymmetricMatrix     *b          ,
                                          const RealArray2D         *x          ,
                                          const RealArray2D         *y          ,
                                          const Boolean              xTranspose ,
                                          const Boolean              yTranspose ,
                                                RealArray2D         *u          ,
                                                RealArray2D         *v          ,
                                                RealArray2D         *w          ,
                                                Status              *status     )
{
    if ( ( self != NULL ) &&
         ( a    != NULL ) &&
         ( b    != NULL ) &&
         ( x    != NULL ) &&
         ( y    != NULL ) &&
         ( u    != NULL ) &&
         ( v    != NULL ) &&
         ( w    != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Boolean isOK ;
        isOK = ( a->extent == b->extent ) ;
        if ( xTranspose ) isOK = isOK && ( ( a->extent == View2D_Rows    ( x ) ) && ( View2D_Columns ( x ) == self->extent ) ) ;
        else              isOK = isOK && ( ( a->extent == View2D_Columns ( x ) ) && ( View2D_Rows    ( x ) == self->extent ) ) ;
        if ( yTranspose ) isOK = isOK && ( ( b->extent == View2D_Columns ( y ) ) && ( View2D_Rows    ( y ) == self->extent ) ) ;
        else              isOK = isOK && ( ( b->extent == View2D_Rows    ( y ) ) && ( View2D_Columns ( y ) == self->extent ) ) ;
        if ( isOK )
        {
            auto RealArray2D uSlice, vSlice ; /* . Use slices for the cases when x or y are not square. */
            /* . A * B in W. */
            SymmetricMatrix_CopyToRealArray2D ( a, u, status ) ;
            SymmetricMatrix_CopyToRealArray2D ( b, v, status ) ;
            RealArray2D_MatrixMultiply ( False, False, 1.0e+00, u, v, 0.0e+00, w, status ) ;
            /* . X(^T) * W in U. */
            RealArray2D_View ( u, 0, 0, self->extent, a->extent, 1, 1, False, &uSlice, status ) ;
            RealArray2D_MatrixMultiply ( xTranspose, False, 1.0e+00, x, w, 0.0e+00, &uSlice, status ) ;
            /* . U * Y(^T) in V. */
            RealArray2D_View ( v, 0, 0, self->extent, self->extent, 1, 1, False, &vSlice, status ) ;
            RealArray2D_MatrixMultiply ( False, yTranspose, 1.0e+00, &uSlice, y, 0.0e+00, &vSlice, status ) ;
            /* . Save the result. */
            AntisymmetricMatrix_CopyFromRealArray2D ( self, &vSlice, False, status ) ;
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copy from a 2D array with antisymmetrization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_CopyFromRealArray2D ( AntisymmetricMatrix *self, const RealArray2D *other, const Boolean scale, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer n = self->extent ;
        if ( ( n != View2D_Rows ( other ) ) || ( n != View2D_Columns ( other ) ) ) Status_Set ( status, Status_NonConformableArrays ) ;
        else
        {
            auto Integer  i, j ;
            for ( i = 0 ; i < n ; i++ )
            {
                for ( j = 0 ; j < i ; j++ ) AntisymmetricMatrix_Item ( self, i, j ) = ( Array2D_Item ( other, i, j ) - Array2D_Item ( other, j, i ) ) ;
            }
            if ( scale ) AntisymmetricMatrix_Scale ( self, 0.5e+00 ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copying.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_CopyTo ( const AntisymmetricMatrix *self, AntisymmetricMatrix *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) && Status_IsOK ( status ) )
    {
        Utilities1D_CopyTo ( self->size, self->data, 1, other->size, other->data, 1, status ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copy to a full 2D array.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_CopyToRealArray2D ( const AntisymmetricMatrix  *self, RealArray2D *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer  n ;
        n = self->extent ;
        if ( ( n != View2D_Rows ( other ) ) || ( n != View2D_Columns ( other ) ) ) Status_Set ( status, Status_NonConformableArrays ) ;
        else
        {
            auto Integer  i, j ;
            auto Real     item ;
            for ( i = 0 ; i < n ; i++ )
            {
                for ( j = 0 ; j < i ; j++ )
                {
                    item = AntisymmetricMatrix_Item ( self, i, j ) ;
                    Array2D_Item ( other, i, j ) =  item ;
                    Array2D_Item ( other, j, i ) = -item ;
                }
                Array2D_Item ( other, i, i ) = 0.0e+00 ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_Deallocate ( AntisymmetricMatrix **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        RealBlock_Deallocate ( &((*self)->block) ) ;
        Memory_Deallocate    (   (*self)         ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Constructor from an extent and, optionally, a block.
!---------------------------------------------------------------------------------------------------------------------------------*/
AntisymmetricMatrix *AntisymmetricMatrix_FromExtentBlock ( const Integer    extent        ,
                                                                 RealBlock *block         ,
                                                           const Boolean    withReference ,
                                                                 Status    *status        )
{
    AntisymmetricMatrix *self = AntisymmetricMatrix_Allocate ( status ) ;
    if ( ( self != NULL ) && ( extent > 0 ) )
    {
        self->extent = extent ;
        self->size   = ( extent * ( extent - 1 ) ) / 2 ;
        if ( block != NULL )
        {
            if ( withReference ) { Array_AssignBlockWithReference    ( self, block, 0 ) ; }
            else                 { Array_AssignBlockWithoutReference ( self, block, 0 ) ; }
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get a column of the matrix.
! . To get a row scale by -1.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_GetColumn ( const AntisymmetricMatrix *self, const Integer  n, RealArray1D *column, Status *status )
{
    if ( ( self != NULL ) && ( column != NULL ) && Status_IsOK ( status ) )
    {
        if ( ( n < 0 ) || ( n >= self->extent ) ) Status_Set ( status, Status_IndexOutOfRange ) ;
        else if ( self->extent != View1D_Extent ( column ) ) Status_Set ( status, Status_NonConformableArrays ) ;
        else
        {
            auto Integer  i ;
            for ( i = 0     ; i < n            ; i++ ) Array1D_Item ( column, i ) = - AntisymmetricMatrix_Item ( self, n, i ) ;
            Array1D_Item ( column, n ) = 0.0e+00 ;
            for ( i = (n+1) ; i < self->extent ; i++ ) Array1D_Item ( column, i ) =   AntisymmetricMatrix_Item ( self, i, n ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real AntisymmetricMatrix_GetItem ( const AntisymmetricMatrix *self, const Integer  i, const Integer  j, Status *status )
{
    Boolean  isOK ;
    Integer  ij   = 0 ;
    Real     sign = 0.0e+00, value = 0.0e+00 ;
    isOK = AntisymmetricMatrix_GetItemIndexAndSign ( self, i, j, &ij, &sign, status ) ;
    if ( isOK ) value = sign * self->data[ij] ;
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the index and sign of an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean AntisymmetricMatrix_GetItemIndexAndSign ( const AntisymmetricMatrix *self   ,
                                                  const Integer              i      ,
                                                  const Integer              j      ,
                                                        Integer             *index  ,
                                                        Real                *sign   ,
                                                        Status              *status )
{
    Boolean isOK = False ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        if ( ( i >= 0 ) && ( i < self->extent ) && ( j >= 0 ) && ( j < self->extent ) )
        {
            if ( i != j )
            {
                if ( i > j )
                {
                    (*index) = ( ( i * ( i - 1 ) ) / 2 + j ) ;
                    (*sign ) =  1.0e+00 ;
                }
                else
                {
                    (*index) = ( ( j * ( j - 1 ) ) / 2 + i ) ;
                    (*sign ) = -1.0e+00 ;
                }
            }
            isOK = True ;
        }
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
    return isOK ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the default iterator.
!---------------------------------------------------------------------------------------------------------------------------------*/
Iterator *AntisymmetricMatrix_MakeIterator ( const AntisymmetricMatrix *self, Status *status )
{
    Iterator *iterator = NULL ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto Iterator1D *iterator1D ; 
        iterator   = Iterator_Allocate   ( status ) ;
        iterator1D = Iterator1D_Allocate ( status ) ;
        Iterator1D_Initialize            ( iterator1D, 0, self->size, 1, status ) ;
        Iterator1D_MakeIterator          ( iterator1D, iterator ) ;
    }
    return iterator ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Printing.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_Print ( const AntisymmetricMatrix *self )
{
    if ( self == NULL ) printf ( "Null antisymmetric matrix.\n" ) ;
    else
    {
        auto Integer  i, j ;
        for ( i = 0 ; i < self->extent ; i++ )
        {
            for ( j = 0 ; j < i ; j++ ) printf ( "%15.10f", AntisymmetricMatrix_Item ( self, i, j ) ) ;
            printf ( "   0.0000000000\n" ) ;
        }
        printf ( "\n" ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Scale the items of a matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_Scale ( AntisymmetricMatrix *self, const Real value )
{
    if ( self != NULL ) Utilities1D_Scale ( self->size, self->data, 1, value ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set the items of a matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_Set ( AntisymmetricMatrix *self, const Real value )
{
    if ( self != NULL ) Utilities1D_Set ( self->size, self->data, 1, value ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_SetItem ( AntisymmetricMatrix *self, const Integer  i, const Integer  j, const Real value, Status *status )
{
    Boolean  isOK ;
    Integer  ij = 0 ;
    Real     sign = 0.0e+00 ;
    isOK = AntisymmetricMatrix_GetItemIndexAndSign ( self, i, j, &ij, &sign, status ) ;
    if ( isOK ) self->data[ij] = sign * value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Transform by a symmetric matrix (S * A * S).
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_SymmetricTransform  ( const AntisymmetricMatrix  *self, const SymmetricMatrix *matrix, AntisymmetricMatrix *result, Status *status )
{
    if ( ( self != NULL ) && ( matrix != NULL ) && ( result != NULL ) && Status_IsOK ( status ) )
    {
        if ( ( self->extent == matrix->extent ) && ( matrix->extent == result->extent ) )
        {
            auto Integer      n ;
            auto RealArray1D *icolumn = NULL, *SAi = NULL, *tcolumn = NULL ;
            n = self->extent ;
            icolumn = RealArray1D_AllocateWithExtent ( n, status ) ;
            SAi     = RealArray1D_AllocateWithExtent ( n, status ) ;
            tcolumn = RealArray1D_AllocateWithExtent ( n, status ) ;
            if ( ( icolumn != NULL ) && ( SAi != NULL ) && ( tcolumn != NULL ) )
            {
                auto Integer  i, k, l ;
                for ( i = 1 ; i < n ; i++ )
                {
                    SymmetricMatrix_GetColumn ( matrix, i, icolumn, NULL ) ;
                    for ( k = 0 ; k < n ; k++ )
                    {
                        AntisymmetricMatrix_GetColumn ( self, k, tcolumn, NULL ) ;
                        Array1D_Item ( SAi, k ) = RealArray1D_Dot ( icolumn, tcolumn, NULL ) ;
                    }
                    for ( l = 0 ; l < i ; l++ )
                    {
                        SymmetricMatrix_GetColumn ( matrix, l, tcolumn, NULL ) ;
                        AntisymmetricMatrix_Item ( result, i, l ) = RealArray1D_Dot ( SAi, tcolumn, NULL ) ;
                    }
                }
            }
            else Status_Set ( status, Status_OutOfMemory ) ;
            RealArray1D_Deallocate ( &icolumn ) ;
            RealArray1D_Deallocate ( &SAi     ) ;
            RealArray1D_Deallocate ( &tcolumn ) ;
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the trace of the product of two matrices.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real AntisymmetricMatrix_TraceOfProduct ( const AntisymmetricMatrix *self, const AntisymmetricMatrix *other, Status *status )
{
    Real value = 0.0e+00 ;
    if ( ( self != NULL ) && ( other != NULL ) && Status_IsOK ( status ) )
    {
        Utilities1D_Dot ( self->size, self->data, 1, other->size, other->data, 1, value, status ) ;
        value *= - 2.0e+00 ;
    }
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Transform by a matrix (B^T * A * B) or its transpose (B * A * B^T).
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_Transform ( const AntisymmetricMatrix *self, const RealArray2D *matrix, const Boolean useTranspose, AntisymmetricMatrix *result, Status *status )
{
    if ( ( self != NULL ) && ( matrix != NULL ) && ( result != NULL ) && Status_IsOK ( status ) )
    {
        auto RealArray1D *BAik = NULL, *column = NULL ;
        BAik   = RealArray1D_AllocateWithExtent ( self->extent, status ) ;
        column = RealArray1D_AllocateWithExtent ( self->extent, status ) ;
        if ( ( BAik != NULL ) && ( column != NULL ) )
        {
            auto Integer     i, k, l ;
            auto RealArray1D iview, lview ;
            /* . (Bij * Ajk * (B^T)kl = Bij * Ajk * Blk). */
            if ( useTranspose )
            {
                if ( ( self->extent == View2D_Columns ( matrix ) ) && ( result->extent == View2D_Rows ( matrix ) ) )
                {
                    for ( i = 0 ; i < View2D_Rows ( matrix ) ; i++ )
                    {
                        RealArray2D_RowView ( matrix, i, False, &iview, NULL ) ;
                        for ( k = 0 ; k < self->extent ; k++ )
                        {
                            AntisymmetricMatrix_GetColumn ( self, k, column, status ) ;
                            Array1D_Item ( BAik, k ) = RealArray1D_Dot ( column, &iview, NULL ) ;
                        }
                        for ( l = 0 ; l < i ; l++ )
                        {
                            RealArray2D_RowView ( matrix, l, False, &lview, NULL ) ;
                            AntisymmetricMatrix_Item ( result, i, l ) = RealArray1D_Dot ( BAik, &lview, NULL ) ;
                        }
                    }
                }
                else Status_Set ( status, Status_NonConformableArrays ) ;
            }
            /* . ((B^T)ij * Ajk * Bkl = Bji * Ajk * Bkl). */
            else
            {
                if ( ( self->extent == View2D_Rows ( matrix ) ) && ( result->extent == View2D_Columns ( matrix ) ) )
                {
                    for ( i = 0 ; i < View2D_Columns ( matrix ) ; i++ )
                    {
                        RealArray2D_ColumnView ( matrix, i, False, &iview, NULL ) ;
                        for ( k = 0 ; k < self->extent ; k++ )
                        {
                            AntisymmetricMatrix_GetColumn ( self, k, column, status ) ;
                            Array1D_Item ( BAik, k ) = RealArray1D_Dot ( column, &iview, NULL ) ;
                        }
                        for ( l = 0 ; l < i ; l++ )
                        {
                            RealArray2D_ColumnView ( matrix, l, False, &lview, NULL ) ;
                            AntisymmetricMatrix_Item ( result, i, l ) = RealArray1D_Dot ( BAik, &lview, NULL ) ;
                        }
                    }
                }
                else Status_Set ( status, Status_NonConformableArrays ) ;
            }
        }
        else Status_Set ( status, Status_OutOfMemory ) ;
        RealArray1D_Deallocate ( &BAik   ) ;
        RealArray1D_Deallocate ( &column ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Transpose the matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void AntisymmetricMatrix_Transpose ( AntisymmetricMatrix *self ) { AntisymmetricMatrix_Scale ( self, -1.0e+00 ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . The view size.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer AntisymmetricMatrix_ViewSize ( const Integer extent )
{
    Integer size = 0 ;
    if ( extent > 0 ) size = ( extent * ( extent - 1 ) ) / 2 ;
    return size ;
}

# undef _UseCBLAS
# undef _UseReal

