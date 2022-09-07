/*==================================================================================================================================
! . Real symmetric matrices.
!=================================================================================================================================*/

# include <stdio.h>
# include <stdlib.h>

# include "Array_Macros.h"
# include "Iterator1D.h"
# include "Memory.h"
# include "NumericalMacros.h"
# include "SymmetricMatrix.h"

# define _UseCBLAS
# define _UseReal
# include "Utilities1D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . The maximum absolute value.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real SymmetricMatrix_AbsoluteMaximum ( const SymmetricMatrix *self )
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
void SymmetricMatrix_Add ( SymmetricMatrix *self, const Real alpha, const SymmetricMatrix *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL )  && Status_IsOK ( status ) )
    {
        Utilities1D_Add ( self->size, self->data, 1, other->size, other->data, 1, alpha, status ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Raw allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
SymmetricMatrix *SymmetricMatrix_Allocate ( Status *status )
{
    SymmetricMatrix *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( SymmetricMatrix ) ;
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
SymmetricMatrix *SymmetricMatrix_AllocateWithExtent ( const Integer extent, Status *status )
{
    SymmetricMatrix *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        auto RealBlock *block = NULL ;
        block = RealBlock_Allocate ( SymmetricMatrix_ViewSize ( extent ), status ) ;
        self  = SymmetricMatrix_FromExtentBlock ( extent, block, True   , status ) ;
        if ( self == NULL ) RealBlock_Deallocate ( &block ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Anticommutator of two symmetric matrices A * B + B * A.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_AnticommutatorSS ( SymmetricMatrix *self, const SymmetricMatrix *a, const SymmetricMatrix *b, Status *status )
{
    if ( ( self != NULL ) && ( a != NULL ) && ( b != NULL ) && Status_IsOK ( status ) )
    {
        if ( ( self->extent == a->extent ) && ( a->extent == b->extent ) )
        {
            auto Integer  i, j, k ;
            auto Real    sum ;
            for ( i = 0 ; i < self->extent ; i++ )
            {
                for ( k = 0 ; k <= i ; k++ )
                {
                    sum = 0.0e+00 ;
                    for ( j = 0         ; j <= k           ; j++ ) sum += ( SymmetricMatrix_Item ( a, i, j ) * SymmetricMatrix_Item ( b, k, j ) + SymmetricMatrix_Item ( b, i, j ) * SymmetricMatrix_Item ( a, k, j ) ) ;
                    for ( j = ( k + 1 ) ; j <= i           ; j++ ) sum += ( SymmetricMatrix_Item ( a, i, j ) * SymmetricMatrix_Item ( b, j, k ) + SymmetricMatrix_Item ( b, i, j ) * SymmetricMatrix_Item ( a, j, k ) ) ;
                    for ( j = ( i + 1 ) ; j < self->extent ; j++ ) sum += ( SymmetricMatrix_Item ( a, j, i ) * SymmetricMatrix_Item ( b, j, k ) + SymmetricMatrix_Item ( b, j, i ) * SymmetricMatrix_Item ( a, j, k ) ) ;
                    SymmetricMatrix_Item ( self, i, k ) = sum ;
                }
            }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deep cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
SymmetricMatrix *SymmetricMatrix_CloneDeep ( const SymmetricMatrix *self, Status *status )
{
    SymmetricMatrix *clone = NULL ;
    if ( self != NULL )
    {
        clone = SymmetricMatrix_AllocateWithExtent ( self->extent, status ) ;
        SymmetricMatrix_CopyTo ( self, clone, status ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Shallow cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
SymmetricMatrix *SymmetricMatrix_CloneShallow ( const SymmetricMatrix *self, Status *status )
{
    SymmetricMatrix *clone = NULL ;
    if ( self != NULL ) clone = SymmetricMatrix_FromExtentBlock ( self->extent, self->block, True, status ) ;
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copy from a full 2D array ensuring symmetry.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_CopyFromRealArray2D ( SymmetricMatrix *self, const RealArray2D *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer  n = self->extent ;
        if ( ( n == View2D_Columns ( other ) ) && ( n == View2D_Rows ( other ) ) )
        {
            auto Integer  i, j ;
            for ( i = 0 ; i < n ; i++ )
            {
                for ( j = 0 ; j <= i ; j++ ) SymmetricMatrix_Item ( self, i, j ) = 0.5e+00 * ( Array2D_Item ( other, i, j ) + Array2D_Item ( other, j, i ) ) ;
            }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copying.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_CopyTo ( const SymmetricMatrix *self, SymmetricMatrix *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) && Status_IsOK ( status ) )
    {
        Utilities1D_CopyTo ( self->size, self->data, 1, other->size, other->data, 1, status ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copy to a full 2D array.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_CopyToRealArray2D ( const SymmetricMatrix *self, RealArray2D *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer n = self->extent ;
        if ( ( n == View2D_Columns ( other ) ) && ( n == View2D_Rows ( other ) ) )
        {
            auto Integer  i, j ;
            auto Real     item ;
            for ( i = 0 ; i < n ; i++ )
            {
                for ( j = 0 ; j < i ; j++ )
                {
                    item = SymmetricMatrix_Item ( self, i, j ) ;
                    Array2D_Item ( other, i, j ) = item ;
                    Array2D_Item ( other, j, i ) = item ;
                }
                Array2D_Item ( other, i, i ) = SymmetricMatrix_Item ( self, i, i ) ;
            }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Deallocate ( SymmetricMatrix **self )
{
    if ( (*self) != NULL )
    {
        RealBlock_Deallocate ( &((*self)->block) ) ;
        Memory_Deallocate    (   (*self)         ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Diagonal of product.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_DiagonalOfProduct ( const SymmetricMatrix *self   ,
                                         const SymmetricMatrix *other  ,
                                               RealArray1D     *result ,
                                               Status          *status )
{
    if ( ( self != NULL ) && ( other != NULL ) && ( result != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer n = self->extent ;
        if ( ( n == other->extent ) && ( n == View1D_Extent ( result ) ) )
        {
            auto Integer i, j ;
            auto Real    d ;
            for ( i = 0 ; i < n ; i++ )
            {
                d = 0.0e+00 ;
                for ( j = 0   ; j <= i ; j++ ) d += ( SymmetricMatrix_Item ( self, i, j ) * SymmetricMatrix_Item ( other, i, j ) ) ;
                for ( j = i+1 ; j <  n ; j++ ) d += ( SymmetricMatrix_Item ( self, j, i ) * SymmetricMatrix_Item ( other, j, i ) ) ;
                Array1D_Item ( result, i ) = d ;
            }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Diagonal of transform.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_DiagonalOfTransform ( const SymmetricMatrix *self         ,
                                           const RealArray2D     *matrix       ,
                                           const Boolean          useTranspose ,
                                                 RealArray1D     *result       ,
                                                 Status          *status       )
{
    if ( ( self != NULL ) && ( matrix != NULL ) && ( result != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer m = View1D_Extent ( result ), n = self->extent ;
        if ( (     useTranspose   && ( m == View2D_Rows ( matrix ) ) && ( n == View2D_Columns ( matrix ) ) ) ||
             ( ( ! useTranspose ) && ( n == View2D_Rows ( matrix ) ) && ( m == View2D_Columns ( matrix ) ) ) )
        {
            auto Integer i ;
            auto RealArray1D mView, *vector ;
            vector = RealArray1D_AllocateWithExtent ( n, status ) ;
            for ( i = 0 ; i < m ; i++ )
            {
                if ( useTranspose ) RealArray2D_RowView    ( matrix, i, False, &mView, status ) ;
                else                RealArray2D_ColumnView ( matrix, i, False, &mView, status ) ;
                SymmetricMatrix_VectorMultiply ( self, &mView, vector, status ) ;
                Array1D_Item ( result, i ) = RealArray1D_Dot ( &mView, vector, status ) ;
            }
            RealArray1D_Deallocate ( &vector ) ;
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Constructor from an extent and, optionally, a block.
!---------------------------------------------------------------------------------------------------------------------------------*/
SymmetricMatrix *SymmetricMatrix_FromExtentBlock ( const Integer    extent        ,
                                                         RealBlock *block         ,
                                                   const Boolean    withReference ,
                                                         Status    *status        )
{
    SymmetricMatrix *self = SymmetricMatrix_Allocate ( status ) ;
    if ( ( self != NULL ) && ( extent > 0 ) )
    {
        self->extent = extent ;
        self->size   = ( extent * ( extent + 1 ) ) / 2 ;
        if ( block != NULL )
        {
            if ( withReference ) { Array_AssignBlockWithReference    ( self, block, 0 ) ; }
            else                 { Array_AssignBlockWithoutReference ( self, block, 0 ) ; }
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get a column (or a row) of the matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_GetColumn ( const SymmetricMatrix *self, const Integer  n, RealArray1D *column, Status *status )
{
    if ( ( self != NULL ) && ( column != NULL ) && Status_IsOK ( status ) )
    {
             if ( self->extent <= n                         ) Status_Set ( status, Status_IndexOutOfRange      ) ;
        else if ( self->extent != View1D_Extent ( column ) ) Status_Set ( status, Status_NonConformableArrays ) ;
        else
        {
            auto Integer  i ;
            for ( i = 0     ; i <= n           ; i++ ) Array1D_Item ( column, i ) = SymmetricMatrix_Item ( self, n, i ) ;
            for ( i = (n+1) ; i < self->extent ; i++ ) Array1D_Item ( column, i ) = SymmetricMatrix_Item ( self, i, n ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real SymmetricMatrix_GetItem ( const SymmetricMatrix *self, const Integer  i, const Integer  j, Status *status )
{
    Real value = 0.0e+00 ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        if ( ( i >= 0 ) && ( i < self->extent ) && ( j >= 0 ) && ( j < self->extent ) )
        {
            if ( i >= j ) value = SymmetricMatrix_Item ( self, i, j ) ;
            else          value = SymmetricMatrix_Item ( self, j, i ) ;
        }
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Add a scalar.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Increment ( SymmetricMatrix *self, const Real value )
{
    if ( self != NULL ) Utilities1D_Increment ( self->size, self->data, 1, value ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copy selected elements of the matrix to a square form.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_IndexedCopyToRealArray2D ( const SymmetricMatrix *self, const IntegerBlock *indices, RealArray2D *target, Status *status )
{
    if ( ( self != NULL ) && ( indices != NULL ) && ( target != NULL ) )
    {
        auto Integer i, j, m, n ;
        auto Real    v ;
        for ( i = 0 ; i < indices->capacity ; i++ )
        {
            m = Block_Item ( indices, i ) ;
            for ( j = 0 ; j <= i ; j++ )
            {
                n = Block_Item  ( indices, j ) ;
                v = SymmetricMatrix_Item ( self, m, n ) ;
                Array2D_Item ( target, i, j ) = v ;
                Array2D_Item ( target, j, i ) = v ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for a diagonal matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean SymmetricMatrix_IsDiagonal ( const SymmetricMatrix *self, const Real tolerance )
{
    Boolean isDiagonal = False ;
    if ( self != NULL )
    {
        auto Integer  i, ij, j ;
        isDiagonal = True ;
        for ( i = 0, ij = 0 ; i < self->extent ; i++, ij++ )
        {
            for ( j = 0 ; j < i ; ij++, j++ )
            {
                if ( fabs ( self->data[ij] ) > tolerance ) { isDiagonal = False ; break ; }
            }
            if ( ! isDiagonal ) break ;
        }
    }
    return isDiagonal ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Constructor from a set of eigenvalues and eigenvectors.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_MakeFromEigensystem (       SymmetricMatrix *self            ,
                                           const Boolean          zeroMatrix      ,
                                           const Integer          numberOfVectors ,
                                           const RealArray1D     *eigenValues     ,
                                           const RealArray2D     *eigenVectors    ,
                                                 Status          *status          )
{
    if ( zeroMatrix ) SymmetricMatrix_Set ( self, 0.0e+00 ) ;
    if ( ( self != NULL ) && ( eigenValues != NULL ) && ( eigenVectors != NULL ) && ( numberOfVectors != 0 ) && Status_IsOK ( status ) )
    {
        auto Integer  extent ;
        extent = View2D_Rows ( eigenVectors ) ;
        if ( ( self->extent    == extent                               ) &&
             ( numberOfVectors <= View1D_Extent  ( eigenValues  ) ) &&
             ( numberOfVectors <= View2D_Columns ( eigenVectors ) ) )
        {
# ifdef USEOPENMP
            #pragma omp parallel shared ( status )
# endif
            {
                auto Integer  i, j, o ;
                auto Real     sum, *work ;
                work = Memory_AllocateArrayOfTypes ( numberOfVectors, Real ) ;
                if ( work != NULL )
                {
# ifdef USEOPENMP
		    #pragma omp for schedule ( dynamic )
# endif
                    for ( i = 0 ; i < self->extent ; i++ )
                    {
                        for ( o = 0 ; o <  numberOfVectors ; o++ ) work[o] = Array1D_Item ( eigenValues, o ) * Array2D_Item ( eigenVectors, i, o ) ;
                        for ( j = 0 ; j <= i               ; j++ )
                        {
                            for ( o = 0, sum = 0.0e+00 ; o < numberOfVectors ; o++ ) sum += work[o] * Array2D_Item ( eigenVectors, j, o ) ;
                            SymmetricMatrix_Item ( self, i, j ) += sum ;
                        }
                    }
                    Memory_Deallocate ( work ) ;
                }
                else Status_Set ( status, Status_OutOfMemory ) ;
            }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the default iterator.
!---------------------------------------------------------------------------------------------------------------------------------*/
Iterator *SymmetricMatrix_MakeIterator ( const SymmetricMatrix *self, Status *status )
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
! . Post-matrix multiply.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_PostMatrixMultiply ( const SymmetricMatrix *self, const RealArray2D *matrix, const Boolean useTranspose, RealArray2D *result, Status *status )
{
    if ( ( self != NULL ) && ( matrix != NULL ) && ( result != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer  i ;
        if ( useTranspose )
        {
            if ( ( self->extent                == View2D_Columns ( matrix ) ) &&
                 ( self->extent                == View2D_Rows    ( result ) ) &&
                 ( View2D_Rows ( matrix ) == View2D_Columns ( result ) ) )
            {
                for ( i = 0 ; i < View2D_Rows ( matrix ) ; i++ )
                {
                    cblas_dspmv ( CblasColMajor, CblasUpper, self->extent, 1.0e+00, self->data, Array2D_ItemPointer ( matrix, i, 0 ), matrix->stride1 ,
                                                                           0.0e+00,             Array2D_ItemPointer ( result, 0, i ), result->stride0 ) ;
                }
            }
            else Status_Set ( status, Status_NonConformableArrays ) ;
        }
        else
        {
            if ( ( self->extent                   == View2D_Rows    ( matrix ) ) &&
                 ( self->extent                   == View2D_Rows    ( result ) ) &&
                 ( View2D_Columns ( matrix ) == View2D_Columns ( result ) ) )
            {
                for ( i = 0 ; i < View2D_Columns ( matrix ) ; i++ )
                {
                    cblas_dspmv ( CblasColMajor, CblasUpper, self->extent, 1.0e+00, self->data, Array2D_ItemPointer ( matrix, 0, i ), matrix->stride0 , 
                                                                           0.0e+00,             Array2D_ItemPointer ( result, 0, i ), result->stride0 ) ;
                }
            }
            else Status_Set ( status, Status_NonConformableArrays ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Pre-matrix multiply.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_PreMatrixMultiply ( const SymmetricMatrix *self, const RealArray2D *matrix, const Boolean useTranspose, RealArray2D *result, Status *status )
{
    if ( ( self != NULL ) && ( matrix != NULL ) && ( result != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer  i ;
        if ( useTranspose )
        {
            if ( ( self->extent              == View2D_Rows    ( matrix ) ) &&
                 ( self->extent              == View2D_Columns ( result ) ) &&
                 ( View2D_Columns ( matrix ) == View2D_Rows    ( result ) ) )
            {
                for ( i = 0 ; i < View2D_Columns ( matrix ) ; i++ )
                {
                    cblas_dspmv ( CblasColMajor, CblasUpper, self->extent, 1.0e+00, self->data, Array2D_ItemPointer ( matrix, 0, i ), matrix->stride0 ,
                                                                           0.0e+00,             Array2D_ItemPointer ( result, i, 0 ), result->stride1 ) ;
                }
            }
            else Status_Set ( status, Status_NonConformableArrays ) ;
        }
        else
        {
            if ( ( self->extent                == View2D_Columns ( matrix ) ) &&
                 ( self->extent                == View2D_Columns ( result ) ) &&
                 ( View2D_Rows ( matrix ) == View2D_Rows    ( result ) ) )
            {
                for ( i = 0 ; i < View2D_Rows ( matrix ) ; i++ )
                {
                    cblas_dspmv ( CblasColMajor, CblasUpper, self->extent, 1.0e+00, self->data, Array2D_ItemPointer ( matrix, i, 0 ), matrix->stride1 ,
                                                                           0.0e+00,             Array2D_ItemPointer ( result, i, 0 ), result->stride1 ) ;
                }
            }
            else Status_Set ( status, Status_NonConformableArrays ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Printing - debugging only.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Print ( const SymmetricMatrix *self )
{
    if ( self == NULL ) printf ( "Null symmetric matrix.\n" ) ;
    else
    {
        auto Integer  i, j ;
        for ( i = 0 ; i < self->extent ; i++ )
        {
            for ( j = 0 ; j <= i ; j++ ) printf ( "%15.10f", SymmetricMatrix_Item ( self, i, j ) ) ;
            printf ( "\n" ) ;
        }
        printf ( "\n" ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Create a projection matrix of the form ( 1 - V V^T ).
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_ProjectionMatrix ( SymmetricMatrix *self, const RealArray2D *vectors, Status *status )
{
    if ( ( self != NULL ) && ( vectors != NULL ) && Status_IsOK ( status ) )
    {

        auto Integer r = View2D_Rows ( vectors ) ;
        if ( SymmetricMatrix_Extent ( self ) == r )
        {
             auto Integer  c, i, ij, j, v ;
             auto Real     sum ;
             c = View2D_Columns ( vectors ) ;
             for ( i = 0, ij = 0 ; i < r ; i++ )
             {
                 for ( j = 0 ; j <= i ; ij++, j++ )
                 {
                     for ( v = 0, sum = 0.0e+00 ; v < c ; v++ ) sum += ( Array2D_Item ( vectors, i, v ) * Array2D_Item ( vectors, j, v ) ) ;
                     self->data[ij] = - sum ;
                     if ( i == j ) self->data[ij] += 1.0e+00 ;
                 }
             }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Project out a set of vectors from a symmetric matrix using ( 1 - V V^T ) S ( 1 - V V^T ).
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_ProjectOut ( SymmetricMatrix *self, const RealArray2D *vectors, SymmetricMatrix *result, Status *status )
{
    if ( ( self != NULL ) && ( vectors != NULL ) && ( result != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer e = SymmetricMatrix_Extent ( self ) ;
        if ( ( e == View2D_Rows           ( vectors ) ) &&
             ( e == SymmetricMatrix_Extent ( result  ) ) )
        {
            auto SymmetricMatrix *projection ;
            projection = SymmetricMatrix_AllocateWithExtent ( e, status ) ;
            SymmetricMatrix_ProjectionMatrix   ( projection, vectors, status ) ;
            SymmetricMatrix_SymmetricTransform ( self, projection, result, status ) ;
            SymmetricMatrix_Deallocate ( &projection ) ;
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Add alpha V V^T to a symmetric matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Raise ( SymmetricMatrix *self, const RealArray2D *vectors, const Real value, Status *status )
{
    if ( ( self != NULL ) && ( vectors != NULL ) && ( value != 0.0e+00 ) && Status_IsOK ( status ) )
    {
        auto Integer  c, r ;
        c = View2D_Columns ( vectors ) ;
        r = View2D_Rows    ( vectors ) ;
        if ( ( self->extent == r ) && ( c > 0 ) )
        {
            auto Integer  i, ij, j, v ;
            auto Real     sum         ;
            for ( i = 0, ij = 0 ; i < r ; i++ )
            {
                for ( j = 0 ; j <= i ; ij++, j++ )
                {
                   for ( v = 0, sum = 0.0e+00 ; v < c ; v++ ) sum += Array2D_Item ( vectors, i, v ) * Array2D_Item ( vectors, j, v ) ;
                   self->data[ij] += value * sum ;
                }
            }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Rank-1 update.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . self needs to be compact. */
void SymmetricMatrix_Rank1Update ( SymmetricMatrix *self, const Real alpha, const RealArray1D *vector, Status *status )
{
    if ( ( self != NULL ) && ( alpha != 0.0e+00 ) && ( vector != NULL ) && Status_IsOK ( status ) )
    {
        if ( self->extent == View1D_Extent ( vector ) )
        {
            cblas_dspr ( CblasColMajor, CblasUpper, self->extent, alpha, Array1D_Data ( vector ), vector->stride, SymmetricMatrix_Data ( self ) ) ;
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Root mean square.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real SymmetricMatrix_RootMeanSquare ( const SymmetricMatrix *self )
{
    Real rms = 0.0e+00 ;
    if ( self != NULL ) { rms = sqrt ( SymmetricMatrix_TraceOfProduct ( self, self, NULL ) ) / ( Real ) ( self->extent * self->extent ) ; }
    return rms ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Scaling.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Scale ( SymmetricMatrix *self, const Real value )
{
    if ( self != NULL ) Utilities1D_Scale ( self->size, self->data, 1, value ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Scale the diagonal elements.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_ScaleDiagonal ( SymmetricMatrix *self, const Real value )
{
    if ( self != NULL )
    {
        auto Integer i ;
        for ( i = 0 ; i < self->extent ; i++ ) { SymmetricMatrix_Item ( self, i, i ) *= value ; }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Scale the off-diagonal elements.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_ScaleOffDiagonal ( SymmetricMatrix *self, const Real value )
{
    if ( self != NULL )
    {
        auto Integer  i, j ;
        for ( i = 0 ; i < self->extent ; i++ )
        {
            for ( j = 0 ; j < i ; j++ ) SymmetricMatrix_Item ( self, i, j ) *= value ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Setting.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Set ( SymmetricMatrix *self, const Real value )
{
    if ( self != NULL ) Utilities1D_Set ( self->size, self->data, 1, value ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set a column (or a row) of the matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_SetColumn ( SymmetricMatrix *self, const Integer  n, const RealArray1D *column, Status *status )
{
    if ( ( self != NULL ) && ( column != NULL ) && Status_IsOK ( status ) )
    {
             if ( self->extent <= n                             ) Status_Set ( status, Status_IndexOutOfRange      ) ;
        else if ( self->extent != View1D_Extent ( column ) ) Status_Set ( status, Status_NonConformableArrays ) ;
        else
        {
            auto Integer  i ;
            for ( i = 0     ; i <= n           ; i++ ) SymmetricMatrix_Item ( self, n, i ) = Array1D_Item ( column, i ) ;
            for ( i = (n+1) ; i < self->extent ; i++ ) SymmetricMatrix_Item ( self, i, n ) = Array1D_Item ( column, i ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_SetItem ( SymmetricMatrix *self, const Integer  i, const Integer  j, const Real value, Status *status )
{
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        if ( ( i >= 0 ) && ( i < self->extent ) && ( j >= 0 ) && ( j < self->extent ) )
        {
            if ( i >= j ) SymmetricMatrix_Item ( self, i, j ) = value ;
            else          SymmetricMatrix_Item ( self, j, i ) = value ;
        }
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Sparsity (as a percentage).
! . Sparsity is defined as the number of items with a magnitude less than a given tolerance.
! . The full N^2 dimension of the matrix is taken into account.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real SymmetricMatrix_Sparsity ( const SymmetricMatrix *self, const Real tolerance )
{
    Real sparsity = 1.0e+00 ;
    if ( ( self != NULL ) && ( self->extent > 0 ) )
    {
        /* . Off-diagonal elements count twice. */
        auto Integer  i, ij, j, n = 0 ;
        for ( i = 0, ij = 0 ; i < self->extent ; i++, ij++ )
        {
            for ( j = 0 ; j < i ; ij++, j++ )
            {
                if ( fabs ( self->data[ij] ) <= tolerance ) n += 2 ;
            }
            if ( fabs ( self->data[ij] ) <= tolerance ) n += 1 ;
        }
        sparsity = ( ( Real ) n ) / ( ( Real ) ( self->extent * self->extent ) ) ;
    }
    return ( 100.0e+00 * sparsity ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Convert self and other to ( self+other ) and ( self-other ).
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_SumDifference ( SymmetricMatrix *self, SymmetricMatrix *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) && Status_IsOK ( status ) )
    {
        if ( self->size == other->size )
        {
            auto Integer  i ;
            auto Real     s, t ;
            for ( i = 0 ; i < self->size ; i++ )
            {
                t = self->data[i] + other->data[i] ;
                s = self->data[i] - other->data[i] ;
                self ->data[i] = t ;
                other->data[i] = s ;
            }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the product of two symmetric matrices (self * other).
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_SymmetricMatrixMultiply ( const SymmetricMatrix *self, const SymmetricMatrix *other, RealArray2D *result, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) && ( result != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer  n = self->extent ;
        if ( ( other->extent == n ) && ( View2D_Rows ( result ) == n ) && ( View2D_Columns ( result ) == n ) )
        {
            auto Integer  i, ij, j, jk, k ;
            auto Real     sum ;
            /* . Aij * Bjk = Cik. */
            for ( i = 0 ; i < n ; i++ )
            {
                for ( k = 0 ; k < n ; k++ )
                {
                    for ( j = 0, sum = 0.0e+00 ; j < n ; j++ )
                    {
                        if ( i < j ) ij = ( j * ( j + 1 ) ) / 2 + i ;
                        else         ij = ( i * ( i + 1 ) ) / 2 + j ;
                        if ( k < j ) jk = ( j * ( j + 1 ) ) / 2 + k ;
                        else         jk = ( k * ( k + 1 ) ) / 2 + j ;
                        sum += self->data[ij] * other->data[jk] ;
                    }
                    Array2D_Item ( result, i, k ) = sum ;
                }
            }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Transform a matrix by another symmetric matrix (i.e. calculate B * A * B).
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_SymmetricTransform ( const SymmetricMatrix *A, const SymmetricMatrix *B, SymmetricMatrix *BAB, Status *status )
{
    if ( ( A != NULL ) && ( B != NULL ) && ( BAB != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer n = A->extent ;
        if ( ( B->extent == n ) && ( BAB->extent == n ) )
        {
            auto RealArray1D *BA ;
            BA = RealArray1D_AllocateWithExtent ( n, status ) ;
            if ( BA != NULL )
            {
                auto Integer  i, j, k, l ;
                auto Real     sum ;
                for ( i = 0 ; i < n ; i++ )
                {
                    for ( k = 0 ; k < n ; k++ )
                    {
                        for ( j = 0, sum = 0.0e+00 ; j < n ; j++ )
                        {
                            sum += SymmetricMatrix_GetItem ( B, i, j, NULL ) * SymmetricMatrix_GetItem ( A, j, k, NULL ) ;
                        }
                        Array1D_Item ( BA, k ) = sum ;
                    }
                    for ( l = 0 ; l <= i ; l++ )
                    {
                        for ( k = 0, sum = 0.0e+00 ; k < n ; k++ )
                        {
                            sum += Array1D_Item ( BA, k ) * SymmetricMatrix_GetItem ( B, k, l, NULL ) ;
                        }
                        SymmetricMatrix_Item ( BAB, i, l ) = sum ;
                    }
                }
                RealArray1D_Deallocate ( &BA ) ;
            }
            else Status_Set ( status, Status_OutOfMemory ) ;
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . The trace.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real SymmetricMatrix_Trace ( const SymmetricMatrix *self )
{
    Real trace = 0.0e+00 ;
    if ( self != NULL )
    {
        auto Integer  i ;
        for ( i = 0 ; i < self->extent ; i++ ) { trace += SymmetricMatrix_Item ( self, i, i ) ; }
    }
    return trace ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . The trace of the product of two symmetric matrices.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real SymmetricMatrix_TraceOfProduct ( const SymmetricMatrix *self, const SymmetricMatrix *other, Status *status )
{
    Real value = 0.0e+00 ;
    if ( ( self != NULL ) && ( other != NULL ) && Status_IsOK ( status ) )
    {
        if ( self->extent == other->extent )
        {
            auto Integer  i, ij, j ;
            for ( i = 0, ij = 0 ; i < self->extent ; i++, ij++ )
            {
                for ( j = 0 ; j < i ; ij++, j++ ) value += self->data[ij] * other->data[ij] ;
                value += 0.5e+00 * self->data[ij] * other->data[ij] ;
            }
            value *= 2.0e+00 ;
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Transform by a matrix (B^T * A * B) or its transpose (B * A * B^T).
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_Transform ( const SymmetricMatrix *self, const RealArray2D *matrix, const Boolean useTranspose, SymmetricMatrix *result, Status *status )
{
    if ( ( self != NULL ) && ( matrix != NULL ) && ( result != NULL ) && Status_IsOK ( status ) )
    {
        auto Real *tmpvec = Memory_AllocateArrayOfTypes ( self->extent, Real ) ;
        if ( tmpvec == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
        else
        {
            auto Integer  i, il, l ;
            Memory_Set ( tmpvec, self->extent, 0.0e+00 ) ;
            if ( useTranspose )
            {
                if ( ( self->extent == View2D_Columns ( matrix ) ) && ( result->extent == View2D_Rows ( matrix ) ) )
                {
                    for ( i = 0, il = 0 ; i < View2D_Rows ( matrix ) ; i++ )
                    {
                        cblas_dspmv ( CblasColMajor, CblasUpper, self->extent, 1.0e+00, self->data, Array2D_ItemPointer ( matrix, i, 0 ), matrix->stride1, 0.0e+00, tmpvec, 1 ) ;
                        for ( l = 0 ; l <= i ; il++, l++ ) result->data[il] = cblas_ddot ( self->extent, tmpvec, 1, Array2D_ItemPointer ( matrix, l, 0 ), matrix->stride1 ) ;
                    }
                }
                else Status_Set ( status, Status_NonConformableArrays ) ;
            }
            else
            {
                if ( ( self->extent == View2D_Rows ( matrix ) ) && ( result->extent == View2D_Columns ( matrix ) ) )
                {
                    for ( i = 0, il = 0 ; i < View2D_Columns ( matrix ) ; i++ )
                    {
                        cblas_dspmv ( CblasColMajor, CblasUpper, self->extent, 1.0e+00, self->data, Array2D_ItemPointer ( matrix, 0, i ), matrix->stride0, 0.0e+00, tmpvec, 1 ) ;
                        for ( l = 0 ; l <= i ; il++, l++ ) result->data[il] = cblas_ddot ( self->extent, tmpvec, 1, Array2D_ItemPointer ( matrix, 0, l ), matrix->stride0 ) ;
                    }
                }
                else Status_Set ( status, Status_NonConformableArrays ) ;
            }
            Memory_Deallocate ( tmpvec ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Update a symmetric matrix using an appropriate updating formula.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _UpdateTolerance 1.0e-08
void SymmetricMatrix_Update ( SymmetricMatrix *self, const RealArray1D *dx, const RealArray1D *dg, const SymmetricMatrixUpdating_Option option, const Real *tolerance, Status *status )
{
    if ( ( self != NULL ) && ( dx != NULL ) && ( dg != NULL ) && Status_IsOK ( status ) )
    {
        if ( ( self->extent == View1D_Extent ( dx ) ) && ( self->extent == View1D_Extent ( dg ) ) )
        {
            auto Real         aa, aaxx, ax, fac, gx, mix, tol, xx ;
	    auto RealArray1D *a = NULL ;
            /* . Get the tolerance. */
            if ( tolerance == NULL ) tol = _UpdateTolerance ;
            else                     tol = (*tolerance)     ;
            /* . Allocate space. */
	    a = RealArray1D_AllocateWithExtent ( self->extent, status ) ;
            if ( a != NULL )
            {
                /* . Get A * dx. */
                cblas_dspmv ( CblasColMajor, CblasUpper, self->extent, 1.0e+00, self->data, Array1D_Data ( dx ), dx->stride, 0.0e+00, Array1D_Data ( a ), a->stride ) ;
                /* . For non-BFGS options get g - A * dx. */
                if ( option != SymmetricMatrixUpdating_BFGS )
                {
                    RealArray1D_Scale ( a, -1.0e+00 ) ;
                    RealArray1D_Add   ( a,  1.0e+00, dg, NULL ) ;
                }
                /* . Calculate some dot products. */
                ax = RealArray1D_Dot ( a, dx, NULL ) ;
                /* . BFGS option. */
                if ( option == SymmetricMatrixUpdating_BFGS )
                {
                    gx = RealArray1D_Dot ( dg, dx, NULL ) ;
	            if ( fabs ( ax ) > tol ) cblas_dspr ( CblasColMajor, CblasUpper, self->extent, ( - 1.0e+00 / ax ), Array1D_Data ( a  ), a->stride , self->data ) ;
	            if ( fabs ( gx ) > tol ) cblas_dspr ( CblasColMajor, CblasUpper, self->extent, (   1.0e+00 / gx ), Array1D_Data ( dg ), dg->stride, self->data ) ;
                }
                /* . Other options. */
                else
                {
                    /* . Calculate some dot products. */
                    aa = RealArray1D_Dot (  a,  a, NULL ) ;
                    xx = RealArray1D_Dot ( dx, dx, NULL ) ;
                    aaxx = aa * xx ;
                    /* . Bofill or MS. */
                    if ( ( option == SymmetricMatrixUpdating_Bofill ) || ( option == SymmetricMatrixUpdating_MS ) )
                    {
                        if ( option == SymmetricMatrixUpdating_Bofill )
                        {
                            if ( fabs ( aaxx ) > tol ) fac = ax / aaxx ;
                            else                       fac = 0.0e+00   ;
                        }
                        else
                        {
                            if ( fabs ( ax  ) > tol ) fac = 1.0e+00 / ax ;
                            else                      fac = 0.0e+00      ;
                        }
                        cblas_dspr ( CblasColMajor, CblasUpper, self->extent, fac, Array1D_Data ( a ), a->stride, self->data ) ;
                    }
                    /* . Bofill or Powell. */
                    if ( ( option == SymmetricMatrixUpdating_Bofill ) || ( option == SymmetricMatrixUpdating_Powell ) )
                    {
                        /* . Get the mixing factor. */
                        if ( ( option == SymmetricMatrixUpdating_Bofill ) && ( fabs ( aaxx ) > tol ) ) mix = 1.0e+00 - ax * ax / aaxx ;
                        else                                                                           mix = 1.0e+00                  ;
                        if ( fabs ( xx * xx ) > tol ) cblas_dspr  ( CblasColMajor, CblasUpper, self->extent, ( - mix * ax / ( xx * xx ) ),
                                                                                                       Array1D_Data ( dx ), dx->stride, self->data ) ;
                        if ( fabs ( xx      ) > tol ) cblas_dspr2 ( CblasColMajor, CblasUpper, self->extent, (   mix     /    xx        ),
                                                                    Array1D_Data ( a ), a->stride, Array1D_Data ( dx ), dx->stride, self->data ) ;
                    }
                }
                /* . Finish up. */
	        RealArray1D_Deallocate ( &a ) ;
            }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}
# undef _UpdateTolerance

/*----------------------------------------------------------------------------------------------------------------------------------
! . Multiply by a vector.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetricMatrix_VectorMultiply ( const SymmetricMatrix *self, const RealArray1D *other, RealArray1D *result, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) && ( result != NULL ) && Status_IsOK ( status ) )
    {
        if ( ( self->extent == View1D_Extent ( other ) ) && ( self->extent == View1D_Extent ( other ) ) )
        {
            cblas_dspmv ( CblasColMajor, CblasUpper, self->extent, 1.0e+00, self->data, Array1D_Data ( other  ), other->stride  ,
                                                                   0.0e+00,             Array1D_Data ( result ), result->stride ) ;
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . The view size.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer SymmetricMatrix_ViewSize ( const Integer extent )
{
    Integer size = 0 ;
    if ( extent > 0 ) size = ( extent * ( extent + 1 ) ) / 2 ;
    return size ;
}

# undef _UseCBLAS
# undef _UseReal
