/*==================================================================================================================================
! . 2-D boolean arrays.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "NumericalMacros.h"
# include "RealArray2D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Basic definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _ArrayDataFormat          "%20.10f"
# define _ArrayDataPerLine         6
# define _ArrayDataType            Real
# define _ArrayDataTypeInitializer 0.0e+00
# define _UseCBLAS
# define _UseReal
# include "Array2D_Body.i"
# undef _ArrayDataFormat
# undef _ArrayDataPerLine
# undef _ArrayDataType
# undef _ArrayDataTypeInitializer
# undef _UseCBLAS
# undef _UseReal

/*----------------------------------------------------------------------------------------------------------------------------------
! . Specific functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
# include "BooleanArray1D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Matrix multiply - diagonal values only.
! . Diagonal must be initialized before entry.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RealArray2D_DiagonalOfProduct ( const RealArray2D  *self       ,
                                     const Boolean       sTranspose ,
                                     const RealArray2D  *other      ,
                                     const Boolean       oTranspose ,
                                           RealArray1D  *diagonal   ,
                                           Status       *status     )
{
    if ( ( self     != NULL ) &&
         ( other    != NULL ) &&
         ( diagonal != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Boolean isOK ;
        auto Integer a, b ;
        if ( sTranspose ) { a = View2D_Columns ( self ) ; b = View2D_Rows    ( self ) ; }
        else              { a = View2D_Rows    ( self ) ; b = View2D_Columns ( self ) ; }
        isOK = ( a == View1D_Extent ( diagonal ) ) ;
        if ( oTranspose ) isOK = isOK && ( ( a == View2D_Rows    ( other ) ) && ( b == View2D_Columns ( other ) ) ) ;
        else              isOK = isOK && ( ( a == View2D_Columns ( other ) ) && ( b == View2D_Rows    ( other ) ) ) ;
        if ( isOK )
        {
            auto Integer     i ;
            auto RealArray1D oView, sView ;
            for ( i = 0 ; i < a ; i++ )
            {
                if ( sTranspose ) RealArray2D_ColumnView ( self , i, False, &sView, NULL ) ;
                else              RealArray2D_RowView    ( self , i, False, &sView, NULL ) ;
                if ( oTranspose ) RealArray2D_RowView    ( other, i, False, &oView, NULL ) ;
                else              RealArray2D_ColumnView ( other, i, False, &oView, NULL ) ;
                Array1D_Item ( diagonal, i ) += RealArray1D_Dot ( &sView, &oView, NULL ) ;
            }
       }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Gram-Schmidt orthogonalize a set of vectors stored columnwise in-place.
! . A modified (as opposed to classical) iterative algorithm is employed.
! . The number of orthogonal vectors is returned (<= old number).
! . There is an option to treat the first numberConstant vectors as already orthogonalized.
! . The tolerance is the size of norm2 per element.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define DEFAULT_TOLERANCE 1.0e-10
Integer RealArray2D_GramSchmidtOrthogonalize ( const RealArray2D *self              ,
                                               const Integer     *maximumIterations ,
                                               const Integer     *numberConstant    ,
                                               const Real        *tolerance         ,
                                                     Status      *status            )
{
    Integer numberOrthogonalized = 0 ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer nStart, nVectors ;
        /* . Get various indices. */
        nVectors = View2D_Columns ( self ) ;
        if ( numberConstant != NULL ) nStart = Maximum ( 0, (*numberConstant) ) ;
        else                          nStart = 0 ;
        if ( nStart < nVectors )
        {
            auto Integer     i, iteration, j, nCurrent, nIterations ;
            auto Real        delta, factor ;
            auto RealArray1D iVector, jVector ;
            /* . Get the number of iterations. */
            if ( maximumIterations != NULL ) nIterations = Maximum ( 1, (*maximumIterations) ) ;
            else                             nIterations = 1 ;
            /* . Get the tolerance for normalization. */
            if ( tolerance == NULL ) delta = DEFAULT_TOLERANCE     ;
            else                     delta = fabs ( (*tolerance) ) ;
            delta *= sqrt ( ( Real ) View2D_Rows ( self ) ) ;
            /* . Loop over vectors to be orthogonalized. */
            for ( i = nCurrent = nStart ; i < nVectors ; i++ )
            {
                RealArray2D_ColumnView ( self, i, False, &iVector, status ) ;
                /* . Loop over iterations. */
                for ( iteration = 0 ; iteration < nIterations ; iteration++ )
                {
	            /* . Loop over vectors to be orthogonalized against. */
                    for ( j = 0 ; j < nCurrent ; j++ )
                    {
                        RealArray2D_ColumnView   ( self, j, False, &jVector, status ) ;
                        factor = RealArray1D_Dot ( &iVector, &jVector, status ) ;
                        RealArray1D_Add          ( &iVector, -factor, &jVector, status ) ;
                    }
                }
                /* . Normalization. If OK, scale and move the vector if necessary. */
                factor = RealArray1D_Norm2 ( &iVector ) ;
                if ( factor > delta )
                {
                    RealArray1D_Scale ( &iVector, 1.0e+00 / factor ) ;
                    if ( nCurrent != i )
                    {
                        RealArray2D_ColumnView ( self, nCurrent, False, &jVector, status ) ;
                        RealArray1D_CopyTo     ( &iVector, &jVector, status ) ;
                    }
                    nCurrent++ ;
                    numberOrthogonalized++ ;
                }
            }
        }
    }
    return numberOrthogonalized ;
}
# undef DEFAULT_TOLERANCE

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for a diagonal matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean RealArray2D_IsDiagonal ( const RealArray2D *self, const Real tolerance )
{
    Boolean isDiagonal = False ;
    if ( self != NULL )
    {
        auto Integer i, j ;
        isDiagonal = True ;
        for ( i = 0 ; i < self->extent0 ; i++ )
        {
            for ( j = 0 ; j < self->extent1 ; j++ )
            {
                if ( ( i != j ) && ( fabs ( Array2D_Item ( self, i, j ) ) > tolerance ) ) { isDiagonal = False ; break ; }
            }
            if ( ! isDiagonal ) break ;
        }
    }
    return isDiagonal ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for orthogonality (check procedure).
!---------------------------------------------------------------------------------------------------------------------------------*/
# define ORTHOGONALITY_TOLERANCE 1.0e-10
Boolean RealArray2D_IsOrthogonal ( const RealArray2D *self, const Real *tolerance, Real *deviation, Status *status )
{
    Boolean isOrthogonal = False ;
    if ( deviation != NULL ) (*deviation) = 0.0e+00 ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer      i ;
        auto Real         absmax, tol ;
        auto RealArray2D *r ;
        /* . Allocate space. */
        r = RealArray2D_AllocateWithExtents ( View2D_Columns ( self ), View2D_Columns ( self ), status ) ;
        if ( r != NULL )
        {
            /* . Get the tolerance. */
            if ( tolerance == NULL ) tol = ORTHOGONALITY_TOLERANCE ;
            else                     tol = (*tolerance) ;
            /* . Get self^T * self - I. */
            RealArray2D_MatrixMultiply ( True, False, 1.0e+00, self, self, 0.0e+00, r, status ) ;
            for ( i = 0 ; i < View2D_Rows ( r ) ; i++ ) Array2D_Item ( r, i, i ) -= 1.0e+00 ;
            /* . Find maximum deviation. */
            absmax       = RealArray2D_AbsoluteMaximum ( r ) ;
            isOrthogonal = ( absmax <= tol ) ;
            if ( deviation != NULL ) (*deviation) = absmax ;
        }
        /* . Finish up. */
        RealArray2D_Deallocate ( &r ) ;
    }
    return isOrthogonal ;
}
# undef ORTHOGONALITY_TOLERANCE

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for a symmetric matrix (check procedure).
!---------------------------------------------------------------------------------------------------------------------------------*/
# define SYMMETRIC_TOLERANCE 1.0e-10
Boolean RealArray2D_IsSymmetric ( const RealArray2D *self, const Real *tolerance, Real *deviation )
{
    Boolean isSymmetric = False ;
    if ( deviation != NULL ) (*deviation) = 0.0e+00 ;
    if ( ( self != NULL ) && View2D_IsSquare ( self ) )
    {
        auto Integer i, j ;
        auto Real    difference, tol ;
        /* . Get the tolerance. */
        if ( tolerance == NULL ) tol = SYMMETRIC_TOLERANCE ;
        else                     tol = (*tolerance) ;
        /* . Loop over indices. */
        difference = 0.0e+00 ;
        for ( i = 0 ; i < View2D_Rows ( self ) ; i++ )
        {
            for ( j = 0 ; j < i ; j++ ) difference = Maximum ( difference, fabs ( Array2D_Item ( self, i, j ) - Array2D_Item ( self, j, i ) ) ) ;
        }
        /* . Find maximum difference. */
        isSymmetric = ( difference <= tol ) ;
        if ( deviation != NULL ) (*deviation) = difference ;
    }
    return isSymmetric ;
}
# undef SYMMETRIC_TOLERANCE

/*----------------------------------------------------------------------------------------------------------------------------------
! . Matrix-matrix multiply. A straight interface to cblas_dgemm. C = alpha A B + beta C (A and B can be transposed).
!---------------------------------------------------------------------------------------------------------------------------------*/
void RealArray2D_MatrixMultiply ( const Boolean      aTranspose ,
                                  const Boolean      bTranspose ,
                                  const Real         alpha      ,
                                  const RealArray2D *a          ,
                                  const RealArray2D *b          ,
                                  const Real         beta       ,
                                        RealArray2D *c          ,
                                        Status      *status     )
{
    if ( ( a != NULL ) && ( b != NULL ) && ( c != NULL ) && Status_IsOK ( status ) )
    {
        /* . All matrices need to be compact in dimension 1. */
        if ( View2D_IsCompact1 ( a ) && View2D_IsCompact1 ( b ) && View2D_IsCompact1 ( c ) )
        {
            auto Integer  k, m, n ;
            m = View2D_Rows    ( c ) ;
            n = View2D_Columns ( c ) ;
	    if ( aTranspose ) k = View2D_Rows    ( a ) ;
	    else              k = View2D_Columns ( a ) ;
            if ( ( ( ! aTranspose ) && ( ! bTranspose ) && ( m == View2D_Rows    ( a ) ) && ( n == View2D_Columns ( b ) ) && ( View2D_Columns ( a ) == View2D_Rows    ( b ) ) ) ||
                 ( ( ! aTranspose ) && (   bTranspose ) && ( m == View2D_Rows    ( a ) ) && ( n == View2D_Rows    ( b ) ) && ( View2D_Columns ( a ) == View2D_Columns ( b ) ) ) ||
                 ( (   aTranspose ) && ( ! bTranspose ) && ( m == View2D_Columns ( a ) ) && ( n == View2D_Columns ( b ) ) && ( View2D_Rows    ( a ) == View2D_Rows    ( b ) ) ) ||
                 ( (   aTranspose ) && (   bTranspose ) && ( m == View2D_Columns ( a ) ) && ( n == View2D_Rows    ( b ) ) && ( View2D_Rows    ( a ) == View2D_Columns ( b ) ) ) )
	    {
                auto enum CBLAS_TRANSPOSE aT, bT ;
                if ( aTranspose ) aT = CblasTrans   ;
                else              aT = CblasNoTrans ;
                if ( bTranspose ) bT = CblasTrans   ;
                else              bT = CblasNoTrans ;
                cblas_dgemm ( CblasRowMajor, aT, bT, m, n, k, alpha, Array_DataPointer ( a ), a->stride0 ,
                                                                     Array_DataPointer ( b ), b->stride0 ,
                                                              beta , Array_DataPointer ( c ), c->stride0 ) ;
            }
            else Status_Set ( status, Status_NonConformableArrays ) ;
        }
        else Status_Set ( status, Status_InvalidArrayOperation ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Project a matrix from a vector.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RealArray2D_ProjectOutOfArray1D ( const RealArray2D *self, RealArray1D *vector, Status *status )
{
    if ( ( self != NULL ) && ( vector != NULL ) && Status_IsOK ( status ) )
    {
        auto RealArray1D *Pv ;
        Pv = RealArray1D_AllocateWithExtent ( View2D_Columns ( self ), status ) ;
        RealArray2D_VectorMultiply ( True ,  1.0e+00, self,  vector, 0.0e+00,     Pv, status ) ;
        RealArray2D_VectorMultiply ( False, -1.0e+00, self,  Pv    , 1.0e+00, vector, status ) ;
        RealArray1D_Deallocate ( &Pv ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Trace.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real RealArray2D_Trace ( const RealArray2D *self, Status *status )
{
    Real trace = 0.0e+00 ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        if ( View2D_IsSquare ( self ) )
        {
            auto Integer i ;
            for ( i = 0 ; i < View2D_Rows ( self ) ; i++ ) trace += Array2D_Item ( self, i, i ) ;
        }
        else Status_Set ( status, Status_InvalidArrayOperation ) ;
    }
    return trace ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Trace of product.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real RealArray2D_TraceOfProduct ( const RealArray2D *self, const RealArray2D *other, Status *status )
{
    Real trace = 0.0e+00 ;
    if ( ( self != NULL ) && ( other != NULL ) && Status_IsOK ( status ) )
    {
        if ( ( View2D_Rows ( self ) == View2D_Columns ( other ) ) && ( View2D_Columns ( self ) == View2D_Rows ( other ) ) )
        {
            auto Integer i ;
            for ( i = 0 ; i < View2D_Rows ( self ) ; i++ ) trace += cblas_ddot ( View2D_Columns ( self ), Array2D_RowPointer    ( self , i ), self->stride1  ,
                                                                                                          Array2D_ColumnPointer ( other, i ), other->stride0 ) ;
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
    return trace ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Create a transposed clone of the matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
RealArray2D *RealArray2D_TransposeClone ( const RealArray2D *self, Status *status )
{
    RealArray2D *clone = NULL ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer c, r ;
        c     = View2D_Columns ( self ) ;
        r     = View2D_Rows    ( self ) ;
        clone = RealArray2D_AllocateWithExtents ( c, r, status ) ;
        if ( clone != NULL )
        {
            auto Integer     i ;
            auto RealArray1D column, row ;
            for ( i = 0 ; i < r ; i++ )
            {
                RealArray2D_ColumnView ( clone, i, False, &column, status ) ;
                RealArray2D_RowView    ( self , i, False, &row   , status ) ;
                RealArray1D_CopyTo ( &row, &column, status ) ;
            }
        }
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Transpose a general matrix in-place (the matrices must be uniform).
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The function pi. */
# define IndexFunction( pIn, r, c, pOut ) { t = pIn / r ; pOut = c * ( pIn - r * t ) + t ; }

void RealArray2D_TransposeGeneral ( RealArray2D *self, Status *status )
{
    if ( ( self != NULL ) && ( View2D_Size ( self ) > 1 ) && Status_IsOK ( status ) )
    {
        if ( View2D_IsUniform ( self ) )
        {
            /*
            !
            ! . From E. G. Cate & D. W. Twigg, "Algorithm 513", ACM Transactions in Math. Software, 3, 104-110, 1977.
            !
            ! . Indexing functions (without offset and stride):
            !
            !   index[i,j] = ( i * c + j ) -> transposed[j,i] = ( j * r + i )
            !
            !   pi    ( index     [i,j] ) = transposed[j,i] => pi    ( a ) = ( a * r ) mod ( c * r )
            !   pi^-1 ( transposed[j,i] ) = index     [i,j] => pi^-1 ( a ) = ( a * c )
            !
            */
            BooleanArray1D *isMoved ;
            Integer         columns, cycles, iWork, last, maximumP0, moved, p0, p1, p2, q0, q1, q2, rows, size, stride, t ;
            Real            B, C ;

            /* . Initialization. */
            columns = View2D_Columns ( self ) ;
            rows    = View2D_Rows    ( self ) ;
            stride  = self->stride1 ;

            /* . Try and allocate space - not necessary but increases efficiency. */
            iWork   = ( columns + rows + 1 ) / 2 ;
            isMoved = BooleanArray1D_AllocateWithExtent ( iWork, NULL ) ;
            if ( isMoved == NULL ) iWork = 0 ;
            else                   BooleanArray1D_Set ( isMoved, False ) ;

            /* . Find the number of items (fixed points) which do not need to be moved. */
            moved = 2 ;
            if ( ( rows > 2 ) && ( columns > 2 ) ) moved += ( Integer_GCD ( rows - 1, columns - 1 ) - 1 ) ;

            /* . Initial values for the search. */
            size = rows * columns ;
            last = size - 1 ;
            p0   = 1 ;

            /* . Loop until all items moved. */
            cycles = 0 ;
            while ( moved < size )
            {
                /* . Find the start of the next cycle. */
                if ( cycles > 0 )
                {
                    while ( True )
                    {
                        maximumP0 = last - p0 ;
                        p0       += 1 ;
                        IndexFunction ( p0, rows, columns, p2 ) ;
                        if ( p0 > maximumP0 ) { Status_Set ( status, Status_AlgorithmError ) ; break ; }
                        else if ( p0 != p2 )
                        {
                            if ( p0 >= iWork )
                            {
                                while ( ( p2 > p0 ) && ( p2 < maximumP0 ) )
                                {
                                    p1 = p2 ;
                                    IndexFunction ( p1, rows, columns, p2 ) ;
                                }
                                if ( p2 == p0 ) break ;
                            }
                            else if ( ! Array1D_Item ( isMoved, p0 ) ) break ;
                        }
                    }
                }

                /* . Rearrange the items of a loop and its companion loop. */
                p1 = p0 ;
                q0 = last - p0 ;
                q1 = q0 ;
                B  = Array2D_ItemByIndex ( self, p1 * stride ) ;
                C  = Array2D_ItemByIndex ( self, q1 * stride ) ;
                while ( True )
                {
                    IndexFunction ( p1, rows, columns, p2 ) ;
                    q2 = last - p2 ;
                    if ( p1 < iWork ) Array1D_Item ( isMoved, p1 ) = True ;
                    if ( q1 < iWork ) Array1D_Item ( isMoved, q1 ) = True ;
                    moved += 2 ;
                    if ( p2 == p0 )
                    {
                        Array2D_ItemByIndex ( self, p1 * stride ) = B ;
                        Array2D_ItemByIndex ( self, q1 * stride ) = C ;
                        break ;
                    }
                    if ( p2 == q0 )
                    {
                        Array2D_ItemByIndex ( self, p1 * stride ) = C ;
                        Array2D_ItemByIndex ( self, q1 * stride ) = B ;
                        break ;
                    }
                    Array2D_ItemByIndex ( self, p1 * stride ) = Array2D_ItemByIndex ( self, p2 * stride ) ;
                    Array2D_ItemByIndex ( self, q1 * stride ) = Array2D_ItemByIndex ( self, q2 * stride ) ;
                    p1 = p2 ;
                    q1 = q2 ;
                }

                /* . Finish up. */
                cycles += 1 ;
                /*printf ( "Loop information ( cycle, tag, start, moved, size ): %5d %5d %5d %5d", cycles, p0, moved, size ) ; */
            }

            /* . Finish up. */
            BooleanArray1D_Deallocate ( &isMoved ) ;

            /* . Reset the view variables. */
            t             = self->extent0 ;
            self->extent0 = self->extent1 ;
            self->extent1 = t             ;
            self->stride0 = t * stride    ;
            self->stride1 = stride        ;
        }
        else Status_Set ( status, Status_InvalidArrayOperation ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Transpose a square matrix in-place.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RealArray2D_TransposeSquare ( RealArray2D *self, Status *status )
{
    if ( ( self != NULL ) && ( View2D_Size ( self ) > 1 ) && Status_IsOK ( status ) )
    {
        /* . Square arrays. */
        if ( View2D_IsSquare ( self ) )
        {
            auto Integer i, ij, j, ji ;
            auto Real    t ;
            /* . Loop over rows and columns. */
            for ( i = 0 ; i < View2D_Rows ( self ) ; i++ )
            {
                for ( j = 0 ; j < i ; j++ )
                {
                    ij = View2D_ItemIndex ( self, i, j ) ;
                    ji = View2D_ItemIndex ( self, j, i ) ;
                    t  = Array2D_ItemByIndex ( self, ij ) ;
                    Array2D_ItemByIndex ( self, ij ) = Array2D_ItemByIndex ( self, ji ) ;
                    Array2D_ItemByIndex ( self, ji ) = t ;
                }
            }
        }
        else Status_Set ( status, Status_InvalidArrayOperation ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Matrix-vector multiply. A straight interface to cblas_dgemv. y = alpha A x + beta y (A can be transposed).
!---------------------------------------------------------------------------------------------------------------------------------*/
void RealArray2D_VectorMultiply ( const Boolean      aTranspose ,
                                  const Real         alpha      ,
                                  const RealArray2D *a          ,
                                  const RealArray1D *x          ,
                                  const Real         beta       ,
                                        RealArray1D *y          ,
                                        Status      *status     )
{
    if ( ( a != NULL ) && ( x != NULL ) && ( y != NULL ) && Status_IsOK ( status ) )
    {
        /* . The matrix needs to be compact in dimension 1. */
        if ( ! View2D_IsCompact1 ( a ) ) Status_Set ( status, Status_InvalidArrayOperation ) ;
        else
        {
            auto Integer m, n ;
            m = View2D_Rows    ( a ) ;
            n = View2D_Columns ( a ) ;
            if ( ( ( ! aTranspose ) && ( n == View1D_Extent ( x ) ) && ( m == View1D_Extent ( y ) ) ) ||
                 ( (   aTranspose ) && ( m == View1D_Extent ( x ) ) && ( n == View1D_Extent ( y ) ) ) )
            {
                auto enum CBLAS_TRANSPOSE aT ;
                if ( aTranspose ) aT = CblasTrans   ;
                else              aT = CblasNoTrans ;
                cblas_dgemv ( CblasRowMajor, aT, m, n, alpha, Array_DataPointer ( a ), a->stride0 ,
                                                              Array_DataPointer ( x ), x->stride  ,
                                                       beta , Array_DataPointer ( y ), y->stride  ) ;
            }
            else Status_Set ( status, Status_NonConformableArrays ) ;
        }
    }
}
