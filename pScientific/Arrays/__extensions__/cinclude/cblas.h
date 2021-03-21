# ifndef _CBLAS
# define _CBLAS

/* . Standard cblas headers - adapted to pDynamo. */

# include <stddef.h>
# include "Integer.h"
# include "Real.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define CBLAS_INDEX size_t

enum CBLAS_ORDER     { CblasRowMajor = 101 , CblasColMajor = 102 } ;
enum CBLAS_TRANSPOSE { CblasNoTrans  = 111 , CblasTrans    = 112 , CblasConjTrans = 113 } ;
enum CBLAS_UPLO	     { CblasUpper    = 121 , CblasLower    = 122 } ;
enum CBLAS_DIAG	     { CblasNonUnit  = 131 , CblasUnit     = 132 } ;
enum CBLAS_SIDE	     { CblasLeft     = 141 , CblasRight    = 142 } ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Level 1 BLAS functions and procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real        cblas_dasum  ( const Integer N, const Real *X, const Integer incX ) ;
void        cblas_daxpy  ( const Integer N, const Real alpha, const Real *X, const Integer incX, Real *Y, const Integer incY ) ;
void        cblas_dcopy  ( const Integer N, const Real *X, const Integer incX, Real *Y, const Integer incY ) ;
Real        cblas_ddot   ( const Integer N, const Real *X, const Integer incX, const Real *Y, const Integer incY ) ;
Real        cblas_dnrm2  ( const Integer N, const Real *X, const Integer incX ) ;
void        cblas_drot   ( const Integer N, Real *X, const Integer incX, Real *Y, const Integer incY, const Real c, const Real s ) ;
void        cblas_drotg  ( Real *a, Real *b, Real *c, Real *s ) ;
void        cblas_drotm  ( const Integer N, Real *X, const Integer incX, Real *Y, const Integer incY, const Real *P ) ;
void        cblas_drotmg ( Real *d1, Real *d2, Real *b1, const Real b2, Real *P ) ;
void        cblas_dscal  ( const Integer N, const Real alpha, Real *X, const Integer incX ) ;
void        cblas_dswap  ( const Integer N, Real *X, const Integer incX, Real *Y, const Integer incY ) ;
CBLAS_INDEX cblas_idamax ( const Integer N, const Real *X, const Integer incX ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Level 2 BLAS functions and procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
void cblas_dgbmv ( const enum CBLAS_ORDER Order,
                   const enum CBLAS_TRANSPOSE TransA, const Integer M, const Integer N,
                   const Integer KL, const Integer KU, const Real alpha,
                   const Real *A, const Integer lda, const Real *X,
                   const Integer incX, const Real beta, Real *Y, const Integer incY ) ;
void cblas_dgemv ( const enum CBLAS_ORDER Order,
                   const enum CBLAS_TRANSPOSE TransA, const Integer M, const Integer N,
                   const Real alpha, const Real *A, const Integer lda,
                   const Real *X, const Integer incX, const Real beta,
                   Real *Y, const Integer incY ) ;
void cblas_dger  ( const enum CBLAS_ORDER Order, const Integer M, const Integer N,
                   const Real alpha, const Real *X, const Integer incX,
                   const Real *Y, const Integer incY, Real *A, const Integer lda ) ;
void cblas_dsbmv ( const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const Integer N, const Integer K, const Real alpha, const Real *A,
                   const Integer lda, const Real *X, const Integer incX,
                   const Real beta, Real *Y, const Integer incY ) ;
void cblas_dspmv ( const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const Integer N, const Real alpha, const Real *Ap,
                   const Real *X, const Integer incX,
                   const Real beta, Real *Y, const Integer incY ) ;
void cblas_dspr  ( const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const Integer N, const Real alpha, const Real *X,
                   const Integer incX, Real *Ap ) ;
void cblas_dspr2 ( const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const Integer N, const Real alpha, const Real *X,
                   const Integer incX, const Real *Y, const Integer incY, Real *A ) ;
void cblas_dsymv ( const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const Integer N, const Real alpha, const Real *A,
                   const Integer lda, const Real *X, const Integer incX,
                   const Real beta, Real *Y, const Integer incY ) ;
void cblas_dsyr  ( const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const Integer N, const Real alpha, const Real *X,
                   const Integer incX, Real *A, const Integer lda ) ;
void cblas_dsyr2 ( const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const Integer N, const Real alpha, const Real *X,
                   const Integer incX, const Real *Y, const Integer incY, Real *A,
                   const Integer lda ) ;
void cblas_dtbmv ( const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                   const Integer N, const Integer K, const Real *A, const Integer lda,
                   Real *X, const Integer incX ) ;
void cblas_dtbsv ( const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                   const Integer N, const Integer K, const Real *A, const Integer lda,
                   Real *X, const Integer incX ) ;
void cblas_dtpmv ( const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                   const Integer N, const Real *Ap, Real *X, const Integer incX ) ;
void cblas_dtpsv ( const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                   const Integer N, const Real *Ap, Real *X, const Integer incX ) ;
void cblas_dtrmv ( const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                   const Integer N, const Real *A, const Integer lda,
                   Real *X, const Integer incX ) ;
void cblas_dtrsv ( const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                   const Integer N, const Real *A, const Integer lda, Real *X,
                   const Integer incX ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Level 3 BLAS functions and procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
void cblas_dgemm ( const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                   const enum CBLAS_TRANSPOSE TransB, const Integer M, const Integer N,
                   const Integer K, const Real alpha, const Real *A,
                   const Integer lda, const Real *B, const Integer ldb,
                   const Real beta, Real *C, const Integer ldc ) ;
void cblas_dsymm ( const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                   const enum CBLAS_UPLO Uplo, const Integer M, const Integer N,
                   const Real alpha, const Real *A, const Integer lda,
                   const Real *B, const Integer ldb, const Real beta,
                  Real *C, const Integer ldc ) ;
void cblas_dsyrk ( const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const enum CBLAS_TRANSPOSE Trans, const Integer N, const Integer K,
                   const Real alpha, const Real *A, const Integer lda,
                   const Real beta, Real *C, const Integer ldc ) ;
void cblas_dsyr2k ( const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                    const enum CBLAS_TRANSPOSE Trans, const Integer N, const Integer K,
                    const Real alpha, const Real *A, const Integer lda,
                    const Real *B, const Integer ldb, const Real beta,
                   Real *C, const Integer ldc ) ;
void cblas_dtrmm ( const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                   const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                   const enum CBLAS_DIAG Diag, const Integer M, const Integer N,
                   const Real alpha, const Real *A, const Integer lda,
                   Real *B, const Integer ldb ) ;
void cblas_dtrsm ( const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                   const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                   const enum CBLAS_DIAG Diag, const Integer M, const Integer N,
                   const Real alpha, const Real *A, const Integer lda,
                   Real *B, const Integer ldb ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Utility functions and procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
void cblas_xerbla ( Integer p, const char *rout, const char *form, ... ) ;

# endif
