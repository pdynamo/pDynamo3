# ifndef _F2CBLAS
# define _F2CBLAS

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define FCHAR char *
# define FINT  const int *

/* . Level 1 BLAS. */
# define F77_drotg      drotg_
# define F77_drotmg     drotmg_
# define F77_drot       drot_
# define F77_drotm      drotm_
# define F77_dswap      dswap_
# define F77_dcopy      dcopy_
# define F77_daxpy      daxpy_
# define F77_idamax_sub idamax_
# define F77_ddot_sub   ddot_
# define F77_dscal      dscal_
# define F77_dnrm2_sub  dnrm2_
# define F77_dasum_sub  dasum_

/* . Level 2 BLAS. */
# define F77_dsymv      dsymv_
# define F77_dsbmv      dsbmv_
# define F77_dspmv      dspmv_
# define F77_dger       dger_
# define F77_dsyr       dsyr_
# define F77_dspr       dspr_
# define F77_dsyr2      dsyr2_
# define F77_dspr2      dspr2_
# define F77_dgemv      dgemv_
# define F77_dgbmv      dgbmv_
# define F77_dtrmv      dtrmv_
# define F77_dtbmv      dtbmv_
# define F77_dtpmv      dtpmv_
# define F77_dtrsv      dtrsv_
# define F77_dtbsv      dtbsv_
# define F77_dtpsv      dtpsv_

/* . Level 3 BLAS. */
# define F77_dgemm      dgemm_
# define F77_dsymm      dsymm_
# define F77_dsyrk      dsyrk_
# define F77_dsyr2k     dsyr2k_
# define F77_dtrmm      dtrmm_
# define F77_dtrsm      dtrsm_

/* . Utility functions. */
# define F77_xerbla     xerbla_

/*----------------------------------------------------------------------------------------------------------------------------------
! . Prototypes.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Level 1 BLAS. */
void F77_drot(FINT, double *, FINT, double *, FINT, const double *, const double *);
void F77_drotg(double *,double *,double *,double *);
void F77_drotm( FINT, double *, FINT, double *, FINT, const double *);
void F77_drotmg(double *,double *,double *,const double *, double *);
void F77_dswap( FINT, double *, FINT, double *, FINT);
void F77_dcopy( FINT, const double *, FINT, double *, FINT);
void F77_daxpy( FINT, const double *, const double *, FINT, double *, FINT);
void F77_dswap( FINT, double *, FINT, double *, FINT);
double F77_ddot_sub( FINT, const double *, FINT, const double *, FINT );
void F77_dscal( FINT, const double *, double *, FINT);
double F77_dnrm2_sub( FINT, const double *, FINT );
void F77_dasum_sub( FINT, const double *, FINT, double *);
int F77_idamax_sub( FINT, const double * , FINT );

/* . Level 2 BLAS. */
void F77_dgemv(FCHAR, FINT, FINT, const double *, const double *, FINT, const double *, FINT, const double *, double *, FINT);
void F77_dgbmv(FCHAR, FINT, FINT, FINT, FINT, const double *,  const double *, FINT, const double *, FINT, const double *, double *, FINT);
void F77_dsymv(FCHAR, FINT, const double *, const double *, FINT, const double *,  FINT, const double *, double *, FINT);
void F77_dsbmv(FCHAR, FINT, FINT, const double *, const double *, FINT, const double *, FINT, const double *, double *, FINT);
void F77_dspmv(FCHAR, FINT, const double *, const double *, const double *, FINT, const double *, double *, FINT);
void F77_dtrmv( FCHAR, FCHAR, FCHAR, FINT, const double *, FINT, double *, FINT);
void F77_dtbmv( FCHAR, FCHAR, FCHAR, FINT, FINT, const double *, FINT, double *, FINT);
void F77_dtrsv( FCHAR, FCHAR, FCHAR, FINT, const double *, FINT, double *, FINT);
void F77_dtbsv( FCHAR, FCHAR, FCHAR, FINT, FINT, const double *, FINT, double *, FINT);
void F77_dtpmv( FCHAR, FCHAR, FCHAR, FINT, const double *, double *, FINT);
void F77_dtpsv( FCHAR, FCHAR, FCHAR, FINT, const double *, double *, FINT);
void F77_dger( FINT, FINT, const double *, const double *, FINT, const double *, FINT, double *, FINT);
void F77_dsyr(FCHAR, FINT, const double *, const double *, FINT, double *, FINT);
void F77_dspr(FCHAR, FINT, const double *, const double *, FINT, double *);
void F77_dspr2(FCHAR, FINT, const double *, const double *, FINT, const double *, FINT,  double *);
void F77_dsyr2(FCHAR, FINT, const double *, const double *, FINT, const double *, FINT,  double *, FINT);

/* . Level 3 BLAS. */
void F77_dgemm(FCHAR, FCHAR, FINT, FINT, FINT, const double *, const double *, FINT, const double *, FINT, const double *, double *, FINT);
void F77_dsymm(FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, const double *, FINT, const double *, double *, FINT);
void F77_dsyrk(FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, const double *, double *, FINT);
void F77_dsyr2k(FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, const double *, FINT, const double *, double *, FINT);
void F77_dtrmm(FCHAR, FCHAR, FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, double *, FINT);
void F77_dtrsm(FCHAR, FCHAR, FCHAR, FCHAR, FINT, FINT, const double *, const double *, FINT, double *, FINT);

/* . Utility functions. */
void F77_xerbla(FCHAR, void *);

# endif
