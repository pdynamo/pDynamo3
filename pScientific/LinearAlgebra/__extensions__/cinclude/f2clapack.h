# ifndef _F2CLAPACK
# define _F2CLAPACK

/* . All double precision procedures are currently listed. */

# include "f2c.h"

/* . For machine constants. */
/* Function */ extern doublereal dlamch_ ( char * );

/* Subroutine */ extern int dbdsdc_(char *uplo, char *compq, integer *n, doublereal *
	d__, doublereal *e, doublereal *u, integer *ldu, doublereal *vt,
	integer *ldvt, doublereal *q, integer *iq, doublereal *work, integer *
	iwork, integer *info);

/* Subroutine */ extern int dbdsqr_(char *uplo, integer *n, integer *ncvt, integer *
	nru, integer *ncc, doublereal *d__, doublereal *e, doublereal *vt,
	integer *ldvt, doublereal *u, integer *ldu, doublereal *c__, integer *
	ldc, doublereal *work, integer *info);

/* Subroutine */ extern int ddisna_(char *job, integer *m, integer *n, doublereal *
	d__, doublereal *sep, integer *info);

/* Subroutine */ extern int dgbbrd_(char *vect, integer *m, integer *n, integer *ncc,
	 integer *kl, integer *ku, doublereal *ab, integer *ldab, doublereal *
	d__, doublereal *e, doublereal *q, integer *ldq, doublereal *pt,
	integer *ldpt, doublereal *c__, integer *ldc, doublereal *work,
	integer *info);

/* Subroutine */ extern int dgbcon_(char *norm, integer *n, integer *kl, integer *ku,
	 doublereal *ab, integer *ldab, integer *ipiv, doublereal *anorm,
	doublereal *rcond, doublereal *work, integer *iwork, integer *info);

/* Subroutine */ extern int dgbequ_(integer *m, integer *n, integer *kl, integer *ku,
	 doublereal *ab, integer *ldab, doublereal *r__, doublereal *c__,
	doublereal *rowcnd, doublereal *colcnd, doublereal *amax, integer *
	info);

/* Subroutine */ extern int dgbrfs_(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, doublereal *ab, integer *ldab, doublereal *afb,
	integer *ldafb, integer *ipiv, doublereal *b, integer *ldb,
	doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr,
	doublereal *work, integer *iwork, integer *info);

/* Subroutine */ extern int dgbsv_(integer *n, integer *kl, integer *ku, integer *
	nrhs, doublereal *ab, integer *ldab, integer *ipiv, doublereal *b,
	integer *ldb, integer *info);

/* Subroutine */ extern int dgbsvx_(char *fact, char *trans, integer *n, integer *kl,
	 integer *ku, integer *nrhs, doublereal *ab, integer *ldab,
	doublereal *afb, integer *ldafb, integer *ipiv, char *equed,
	doublereal *r__, doublereal *c__, doublereal *b, integer *ldb,
	doublereal *x, integer *ldx, doublereal *rcond, doublereal *ferr,
	doublereal *berr, doublereal *work, integer *iwork, integer *info);

/* Subroutine */ extern int dgbtf2_(integer *m, integer *n, integer *kl, integer *ku,
	 doublereal *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ extern int dgbtrf_(integer *m, integer *n, integer *kl, integer *ku,
	 doublereal *ab, integer *ldab, integer *ipiv, integer *info);

/* Subroutine */ extern int dgbtrs_(char *trans, integer *n, integer *kl, integer *
	ku, integer *nrhs, doublereal *ab, integer *ldab, integer *ipiv,
	doublereal *b, integer *ldb, integer *info);

/* Subroutine */ extern int dgebak_(char *job, char *side, integer *n, integer *ilo,
	integer *ihi, doublereal *scale, integer *m, doublereal *v, integer *
	ldv, integer *info);

/* Subroutine */ extern int dgebal_(char *job, integer *n, doublereal *a, integer *
	lda, integer *ilo, integer *ihi, doublereal *scale, integer *info);

/* Subroutine */ extern int dgebd2_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *d__, doublereal *e, doublereal *tauq, doublereal *
	taup, doublereal *work, integer *info);

/* Subroutine */ extern int dgebrd_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *d__, doublereal *e, doublereal *tauq, doublereal *
	taup, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ extern int dgecon_(char *norm, integer *n, doublereal *a, integer *
	lda, doublereal *anorm, doublereal *rcond, doublereal *work, integer *
	iwork, integer *info);

/* Subroutine */ extern int dgeequ_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *r__, doublereal *c__, doublereal *rowcnd, doublereal
	*colcnd, doublereal *amax, integer *info);

/* Subroutine */ extern int dgees_(char *jobvs, char *sort, L_fp select, integer *n,
	doublereal *a, integer *lda, integer *sdim, doublereal *wr,
	doublereal *wi, doublereal *vs, integer *ldvs, doublereal *work,
	integer *lwork, logical *bwork, integer *info);

/* Subroutine */ extern int dgeesx_(char *jobvs, char *sort, L_fp select, char *
	sense, integer *n, doublereal *a, integer *lda, integer *sdim,
	doublereal *wr, doublereal *wi, doublereal *vs, integer *ldvs,
	doublereal *rconde, doublereal *rcondv, doublereal *work, integer *
	lwork, integer *iwork, integer *liwork, logical *bwork, integer *info);

/* Subroutine */ extern int dgeev_(char *jobvl, char *jobvr, integer *n, doublereal *
	a, integer *lda, doublereal *wr, doublereal *wi, doublereal *vl,
	integer *ldvl, doublereal *vr, integer *ldvr, doublereal *work,
	integer *lwork, integer *info);

/* Subroutine */ extern int dgeevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, integer *n, doublereal *a, integer *lda, doublereal *wr,
	doublereal *wi, doublereal *vl, integer *ldvl, doublereal *vr,
	integer *ldvr, integer *ilo, integer *ihi, doublereal *scale,
	doublereal *abnrm, doublereal *rconde, doublereal *rcondv, doublereal
	*work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ extern int dgegs_(char *jobvsl, char *jobvsr, integer *n,
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	alphar, doublereal *alphai, doublereal *beta, doublereal *vsl,
	integer *ldvsl, doublereal *vsr, integer *ldvsr, doublereal *work,
	integer *lwork, integer *info);

/* Subroutine */ extern int dgegv_(char *jobvl, char *jobvr, integer *n, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *alphar,
	doublereal *alphai, doublereal *beta, doublereal *vl, integer *ldvl,
	doublereal *vr, integer *ldvr, doublereal *work, integer *lwork,
	integer *info);

/* Subroutine */ extern int dgehd2_(integer *n, integer *ilo, integer *ihi,
	doublereal *a, integer *lda, doublereal *tau, doublereal *work,
	integer *info);

/* Subroutine */ extern int dgehrd_(integer *n, integer *ilo, integer *ihi,
	doublereal *a, integer *lda, doublereal *tau, doublereal *work,
	integer *lwork, integer *info);

/* Subroutine */ extern int dgelq2_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *info);

/* Subroutine */ extern int dgelqf_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ extern int dgels_(char *trans, integer *m, integer *n, integer *
	nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb,
	doublereal *work, integer *lwork, integer *info);

/* Subroutine */ extern int dgelsd_(integer *m, integer *n, integer *nrhs,
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	s, doublereal *rcond, integer *rank, doublereal *work, integer *lwork,
	 integer *iwork, integer *info);

/* Subroutine */ extern int dgelss_(integer *m, integer *n, integer *nrhs,
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	s, doublereal *rcond, integer *rank, doublereal *work, integer *lwork,
	 integer *info);

/* Subroutine */ extern int dgelsx_(integer *m, integer *n, integer *nrhs,
	doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	jpvt, doublereal *rcond, integer *rank, doublereal *work, integer *
	info);

/* Subroutine */ extern int dgelsy_(integer *m, integer *n, integer *nrhs,
	doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	jpvt, doublereal *rcond, integer *rank, doublereal *work, integer *
	lwork, integer *info);

/* Subroutine */ extern int dgeql2_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *info);

/* Subroutine */ extern int dgeqlf_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ extern int dgeqp3_(integer *m, integer *n, doublereal *a, integer *
	lda, integer *jpvt, doublereal *tau, doublereal *work, integer *lwork,
	 integer *info);

/* Subroutine */ extern int dgeqpf_(integer *m, integer *n, doublereal *a, integer *
	lda, integer *jpvt, doublereal *tau, doublereal *work, integer *info);

/* Subroutine */ extern int dgeqr2_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *info);

/* Subroutine */ extern int dgeqrf_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ extern int dgerfs_(char *trans, integer *n, integer *nrhs,
	doublereal *a, integer *lda, doublereal *af, integer *ldaf, integer *
	ipiv, doublereal *b, integer *ldb, doublereal *x, integer *ldx,
	doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork,
	integer *info);

/* Subroutine */ extern int dgerq2_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *info);

/* Subroutine */ extern int dgerqf_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ extern int dgesc2_(integer *n, doublereal *a, integer *lda,
	doublereal *rhs, integer *ipiv, integer *jpiv, doublereal *scale);

/* Subroutine */ extern int dgesdd_(char *jobz, integer *m, integer *n, doublereal *
	a, integer *lda, doublereal *s, doublereal *u, integer *ldu,
	doublereal *vt, integer *ldvt, doublereal *work, integer *lwork,
	integer *iwork, integer *info);

/* Subroutine */ extern int dgesv_(integer *n, integer *nrhs, doublereal *a, integer
	*lda, integer *ipiv, doublereal *b, integer *ldb, integer *info);

/* Subroutine */ extern int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n,
	doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *
	ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork,
	integer *info);

/* Subroutine */ extern int dgesvx_(char *fact, char *trans, integer *n, integer *
	nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf,
	integer *ipiv, char *equed, doublereal *r__, doublereal *c__,
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *
	iwork, integer *info);

/* Subroutine */ extern int dgetc2_(integer *n, doublereal *a, integer *lda, integer
	*ipiv, integer *jpiv, integer *info);

/* Subroutine */ extern int dgetf2_(integer *m, integer *n, doublereal *a, integer *
	lda, integer *ipiv, integer *info);

/* Subroutine */ extern int dgetrf_(integer *m, integer *n, doublereal *a, integer *
	lda, integer *ipiv, integer *info);

/* Subroutine */ extern int dgetri_(integer *n, doublereal *a, integer *lda, integer
	*ipiv, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ extern int dgetrs_(char *trans, integer *n, integer *nrhs,
	doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *
	ldb, integer *info);

/* Subroutine */ extern int dggbak_(char *job, char *side, integer *n, integer *ilo,
	integer *ihi, doublereal *lscale, doublereal *rscale, integer *m,
	doublereal *v, integer *ldv, integer *info);

/* Subroutine */ extern int dggbal_(char *job, integer *n, doublereal *a, integer *
	lda, doublereal *b, integer *ldb, integer *ilo, integer *ihi,
	doublereal *lscale, doublereal *rscale, doublereal *work, integer *
	info);

/* Subroutine */ extern int dgges_(char *jobvsl, char *jobvsr, char *sort, L_fp
	delctg, integer *n, doublereal *a, integer *lda, doublereal *b,
	integer *ldb, integer *sdim, doublereal *alphar, doublereal *alphai,
	doublereal *beta, doublereal *vsl, integer *ldvsl, doublereal *vsr,
	integer *ldvsr, doublereal *work, integer *lwork, logical *bwork,
	integer *info);

/* Subroutine */ extern int dggesx_(char *jobvsl, char *jobvsr, char *sort, L_fp
	delctg, char *sense, integer *n, doublereal *a, integer *lda,
	doublereal *b, integer *ldb, integer *sdim, doublereal *alphar,
	doublereal *alphai, doublereal *beta, doublereal *vsl, integer *ldvsl,
	 doublereal *vsr, integer *ldvsr, doublereal *rconde, doublereal *
	rcondv, doublereal *work, integer *lwork, integer *iwork, integer *
	liwork, logical *bwork, integer *info);

/* Subroutine */ extern int dggev_(char *jobvl, char *jobvr, integer *n, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *alphar,
	doublereal *alphai, doublereal *beta, doublereal *vl, integer *ldvl,
	doublereal *vr, integer *ldvr, doublereal *work, integer *lwork,
	integer *info);

/* Subroutine */ extern int dggevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, integer *n, doublereal *a, integer *lda, doublereal *b,
	integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *
	beta, doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr,
	integer *ilo, integer *ihi, doublereal *lscale, doublereal *rscale,
	doublereal *abnrm, doublereal *bbnrm, doublereal *rconde, doublereal *
	rcondv, doublereal *work, integer *lwork, integer *iwork, logical *
	bwork, integer *info);

/* Subroutine */ extern int dggglm_(integer *n, integer *m, integer *p, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *d__,
	doublereal *x, doublereal *y, doublereal *work, integer *lwork,
	integer *info);

/* Subroutine */ extern int dgghrd_(char *compq, char *compz, integer *n, integer *
	ilo, integer *ihi, doublereal *a, integer *lda, doublereal *b,
	integer *ldb, doublereal *q, integer *ldq, doublereal *z__, integer *
	ldz, integer *info);

/* Subroutine */ extern int dgglse_(integer *m, integer *n, integer *p, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *c__,
	doublereal *d__, doublereal *x, doublereal *work, integer *lwork,
	integer *info);

/* Subroutine */ extern int dggqrf_(integer *n, integer *m, integer *p, doublereal *
	a, integer *lda, doublereal *taua, doublereal *b, integer *ldb,
	doublereal *taub, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ extern int dggrqf_(integer *m, integer *p, integer *n, doublereal *
	a, integer *lda, doublereal *taua, doublereal *b, integer *ldb,
	doublereal *taub, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ extern int dggsvd_(char *jobu, char *jobv, char *jobq, integer *m,
	integer *n, integer *p, integer *k, integer *l, doublereal *a,
	integer *lda, doublereal *b, integer *ldb, doublereal *alpha,
	doublereal *beta, doublereal *u, integer *ldu, doublereal *v, integer
	*ldv, doublereal *q, integer *ldq, doublereal *work, integer *iwork,
	integer *info);

/* Subroutine */ extern int dggsvp_(char *jobu, char *jobv, char *jobq, integer *m,
	integer *p, integer *n, doublereal *a, integer *lda, doublereal *b,
	integer *ldb, doublereal *tola, doublereal *tolb, integer *k, integer
	*l, doublereal *u, integer *ldu, doublereal *v, integer *ldv,
	doublereal *q, integer *ldq, integer *iwork, doublereal *tau,
	doublereal *work, integer *info);

/* Subroutine */ extern int dgtcon_(char *norm, integer *n, doublereal *dl,
	doublereal *d__, doublereal *du, doublereal *du2, integer *ipiv,
	doublereal *anorm, doublereal *rcond, doublereal *work, integer *
	iwork, integer *info);

/* Subroutine */ extern int dgtrfs_(char *trans, integer *n, integer *nrhs,
	doublereal *dl, doublereal *d__, doublereal *du, doublereal *dlf,
	doublereal *df, doublereal *duf, doublereal *du2, integer *ipiv,
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	ferr, doublereal *berr, doublereal *work, integer *iwork, integer *
	info);

/* Subroutine */ extern int dgtsv_(integer *n, integer *nrhs, doublereal *dl,
	doublereal *d__, doublereal *du, doublereal *b, integer *ldb, integer
	*info);

/* Subroutine */ extern int dgtsvx_(char *fact, char *trans, integer *n, integer *
	nrhs, doublereal *dl, doublereal *d__, doublereal *du, doublereal *
	dlf, doublereal *df, doublereal *duf, doublereal *du2, integer *ipiv,
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *
	iwork, integer *info);

/* Subroutine */ extern int dgttrf_(integer *n, doublereal *dl, doublereal *d__,
	doublereal *du, doublereal *du2, integer *ipiv, integer *info);

/* Subroutine */ extern int dgttrs_(char *trans, integer *n, integer *nrhs,
	doublereal *dl, doublereal *d__, doublereal *du, doublereal *du2,
	integer *ipiv, doublereal *b, integer *ldb, integer *info);

/* Subroutine */ extern int dgtts2_(integer *itrans, integer *n, integer *nrhs,
	doublereal *dl, doublereal *d__, doublereal *du, doublereal *du2,
	integer *ipiv, doublereal *b, integer *ldb);

/* Subroutine */ extern int dhgeqz_(char *job, char *compq, char *compz, integer *n,
	integer *ilo, integer *ihi, doublereal *a, integer *lda, doublereal *
	b, integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *
	beta, doublereal *q, integer *ldq, doublereal *z__, integer *ldz,
	doublereal *work, integer *lwork, integer *info);

/* Subroutine */ extern int dhsein_(char *side, char *eigsrc, char *initv, logical *
	select, integer *n, doublereal *h__, integer *ldh, doublereal *wr,
	doublereal *wi, doublereal *vl, integer *ldvl, doublereal *vr,
	integer *ldvr, integer *mm, integer *m, doublereal *work, integer *
	ifaill, integer *ifailr, integer *info);

/* Subroutine */ extern int dhseqr_(char *job, char *compz, integer *n, integer *ilo,
	 integer *ihi, doublereal *h__, integer *ldh, doublereal *wr,
	doublereal *wi, doublereal *z__, integer *ldz, doublereal *work,
	integer *lwork, integer *info);

/* Subroutine */ extern int dlabad_(doublereal *small, doublereal *large);

/* Subroutine */ extern int dlabrd_(integer *m, integer *n, integer *nb, doublereal *
	a, integer *lda, doublereal *d__, doublereal *e, doublereal *tauq,
	doublereal *taup, doublereal *x, integer *ldx, doublereal *y, integer
	*ldy);

/* Subroutine */ extern int dlacon_(integer *n, doublereal *v, doublereal *x,
	integer *isgn, doublereal *est, integer *kase);

/* Subroutine */ extern int dlacpy_(char *uplo, integer *m, integer *n, doublereal *
	a, integer *lda, doublereal *b, integer *ldb);

/* Subroutine */ extern int dladiv_(doublereal *a, doublereal *b, doublereal *c__,
	doublereal *d__, doublereal *p, doublereal *q);

/* Subroutine */ extern int dlae2_(doublereal *a, doublereal *b, doublereal *c__,
	doublereal *rt1, doublereal *rt2);

/* Subroutine */ extern int dlaebz_(integer *ijob, integer *nitmax, integer *n,
	integer *mmax, integer *minp, integer *nbmin, doublereal *abstol,
	doublereal *reltol, doublereal *pivmin, doublereal *d__, doublereal *
	e, doublereal *e2, integer *nval, doublereal *ab, doublereal *c__,
	integer *mout, integer *nab, doublereal *work, integer *iwork,
	integer *info);

/* Subroutine */ extern int dlaed0_(integer *icompq, integer *qsiz, integer *n,
	doublereal *d__, doublereal *e, doublereal *q, integer *ldq,
	doublereal *qstore, integer *ldqs, doublereal *work, integer *iwork,
	integer *info);

/* Subroutine */ extern int dlaed1_(integer *n, doublereal *d__, doublereal *q,
	integer *ldq, integer *indxq, doublereal *rho, integer *cutpnt,
	doublereal *work, integer *iwork, integer *info);

/* Subroutine */ extern int dlaed2_(integer *k, integer *n, integer *n1, doublereal *
	d__, doublereal *q, integer *ldq, integer *indxq, doublereal *rho,
	doublereal *z__, doublereal *dlamda, doublereal *w, doublereal *q2,
	integer *indx, integer *indxc, integer *indxp, integer *coltyp,
	integer *info);

/* Subroutine */ extern int dlaed3_(integer *k, integer *n, integer *n1, doublereal *
	d__, doublereal *q, integer *ldq, doublereal *rho, doublereal *dlamda,
	 doublereal *q2, integer *indx, integer *ctot, doublereal *w,
	doublereal *s, integer *info);

/* Subroutine */ extern int dlaed4_(integer *n, integer *i__, doublereal *d__,
	doublereal *z__, doublereal *delta, doublereal *rho, doublereal *dlam,
	 integer *info);

/* Subroutine */ extern int dlaed5_(integer *i__, doublereal *d__, doublereal *z__,
	doublereal *delta, doublereal *rho, doublereal *dlam);

/* Subroutine */ extern int dlaed6_(integer *kniter, logical *orgati, doublereal *
	rho, doublereal *d__, doublereal *z__, doublereal *finit, doublereal *
	tau, integer *info);

/* Subroutine */ extern int dlaed7_(integer *icompq, integer *n, integer *qsiz,
	integer *tlvls, integer *curlvl, integer *curpbm, doublereal *d__,
	doublereal *q, integer *ldq, integer *indxq, doublereal *rho, integer
	*cutpnt, doublereal *qstore, integer *qptr, integer *prmptr, integer *
	perm, integer *givptr, integer *givcol, doublereal *givnum,
	doublereal *work, integer *iwork, integer *info);

/* Subroutine */ extern int dlaed8_(integer *icompq, integer *k, integer *n, integer
	*qsiz, doublereal *d__, doublereal *q, integer *ldq, integer *indxq,
	doublereal *rho, integer *cutpnt, doublereal *z__, doublereal *dlamda,
	 doublereal *q2, integer *ldq2, doublereal *w, integer *perm, integer
	*givptr, integer *givcol, doublereal *givnum, integer *indxp, integer
	*indx, integer *info);

/* Subroutine */ extern int dlaed9_(integer *k, integer *kstart, integer *kstop,
	integer *n, doublereal *d__, doublereal *q, integer *ldq, doublereal *
	rho, doublereal *dlamda, doublereal *w, doublereal *s, integer *lds,
	integer *info);

/* Subroutine */ extern int dlaeda_(integer *n, integer *tlvls, integer *curlvl,
	integer *curpbm, integer *prmptr, integer *perm, integer *givptr,
	integer *givcol, doublereal *givnum, doublereal *q, integer *qptr,
	doublereal *z__, doublereal *ztemp, integer *info);

/* Subroutine */ extern int dlaein_(logical *rightv, logical *noinit, integer *n,
	doublereal *h__, integer *ldh, doublereal *wr, doublereal *wi,
	doublereal *vr, doublereal *vi, doublereal *b, integer *ldb,
	doublereal *work, doublereal *eps3, doublereal *smlnum, doublereal *
	bignum, integer *info);

/* Subroutine */ extern int dlaev2_(doublereal *a, doublereal *b, doublereal *c__,
	doublereal *rt1, doublereal *rt2, doublereal *cs1, doublereal *sn1);

/* Subroutine */ extern int dlaexc_(logical *wantq, integer *n, doublereal *t,
	integer *ldt, doublereal *q, integer *ldq, integer *j1, integer *n1,
	integer *n2, doublereal *work, integer *info);

/* Subroutine */ extern int dlag2_(doublereal *a, integer *lda, doublereal *b,
	integer *ldb, doublereal *safmin, doublereal *scale1, doublereal *
	scale2, doublereal *wr1, doublereal *wr2, doublereal *wi);

/* Subroutine */ extern int dlags2_(logical *upper, doublereal *a1, doublereal *a2,
	doublereal *a3, doublereal *b1, doublereal *b2, doublereal *b3,
	doublereal *csu, doublereal *snu, doublereal *csv, doublereal *snv,
	doublereal *csq, doublereal *snq);

/* Subroutine */ extern int dlagtf_(integer *n, doublereal *a, doublereal *lambda,
	doublereal *b, doublereal *c__, doublereal *tol, doublereal *d__,
	integer *in, integer *info);

/* Subroutine */ extern int dlagtm_(char *trans, integer *n, integer *nrhs,
	doublereal *alpha, doublereal *dl, doublereal *d__, doublereal *du,
	doublereal *x, integer *ldx, doublereal *beta, doublereal *b, integer
	*ldb);

/* Subroutine */ extern int dlagts_(integer *job, integer *n, doublereal *a,
	doublereal *b, doublereal *c__, doublereal *d__, integer *in,
	doublereal *y, doublereal *tol, integer *info);

/* Subroutine */ extern int dlagv2_(doublereal *a, integer *lda, doublereal *b,
	integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *
	beta, doublereal *csl, doublereal *snl, doublereal *csr, doublereal *
	snr);

/* Subroutine */ extern int dlahqr_(logical *wantt, logical *wantz, integer *n,
	integer *ilo, integer *ihi, doublereal *h__, integer *ldh, doublereal
	*wr, doublereal *wi, integer *iloz, integer *ihiz, doublereal *z__,
	integer *ldz, integer *info);

/* Subroutine */ extern int dlahrd_(integer *n, integer *k, integer *nb, doublereal *
	a, integer *lda, doublereal *tau, doublereal *t, integer *ldt,
	doublereal *y, integer *ldy);

/* Subroutine */ extern int dlaic1_(integer *job, integer *j, doublereal *x,
	doublereal *sest, doublereal *w, doublereal *gamma, doublereal *
	sestpr, doublereal *s, doublereal *c__);

/* Subroutine */ extern int dlaln2_(logical *ltrans, integer *na, integer *nw,
	doublereal *smin, doublereal *ca, doublereal *a, integer *lda,
	doublereal *d1, doublereal *d2, doublereal *b, integer *ldb,
	doublereal *wr, doublereal *wi, doublereal *x, integer *ldx,
	doublereal *scale, doublereal *xnorm, integer *info);

/* Subroutine */ extern int dlals0_(integer *icompq, integer *nl, integer *nr,
	integer *sqre, integer *nrhs, doublereal *b, integer *ldb, doublereal
	*bx, integer *ldbx, integer *perm, integer *givptr, integer *givcol,
	integer *ldgcol, doublereal *givnum, integer *ldgnum, doublereal *
	poles, doublereal *difl, doublereal *difr, doublereal *z__, integer *
	k, doublereal *c__, doublereal *s, doublereal *work, integer *info);

/* Subroutine */ extern int dlalsa_(integer *icompq, integer *smlsiz, integer *n,
	integer *nrhs, doublereal *b, integer *ldb, doublereal *bx, integer *
	ldbx, doublereal *u, integer *ldu, doublereal *vt, integer *k,
	doublereal *difl, doublereal *difr, doublereal *z__, doublereal *
	poles, integer *givptr, integer *givcol, integer *ldgcol, integer *
	perm, doublereal *givnum, doublereal *c__, doublereal *s, doublereal *
	work, integer *iwork, integer *info);

/* Subroutine */ extern int dlalsd_(char *uplo, integer *smlsiz, integer *n, integer
	*nrhs, doublereal *d__, doublereal *e, doublereal *b, integer *ldb,
	doublereal *rcond, integer *rank, doublereal *work, integer *iwork,
	integer *info);

/* Subroutine */ extern int dlamc1_(integer *beta, integer *t, logical *rnd, logical
	*ieee1);

/* Subroutine */ extern int dlamc2_(integer *beta, integer *t, logical *rnd,
	doublereal *eps, integer *emin, doublereal *rmin, integer *emax,
	doublereal *rmax);

/* Subroutine */ extern int dlamc4_(integer *emin, doublereal *start, integer *base);

/* Subroutine */ extern int dlamc5_(integer *beta, integer *p, integer *emin,
	logical *ieee, integer *emax, doublereal *rmax);

/* Subroutine */ extern int dlamrg_(integer *n1, integer *n2, doublereal *a, integer
	*dtrd1, integer *dtrd2, integer *index);

/* Subroutine */ extern int dlanv2_(doublereal *a, doublereal *b, doublereal *c__,
	doublereal *d__, doublereal *rt1r, doublereal *rt1i, doublereal *rt2r,
	 doublereal *rt2i, doublereal *cs, doublereal *sn);

/* Subroutine */ extern int dlapll_(integer *n, doublereal *x, integer *incx,
	doublereal *y, integer *incy, doublereal *ssmin);

/* Subroutine */ extern int dlapmt_(logical *forwrd, integer *m, integer *n,
	doublereal *x, integer *ldx, integer *k);

/* Subroutine */ extern int dlaqgb_(integer *m, integer *n, integer *kl, integer *ku,
	 doublereal *ab, integer *ldab, doublereal *r__, doublereal *c__,
	doublereal *rowcnd, doublereal *colcnd, doublereal *amax, char *equed);

/* Subroutine */ extern int dlaqge_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *r__, doublereal *c__, doublereal *rowcnd, doublereal
	*colcnd, doublereal *amax, char *equed);

/* Subroutine */ extern int dlaqp2_(integer *m, integer *n, integer *offset,
	doublereal *a, integer *lda, integer *jpvt, doublereal *tau,
	doublereal *vn1, doublereal *vn2, doublereal *work);

/* Subroutine */ extern int dlaqps_(integer *m, integer *n, integer *offset, integer
	*nb, integer *kb, doublereal *a, integer *lda, integer *jpvt,
	doublereal *tau, doublereal *vn1, doublereal *vn2, doublereal *auxv,
	doublereal *f, integer *ldf);

/* Subroutine */ extern int dlaqsb_(char *uplo, integer *n, integer *kd, doublereal *
	ab, integer *ldab, doublereal *s, doublereal *scond, doublereal *amax,
	 char *equed);

/* Subroutine */ extern int dlaqsp_(char *uplo, integer *n, doublereal *ap,
	doublereal *s, doublereal *scond, doublereal *amax, char *equed);

/* Subroutine */ extern int dlaqsy_(char *uplo, integer *n, doublereal *a, integer *
	lda, doublereal *s, doublereal *scond, doublereal *amax, char *equed);

/* Subroutine */ extern int dlaqtr_(logical *ltran, logical *lreal, integer *n,
	doublereal *t, integer *ldt, doublereal *b, doublereal *w, doublereal
	*scale, doublereal *x, doublereal *work, integer *info);

/* Subroutine */ extern int dlar1v_(integer *n, integer *b1, integer *bn, doublereal
	*sigma, doublereal *d__, doublereal *l, doublereal *ld, doublereal *
	lld, doublereal *gersch, doublereal *z__, doublereal *ztz, doublereal
	*mingma, integer *r__, integer *isuppz, doublereal *work);

/* Subroutine */ extern int dlar2v_(integer *n, doublereal *x, doublereal *y,
	doublereal *z__, integer *incx, doublereal *c__, doublereal *s,
	integer *incc);

/* Subroutine */ extern int dlarf_(char *side, integer *m, integer *n, doublereal *v,
	 integer *incv, doublereal *tau, doublereal *c__, integer *ldc,
	doublereal *work);

/* Subroutine */ extern int dlarfb_(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, doublereal *v, integer *
	ldv, doublereal *t, integer *ldt, doublereal *c__, integer *ldc,
	doublereal *work, integer *ldwork);

/* Subroutine */ extern int dlarfg_(integer *n, doublereal *alpha, doublereal *x,
	integer *incx, doublereal *tau);

/* Subroutine */ extern int dlarft_(char *direct, char *storev, integer *n, integer *
	k, doublereal *v, integer *ldv, doublereal *tau, doublereal *t,
	integer *ldt);

/* Subroutine */ extern int dlarfx_(char *side, integer *m, integer *n, doublereal *
	v, doublereal *tau, doublereal *c__, integer *ldc, doublereal *work);

/* Subroutine */ extern int dlargv_(integer *n, doublereal *x, integer *incx,
	doublereal *y, integer *incy, doublereal *c__, integer *incc);

/* Subroutine */ extern int dlarnv_(integer *idist, integer *iseed, integer *n,
	doublereal *x);

/* Subroutine */ extern int dlarrb_(integer *n, doublereal *d__, doublereal *l,
	doublereal *ld, doublereal *lld, integer *ifirst, integer *ilast,
	doublereal *sigma, doublereal *reltol, doublereal *w, doublereal *
	wgap, doublereal *werr, doublereal *work, integer *iwork, integer *
	info);

/* Subroutine */ extern int dlarre_(integer *n, doublereal *d__, doublereal *e,
	doublereal *tol, integer *nsplit, integer *isplit, integer *m,
	doublereal *w, doublereal *woff, doublereal *gersch, doublereal *work,
	 integer *info);

/* Subroutine */ extern int dlarrf_(integer *n, doublereal *d__, doublereal *l,
	doublereal *ld, doublereal *lld, integer *ifirst, integer *ilast,
	doublereal *w, doublereal *dplus, doublereal *lplus, doublereal *work,
	 integer *iwork, integer *info);

/* Subroutine */ extern int dlarrv_(integer *n, doublereal *d__, doublereal *l,
	integer *isplit, integer *m, doublereal *w, integer *iblock,
	doublereal *gersch, doublereal *tol, doublereal *z__, integer *ldz,
	integer *isuppz, doublereal *work, integer *iwork, integer *info);

/* Subroutine */ extern int dlartg_(doublereal *f, doublereal *g, doublereal *cs,
	doublereal *sn, doublereal *r__);

/* Subroutine */ extern int dlartv_(integer *n, doublereal *x, integer *incx,
	doublereal *y, integer *incy, doublereal *c__, doublereal *s, integer
	*incc);

/* Subroutine */ extern int dlaruv_(integer *iseed, integer *n, doublereal *x);

/* Subroutine */ extern int dlarz_(char *side, integer *m, integer *n, integer *l,
	doublereal *v, integer *incv, doublereal *tau, doublereal *c__,
	integer *ldc, doublereal *work);

/* Subroutine */ extern int dlarzb_(char *side, char *trans, char *direct, char *
	storev, integer *m, integer *n, integer *k, integer *l, doublereal *v,
	 integer *ldv, doublereal *t, integer *ldt, doublereal *c__, integer *
	ldc, doublereal *work, integer *ldwork);

/* Subroutine */ extern int dlarzt_(char *direct, char *storev, integer *n, integer *
	k, doublereal *v, integer *ldv, doublereal *tau, doublereal *t,
	integer *ldt);

/* Subroutine */ extern int dlas2_(doublereal *f, doublereal *g, doublereal *h__,
	doublereal *ssmin, doublereal *ssmax);

/* Subroutine */ extern int dlascl_(char *type__, integer *kl, integer *ku,
	doublereal *cfrom, doublereal *cto, integer *m, integer *n,
	doublereal *a, integer *lda, integer *info);

/* Subroutine */ extern int dlasd0_(integer *n, integer *sqre, doublereal *d__,
	doublereal *e, doublereal *u, integer *ldu, doublereal *vt, integer *
	ldvt, integer *smlsiz, integer *iwork, doublereal *work, integer *
	info);

/* Subroutine */ extern int dlasd1_(integer *nl, integer *nr, integer *sqre,
	doublereal *d__, doublereal *alpha, doublereal *beta, doublereal *u,
	integer *ldu, doublereal *vt, integer *ldvt, integer *idxq, integer *
	iwork, doublereal *work, integer *info);

/* Subroutine */ extern int dlasd2_(integer *nl, integer *nr, integer *sqre, integer
	*k, doublereal *d__, doublereal *z__, doublereal *alpha, doublereal *
	beta, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt,
	doublereal *dsigma, doublereal *u2, integer *ldu2, doublereal *vt2,
	integer *ldvt2, integer *idxp, integer *idx, integer *idxc, integer *
	idxq, integer *coltyp, integer *info);

/* Subroutine */ extern int dlasd3_(integer *nl, integer *nr, integer *sqre, integer
	*k, doublereal *d__, doublereal *q, integer *ldq, doublereal *dsigma,
	doublereal *u, integer *ldu, doublereal *u2, integer *ldu2,
	doublereal *vt, integer *ldvt, doublereal *vt2, integer *ldvt2,
	integer *idxc, integer *ctot, doublereal *z__, integer *info);

/* Subroutine */ extern int dlasd4_(integer *n, integer *i__, doublereal *d__,
	doublereal *z__, doublereal *delta, doublereal *rho, doublereal *
	sigma, doublereal *work, integer *info);

/* Subroutine */ extern int dlasd5_(integer *i__, doublereal *d__, doublereal *z__,
	doublereal *delta, doublereal *rho, doublereal *dsigma, doublereal *
	work);

/* Subroutine */ extern int dlasd6_(integer *icompq, integer *nl, integer *nr,
	integer *sqre, doublereal *d__, doublereal *vf, doublereal *vl,
	doublereal *alpha, doublereal *beta, integer *idxq, integer *perm,
	integer *givptr, integer *givcol, integer *ldgcol, doublereal *givnum,
	 integer *ldgnum, doublereal *poles, doublereal *difl, doublereal *
	difr, doublereal *z__, integer *k, doublereal *c__, doublereal *s,
	doublereal *work, integer *iwork, integer *info);

/* Subroutine */ extern int dlasd7_(integer *icompq, integer *nl, integer *nr,
	integer *sqre, integer *k, doublereal *d__, doublereal *z__,
	doublereal *zw, doublereal *vf, doublereal *vfw, doublereal *vl,
	doublereal *vlw, doublereal *alpha, doublereal *beta, doublereal *
	dsigma, integer *idx, integer *idxp, integer *idxq, integer *perm,
	integer *givptr, integer *givcol, integer *ldgcol, doublereal *givnum,
	 integer *ldgnum, doublereal *c__, doublereal *s, integer *info);

/* Subroutine */ extern int dlasd8_(integer *icompq, integer *k, doublereal *d__,
	doublereal *z__, doublereal *vf, doublereal *vl, doublereal *difl,
	doublereal *difr, integer *lddifr, doublereal *dsigma, doublereal *
	work, integer *info);

/* Subroutine */ extern int dlasd9_(integer *icompq, integer *ldu, integer *k,
	doublereal *d__, doublereal *z__, doublereal *vf, doublereal *vl,
	doublereal *difl, doublereal *difr, doublereal *dsigma, doublereal *
	work, integer *info);

/* Subroutine */ extern int dlasda_(integer *icompq, integer *smlsiz, integer *n,
	integer *sqre, doublereal *d__, doublereal *e, doublereal *u, integer
	*ldu, doublereal *vt, integer *k, doublereal *difl, doublereal *difr,
	doublereal *z__, doublereal *poles, integer *givptr, integer *givcol,
	integer *ldgcol, integer *perm, doublereal *givnum, doublereal *c__,
	doublereal *s, doublereal *work, integer *iwork, integer *info);

/* Subroutine */ extern int dlasdq_(char *uplo, integer *sqre, integer *n, integer *
	ncvt, integer *nru, integer *ncc, doublereal *d__, doublereal *e,
	doublereal *vt, integer *ldvt, doublereal *u, integer *ldu,
	doublereal *c__, integer *ldc, doublereal *work, integer *info);

/* Subroutine */ extern int dlasdt_(integer *n, integer *lvl, integer *nd, integer *
	inode, integer *ndiml, integer *ndimr, integer *msub);

/* Subroutine */ extern int dlaset_(char *uplo, integer *m, integer *n, doublereal *
	alpha, doublereal *beta, doublereal *a, integer *lda);

/* Subroutine */ extern int dlasq1_(integer *n, doublereal *d__, doublereal *e,
	doublereal *work, integer *info);

/* Subroutine */ extern int dlasq2_(integer *n, doublereal *z__, integer *info);

/* Subroutine */ extern int dlasq3_(integer *i0, integer *n0, doublereal *z__,
	integer *pp, doublereal *dmin__, doublereal *sigma, doublereal *desig,
	 doublereal *qmax, integer *nfail, integer *iter, integer *ndiv,
	logical *ieee);

/* Subroutine */ extern int dlasq4_(integer *i0, integer *n0, doublereal *z__,
	integer *pp, integer *n0in, doublereal *dmin__, doublereal *dmin1,
	doublereal *dmin2, doublereal *dn, doublereal *dn1, doublereal *dn2,
	doublereal *tau, integer *ttype);

/* Subroutine */ extern int dlasq5_(integer *i0, integer *n0, doublereal *z__,
	integer *pp, doublereal *tau, doublereal *dmin__, doublereal *dmin1,
	doublereal *dmin2, doublereal *dn, doublereal *dnm1, doublereal *dnm2,
	 logical *ieee);

/* Subroutine */ extern int dlasq6_(integer *i0, integer *n0, doublereal *z__,
	integer *pp, doublereal *dmin__, doublereal *dmin1, doublereal *dmin2,
	 doublereal *dn, doublereal *dnm1, doublereal *dnm2);

/* Subroutine */ extern int dlasr_(char *side, char *pivot, char *direct, integer *m,
	 integer *n, doublereal *c__, doublereal *s, doublereal *a, integer *
	lda);

/* Subroutine */ extern int dlasrt_(char *id, integer *n, doublereal *d__, integer *
	info);

/* Subroutine */ extern int dlassq_(integer *n, doublereal *x, integer *incx,
	doublereal *scale, doublereal *sumsq);

/* Subroutine */ extern int dlasv2_(doublereal *f, doublereal *g, doublereal *h__,
	doublereal *ssmin, doublereal *ssmax, doublereal *snr, doublereal *
	csr, doublereal *snl, doublereal *csl);

/* Subroutine */ extern int dlaswp_(integer *n, doublereal *a, integer *lda, integer
	*k1, integer *k2, integer *ipiv, integer *incx);

/* Subroutine */ extern int dlasy2_(logical *ltranl, logical *ltranr, integer *isgn,
	integer *n1, integer *n2, doublereal *tl, integer *ldtl, doublereal *
	tr, integer *ldtr, doublereal *b, integer *ldb, doublereal *scale,
	doublereal *x, integer *ldx, doublereal *xnorm, integer *info);

/* Subroutine */ extern int dlasyf_(char *uplo, integer *n, integer *nb, integer *kb,
	 doublereal *a, integer *lda, integer *ipiv, doublereal *w, integer *
	ldw, integer *info);

/* Subroutine */ extern int dlatbs_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, integer *kd, doublereal *ab, integer *ldab,
	doublereal *x, doublereal *scale, doublereal *cnorm, integer *info);

/* Subroutine */ extern int dlatdf_(integer *ijob, integer *n, doublereal *z__,
	integer *ldz, doublereal *rhs, doublereal *rdsum, doublereal *rdscal,
	integer *ipiv, integer *jpiv);

/* Subroutine */ extern int dlatps_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, doublereal *ap, doublereal *x, doublereal *scale,
	doublereal *cnorm, integer *info);

/* Subroutine */ extern int dlatrd_(char *uplo, integer *n, integer *nb, doublereal *
	a, integer *lda, doublereal *e, doublereal *tau, doublereal *w,
	integer *ldw);

/* Subroutine */ extern int dlatrs_(char *uplo, char *trans, char *diag, char *
	normin, integer *n, doublereal *a, integer *lda, doublereal *x,
	doublereal *scale, doublereal *cnorm, integer *info);

/* Subroutine */ extern int dlatrz_(integer *m, integer *n, integer *l, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work);

/* Subroutine */ extern int dlatzm_(char *side, integer *m, integer *n, doublereal *
	v, integer *incv, doublereal *tau, doublereal *c1, doublereal *c2,
	integer *ldc, doublereal *work);

/* Subroutine */ extern int dlauu2_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *info);

/* Subroutine */ extern int dlauum_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *info);

/* Subroutine */ extern int dopgtr_(char *uplo, integer *n, doublereal *ap,
	doublereal *tau, doublereal *q, integer *ldq, doublereal *work,
	integer *info);

/* Subroutine */ extern int dopmtr_(char *side, char *uplo, char *trans, integer *m,
	integer *n, doublereal *ap, doublereal *tau, doublereal *c__, integer
	*ldc, doublereal *work, integer *info);

/* Subroutine */ extern int dorg2l_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *info);

/* Subroutine */ extern int dorg2r_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *info);

/* Subroutine */ extern int dorgbr_(char *vect, integer *m, integer *n, integer *k,
	doublereal *a, integer *lda, doublereal *tau, doublereal *work,
	integer *lwork, integer *info);

/* Subroutine */ extern int dorghr_(integer *n, integer *ilo, integer *ihi,
	doublereal *a, integer *lda, doublereal *tau, doublereal *work,
	integer *lwork, integer *info);

/* Subroutine */ extern int dorgl2_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *info);

/* Subroutine */ extern int dorglq_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *lwork,
	integer *info);

/* Subroutine */ extern int dorgql_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *lwork,
	integer *info);

/* Subroutine */ extern int dorgqr_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *lwork,
	integer *info);

/* Subroutine */ extern int dorgr2_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *info);

/* Subroutine */ extern int dorgrq_(integer *m, integer *n, integer *k, doublereal *
	a, integer *lda, doublereal *tau, doublereal *work, integer *lwork,
	integer *info);

/* Subroutine */ extern int dorgtr_(char *uplo, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ extern int dorm2l_(char *side, char *trans, integer *m, integer *n,
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *info);

/* Subroutine */ extern int dorm2r_(char *side, char *trans, integer *m, integer *n,
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *info);

/* Subroutine */ extern int dormbr_(char *vect, char *side, char *trans, integer *m,
	integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau,
	doublereal *c__, integer *ldc, doublereal *work, integer *lwork,
	integer *info);

/* Subroutine */ extern int dormhr_(char *side, char *trans, integer *m, integer *n,
	integer *ilo, integer *ihi, doublereal *a, integer *lda, doublereal *
	tau, doublereal *c__, integer *ldc, doublereal *work, integer *lwork,
	integer *info);

/* Subroutine */ extern int dorml2_(char *side, char *trans, integer *m, integer *n,
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *info);

/* Subroutine */ extern int dormlq_(char *side, char *trans, integer *m, integer *n,
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ extern int dormql_(char *side, char *trans, integer *m, integer *n,
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ extern int dormqr_(char *side, char *trans, integer *m, integer *n,
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ extern int dormr2_(char *side, char *trans, integer *m, integer *n,
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *info);

/* Subroutine */ extern int dormr3_(char *side, char *trans, integer *m, integer *n,
	integer *k, integer *l, doublereal *a, integer *lda, doublereal *tau,
	doublereal *c__, integer *ldc, doublereal *work, integer *info);

/* Subroutine */ extern int dormrq_(char *side, char *trans, integer *m, integer *n,
	integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ extern int dormrz_(char *side, char *trans, integer *m, integer *n,
	integer *k, integer *l, doublereal *a, integer *lda, doublereal *tau,
	doublereal *c__, integer *ldc, doublereal *work, integer *lwork,
	integer *info);

/* Subroutine */ extern int dormtr_(char *side, char *uplo, char *trans, integer *m,
	integer *n, doublereal *a, integer *lda, doublereal *tau, doublereal *
	c__, integer *ldc, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ extern int dpbcon_(char *uplo, integer *n, integer *kd, doublereal *
	ab, integer *ldab, doublereal *anorm, doublereal *rcond, doublereal *
	work, integer *iwork, integer *info);

/* Subroutine */ extern int dpbequ_(char *uplo, integer *n, integer *kd, doublereal *
	ab, integer *ldab, doublereal *s, doublereal *scond, doublereal *amax,
	 integer *info);

/* Subroutine */ extern int dpbrfs_(char *uplo, integer *n, integer *kd, integer *
	nrhs, doublereal *ab, integer *ldab, doublereal *afb, integer *ldafb,
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	ferr, doublereal *berr, doublereal *work, integer *iwork, integer *
	info);

/* Subroutine */ extern int dpbstf_(char *uplo, integer *n, integer *kd, doublereal *
	ab, integer *ldab, integer *info);

/* Subroutine */ extern int dpbsv_(char *uplo, integer *n, integer *kd, integer *
	nrhs, doublereal *ab, integer *ldab, doublereal *b, integer *ldb,
	integer *info);

/* Subroutine */ extern int dpbsvx_(char *fact, char *uplo, integer *n, integer *kd,
	integer *nrhs, doublereal *ab, integer *ldab, doublereal *afb,
	integer *ldafb, char *equed, doublereal *s, doublereal *b, integer *
	ldb, doublereal *x, integer *ldx, doublereal *rcond, doublereal *ferr,
	 doublereal *berr, doublereal *work, integer *iwork, integer *info);

/* Subroutine */ extern int dpbtf2_(char *uplo, integer *n, integer *kd, doublereal *
	ab, integer *ldab, integer *info);

/* Subroutine */ extern int dpbtrf_(char *uplo, integer *n, integer *kd, doublereal *
	ab, integer *ldab, integer *info);

/* Subroutine */ extern int dpbtrs_(char *uplo, integer *n, integer *kd, integer *
	nrhs, doublereal *ab, integer *ldab, doublereal *b, integer *ldb,
	integer *info);

/* Subroutine */ extern int dpocon_(char *uplo, integer *n, doublereal *a, integer *
	lda, doublereal *anorm, doublereal *rcond, doublereal *work, integer *
	iwork, integer *info);

/* Subroutine */ extern int dpoequ_(integer *n, doublereal *a, integer *lda,
	doublereal *s, doublereal *scond, doublereal *amax, integer *info);

/* Subroutine */ extern int dporfs_(char *uplo, integer *n, integer *nrhs,
	doublereal *a, integer *lda, doublereal *af, integer *ldaf,
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	ferr, doublereal *berr, doublereal *work, integer *iwork, integer *
	info);

/* Subroutine */ extern int dposv_(char *uplo, integer *n, integer *nrhs, doublereal
	*a, integer *lda, doublereal *b, integer *ldb, integer *info);

/* Subroutine */ extern int dposvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf,
	char *equed, doublereal *s, doublereal *b, integer *ldb, doublereal *
	x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *
	berr, doublereal *work, integer *iwork, integer *info);

/* Subroutine */ extern int dpotf2_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *info);

/* Subroutine */ extern int dpotrf_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *info);

/* Subroutine */ extern int dpotri_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *info);

/* Subroutine */ extern int dpotrs_(char *uplo, integer *n, integer *nrhs,
	doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	info);

/* Subroutine */ extern int dppcon_(char *uplo, integer *n, doublereal *ap,
	doublereal *anorm, doublereal *rcond, doublereal *work, integer *
	iwork, integer *info);

/* Subroutine */ extern int dppequ_(char *uplo, integer *n, doublereal *ap,
	doublereal *s, doublereal *scond, doublereal *amax, integer *info);

/* Subroutine */ extern int dpprfs_(char *uplo, integer *n, integer *nrhs,
	doublereal *ap, doublereal *afp, doublereal *b, integer *ldb,
	doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr,
	doublereal *work, integer *iwork, integer *info);

/* Subroutine */ extern int dppsv_(char *uplo, integer *n, integer *nrhs, doublereal
	*ap, doublereal *b, integer *ldb, integer *info);

/* Subroutine */ extern int dppsvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublereal *ap, doublereal *afp, char *equed, doublereal *s,
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *
	iwork, integer *info);

/* Subroutine */ extern int dpptrf_(char *uplo, integer *n, doublereal *ap, integer *
	info);

/* Subroutine */ extern int dpptri_(char *uplo, integer *n, doublereal *ap, integer *
	info);

/* Subroutine */ extern int dpptrs_(char *uplo, integer *n, integer *nrhs,
	doublereal *ap, doublereal *b, integer *ldb, integer *info);

/* Subroutine */ extern int dptcon_(integer *n, doublereal *d__, doublereal *e,
	doublereal *anorm, doublereal *rcond, doublereal *work, integer *info);

/* Subroutine */ extern int dpteqr_(char *compz, integer *n, doublereal *d__,
	doublereal *e, doublereal *z__, integer *ldz, doublereal *work,
	integer *info);

/* Subroutine */ extern int dptrfs_(integer *n, integer *nrhs, doublereal *d__,
	doublereal *e, doublereal *df, doublereal *ef, doublereal *b, integer
	*ldb, doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr,
	 doublereal *work, integer *info);

/* Subroutine */ extern int dptsv_(integer *n, integer *nrhs, doublereal *d__,
	doublereal *e, doublereal *b, integer *ldb, integer *info);

/* Subroutine */ extern int dptsvx_(char *fact, integer *n, integer *nrhs,
	doublereal *d__, doublereal *e, doublereal *df, doublereal *ef,
	doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *
	rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *
	info);

/* Subroutine */ extern int dpttrf_(integer *n, doublereal *d__, doublereal *e,
	integer *info);

/* Subroutine */ extern int dpttrs_(integer *n, integer *nrhs, doublereal *d__,
	doublereal *e, doublereal *b, integer *ldb, integer *info);

/* Subroutine */ extern int dptts2_(integer *n, integer *nrhs, doublereal *d__,
	doublereal *e, doublereal *b, integer *ldb);

/* Subroutine */ extern int drscl_(integer *n, doublereal *sa, doublereal *sx,
	integer *incx);

/* Subroutine */ extern int dsbev_(char *jobz, char *uplo, integer *n, integer *kd,
	doublereal *ab, integer *ldab, doublereal *w, doublereal *z__,
	integer *ldz, doublereal *work, integer *info);

/* Subroutine */ extern int dsbevd_(char *jobz, char *uplo, integer *n, integer *kd,
	doublereal *ab, integer *ldab, doublereal *w, doublereal *z__,
	integer *ldz, doublereal *work, integer *lwork, integer *iwork,
	integer *liwork, integer *info);

/* Subroutine */ extern int dsbevx_(char *jobz, char *range, char *uplo, integer *n,
	integer *kd, doublereal *ab, integer *ldab, doublereal *q, integer *
	ldq, doublereal *vl, doublereal *vu, integer *il, integer *iu,
	doublereal *abstol, integer *m, doublereal *w, doublereal *z__,
	integer *ldz, doublereal *work, integer *iwork, integer *ifail,
	integer *info);

/* Subroutine */ extern int dsbgst_(char *vect, char *uplo, integer *n, integer *ka,
	integer *kb, doublereal *ab, integer *ldab, doublereal *bb, integer *
	ldbb, doublereal *x, integer *ldx, doublereal *work, integer *info);

/* Subroutine */ extern int dsbgv_(char *jobz, char *uplo, integer *n, integer *ka,
	integer *kb, doublereal *ab, integer *ldab, doublereal *bb, integer *
	ldbb, doublereal *w, doublereal *z__, integer *ldz, doublereal *work,
	integer *info);

/* Subroutine */ extern int dsbgvd_(char *jobz, char *uplo, integer *n, integer *ka,
	integer *kb, doublereal *ab, integer *ldab, doublereal *bb, integer *
	ldbb, doublereal *w, doublereal *z__, integer *ldz, doublereal *work,
	integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ extern int dsbgvx_(char *jobz, char *range, char *uplo, integer *n,
	integer *ka, integer *kb, doublereal *ab, integer *ldab, doublereal *
	bb, integer *ldbb, doublereal *q, integer *ldq, doublereal *vl,
	doublereal *vu, integer *il, integer *iu, doublereal *abstol, integer
	*m, doublereal *w, doublereal *z__, integer *ldz, doublereal *work,
	integer *iwork, integer *ifail, integer *info);

/* Subroutine */ extern int dsbtrd_(char *vect, char *uplo, integer *n, integer *kd,
	doublereal *ab, integer *ldab, doublereal *d__, doublereal *e,
	doublereal *q, integer *ldq, doublereal *work, integer *info);

/* Subroutine */ extern int dspcon_(char *uplo, integer *n, doublereal *ap, integer *
	ipiv, doublereal *anorm, doublereal *rcond, doublereal *work, integer
	*iwork, integer *info);

/* Subroutine */ extern int dspev_(char *jobz, char *uplo, integer *n, doublereal *
	ap, doublereal *w, doublereal *z__, integer *ldz, doublereal *work,
	integer *info);

/* Subroutine */ extern int dspevd_(char *jobz, char *uplo, integer *n, doublereal *
	ap, doublereal *w, doublereal *z__, integer *ldz, doublereal *work,
	integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ extern int dspevx_(char *jobz, char *range, char *uplo, integer *n,
	doublereal *ap, doublereal *vl, doublereal *vu, integer *il, integer *
	iu, doublereal *abstol, integer *m, doublereal *w, doublereal *z__,
	integer *ldz, doublereal *work, integer *iwork, integer *ifail,
	integer *info);

/* Subroutine */ extern int dspgst_(integer *itype, char *uplo, integer *n,
	doublereal *ap, doublereal *bp, integer *info);

/* Subroutine */ extern int dspgv_(integer *itype, char *jobz, char *uplo, integer *
	n, doublereal *ap, doublereal *bp, doublereal *w, doublereal *z__,
	integer *ldz, doublereal *work, integer *info);

/* Subroutine */ extern int dspgvd_(integer *itype, char *jobz, char *uplo, integer *
	n, doublereal *ap, doublereal *bp, doublereal *w, doublereal *z__,
	integer *ldz, doublereal *work, integer *lwork, integer *iwork,
	integer *liwork, integer *info);

/* Subroutine */ extern int dspgvx_(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, doublereal *ap, doublereal *bp, doublereal *vl,
	doublereal *vu, integer *il, integer *iu, doublereal *abstol, integer
	*m, doublereal *w, doublereal *z__, integer *ldz, doublereal *work,
	integer *iwork, integer *ifail, integer *info);

/* Subroutine */ extern int dsprfs_(char *uplo, integer *n, integer *nrhs,
	doublereal *ap, doublereal *afp, integer *ipiv, doublereal *b,
	integer *ldb, doublereal *x, integer *ldx, doublereal *ferr,
	doublereal *berr, doublereal *work, integer *iwork, integer *info);

/* Subroutine */ extern int dspsv_(char *uplo, integer *n, integer *nrhs, doublereal
	*ap, integer *ipiv, doublereal *b, integer *ldb, integer *info);

/* Subroutine */ extern int dspsvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublereal *ap, doublereal *afp, integer *ipiv, doublereal *b,
	integer *ldb, doublereal *x, integer *ldx, doublereal *rcond,
	doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork,
	integer *info);

/* Subroutine */ extern int dsptrd_(char *uplo, integer *n, doublereal *ap,
	doublereal *d__, doublereal *e, doublereal *tau, integer *info);

/* Subroutine */ extern int dsptrf_(char *uplo, integer *n, doublereal *ap, integer *
	ipiv, integer *info);

/* Subroutine */ extern int dsptri_(char *uplo, integer *n, doublereal *ap, integer *
	ipiv, doublereal *work, integer *info);

/* Subroutine */ extern int dsptrs_(char *uplo, integer *n, integer *nrhs,
	doublereal *ap, integer *ipiv, doublereal *b, integer *ldb, integer *
	info);

/* Subroutine */ extern int dstebz_(char *range, char *order, integer *n, doublereal
	*vl, doublereal *vu, integer *il, integer *iu, doublereal *abstol,
	doublereal *d__, doublereal *e, integer *m, integer *nsplit,
	doublereal *w, integer *iblock, integer *isplit, doublereal *work,
	integer *iwork, integer *info);

/* Subroutine */ extern int dstedc_(char *compz, integer *n, doublereal *d__,
	doublereal *e, doublereal *z__, integer *ldz, doublereal *work,
	integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ extern int dstegr_(char *jobz, char *range, integer *n, doublereal *
	d__, doublereal *e, doublereal *vl, doublereal *vu, integer *il,
	integer *iu, doublereal *abstol, integer *m, doublereal *w,
	doublereal *z__, integer *ldz, integer *isuppz, doublereal *work,
	integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ extern int dstein_(integer *n, doublereal *d__, doublereal *e,
	integer *m, doublereal *w, integer *iblock, integer *isplit,
	doublereal *z__, integer *ldz, doublereal *work, integer *iwork,
	integer *ifail, integer *info);

/* Subroutine */ extern int dsteqr_(char *compz, integer *n, doublereal *d__,
	doublereal *e, doublereal *z__, integer *ldz, doublereal *work,
	integer *info);

/* Subroutine */ extern int dsterf_(integer *n, doublereal *d__, doublereal *e,
	integer *info);

/* Subroutine */ extern int dstev_(char *jobz, integer *n, doublereal *d__,
	doublereal *e, doublereal *z__, integer *ldz, doublereal *work,
	integer *info);

/* Subroutine */ extern int dstevd_(char *jobz, integer *n, doublereal *d__,
	doublereal *e, doublereal *z__, integer *ldz, doublereal *work,
	integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ extern int dstevr_(char *jobz, char *range, integer *n, doublereal *
	d__, doublereal *e, doublereal *vl, doublereal *vu, integer *il,
	integer *iu, doublereal *abstol, integer *m, doublereal *w,
	doublereal *z__, integer *ldz, integer *isuppz, doublereal *work,
	integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ extern int dstevx_(char *jobz, char *range, integer *n, doublereal *
	d__, doublereal *e, doublereal *vl, doublereal *vu, integer *il,
	integer *iu, doublereal *abstol, integer *m, doublereal *w,
	doublereal *z__, integer *ldz, doublereal *work, integer *iwork,
	integer *ifail, integer *info);

/* Subroutine */ extern int dsycon_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *ipiv, doublereal *anorm, doublereal *rcond, doublereal *
	work, integer *iwork, integer *info);

/* Subroutine */ extern int dsyev_(char *jobz, char *uplo, integer *n, doublereal *a,
	 integer *lda, doublereal *w, doublereal *work, integer *lwork,
	integer *info);

/* Subroutine */ extern int dsyevd_(char *jobz, char *uplo, integer *n, doublereal *
	a, integer *lda, doublereal *w, doublereal *work, integer *lwork,
	integer *iwork, integer *liwork, integer *info);

/* Subroutine */ extern int dsyevr_(char *jobz, char *range, char *uplo, integer *n,
	doublereal *a, integer *lda, doublereal *vl, doublereal *vu, integer *
	il, integer *iu, doublereal *abstol, integer *m, doublereal *w,
	doublereal *z__, integer *ldz, integer *isuppz, doublereal *work,
	integer *lwork, integer *iwork, integer *liwork, integer *info);

/* Subroutine */ extern int dsyevx_(char *jobz, char *range, char *uplo, integer *n,
	doublereal *a, integer *lda, doublereal *vl, doublereal *vu, integer *
	il, integer *iu, doublereal *abstol, integer *m, doublereal *w,
	doublereal *z__, integer *ldz, doublereal *work, integer *lwork,
	integer *iwork, integer *ifail, integer *info);

/* Subroutine */ extern int dsygs2_(integer *itype, char *uplo, integer *n,
	doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	info);

/* Subroutine */ extern int dsygst_(integer *itype, char *uplo, integer *n,
	doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *
	info);

/* Subroutine */ extern int dsygv_(integer *itype, char *jobz, char *uplo, integer *
	n, doublereal *a, integer *lda, doublereal *b, integer *ldb,
	doublereal *w, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ extern int dsygvd_(integer *itype, char *jobz, char *uplo, integer *
	n, doublereal *a, integer *lda, doublereal *b, integer *ldb,
	doublereal *w, doublereal *work, integer *lwork, integer *iwork,
	integer *liwork, integer *info);

/* Subroutine */ extern int dsygvx_(integer *itype, char *jobz, char *range, char *
	uplo, integer *n, doublereal *a, integer *lda, doublereal *b, integer
	*ldb, doublereal *vl, doublereal *vu, integer *il, integer *iu,
	doublereal *abstol, integer *m, doublereal *w, doublereal *z__,
	integer *ldz, doublereal *work, integer *lwork, integer *iwork,
	integer *ifail, integer *info);

/* Subroutine */ extern int dsyrfs_(char *uplo, integer *n, integer *nrhs,
	doublereal *a, integer *lda, doublereal *af, integer *ldaf, integer *
	ipiv, doublereal *b, integer *ldb, doublereal *x, integer *ldx,
	doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork,
	integer *info);

/* Subroutine */ extern int dsysv_(char *uplo, integer *n, integer *nrhs, doublereal
	*a, integer *lda, integer *ipiv, doublereal *b, integer *ldb,
	doublereal *work, integer *lwork, integer *info);

/* Subroutine */ extern int dsysvx_(char *fact, char *uplo, integer *n, integer *
	nrhs, doublereal *a, integer *lda, doublereal *af, integer *ldaf,
	integer *ipiv, doublereal *b, integer *ldb, doublereal *x, integer *
	ldx, doublereal *rcond, doublereal *ferr, doublereal *berr,
	doublereal *work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ extern int dsytd2_(char *uplo, integer *n, doublereal *a, integer *
	lda, doublereal *d__, doublereal *e, doublereal *tau, integer *info);

/* Subroutine */ extern int dsytf2_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *ipiv, integer *info);

/* Subroutine */ extern int dsytrd_(char *uplo, integer *n, doublereal *a, integer *
	lda, doublereal *d__, doublereal *e, doublereal *tau, doublereal *
	work, integer *lwork, integer *info);

/* Subroutine */ extern int dsytrf_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *ipiv, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ extern int dsytri_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *ipiv, doublereal *work, integer *info);

/* Subroutine */ extern int dsytrs_(char *uplo, integer *n, integer *nrhs,
	doublereal *a, integer *lda, integer *ipiv, doublereal *b, integer *
	ldb, integer *info);

/* Subroutine */ extern int dtbcon_(char *norm, char *uplo, char *diag, integer *n,
	integer *kd, doublereal *ab, integer *ldab, doublereal *rcond,
	doublereal *work, integer *iwork, integer *info);

/* Subroutine */ extern int dtbrfs_(char *uplo, char *trans, char *diag, integer *n,
	integer *kd, integer *nrhs, doublereal *ab, integer *ldab, doublereal
	*b, integer *ldb, doublereal *x, integer *ldx, doublereal *ferr,
	doublereal *berr, doublereal *work, integer *iwork, integer *info);

/* Subroutine */ extern int dtbtrs_(char *uplo, char *trans, char *diag, integer *n,
	integer *kd, integer *nrhs, doublereal *ab, integer *ldab, doublereal
	*b, integer *ldb, integer *info);

/* Subroutine */ extern int dtgevc_(char *side, char *howmny, logical *select,
	integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb,
	doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr, integer
	*mm, integer *m, doublereal *work, integer *info);

/* Subroutine */ extern int dtgex2_(logical *wantq, logical *wantz, integer *n,
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	q, integer *ldq, doublereal *z__, integer *ldz, integer *j1, integer *
	n1, integer *n2, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ extern int dtgexc_(logical *wantq, logical *wantz, integer *n,
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	q, integer *ldq, doublereal *z__, integer *ldz, integer *ifst,
	integer *ilst, doublereal *work, integer *lwork, integer *info);

/* Subroutine */ extern int dtgsen_(integer *ijob, logical *wantq, logical *wantz,
	logical *select, integer *n, doublereal *a, integer *lda, doublereal *
	b, integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *
	beta, doublereal *q, integer *ldq, doublereal *z__, integer *ldz,
	integer *m, doublereal *pl, doublereal *pr, doublereal *dif,
	doublereal *work, integer *lwork, integer *iwork, integer *liwork,
	integer *info);

/* Subroutine */ extern int dtgsja_(char *jobu, char *jobv, char *jobq, integer *m,
	integer *p, integer *n, integer *k, integer *l, doublereal *a,
	integer *lda, doublereal *b, integer *ldb, doublereal *tola,
	doublereal *tolb, doublereal *alpha, doublereal *beta, doublereal *u,
	integer *ldu, doublereal *v, integer *ldv, doublereal *q, integer *
	ldq, doublereal *work, integer *ncycle, integer *info);

/* Subroutine */ extern int dtgsna_(char *job, char *howmny, logical *select,
	integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb,
	doublereal *vl, integer *ldvl, doublereal *vr, integer *ldvr,
	doublereal *s, doublereal *dif, integer *mm, integer *m, doublereal *
	work, integer *lwork, integer *iwork, integer *info);

/* Subroutine */ extern int dtgsy2_(char *trans, integer *ijob, integer *m, integer *
	n, doublereal *a, integer *lda, doublereal *b, integer *ldb,
	doublereal *c__, integer *ldc, doublereal *d__, integer *ldd,
	doublereal *e, integer *lde, doublereal *f, integer *ldf, doublereal *
	scale, doublereal *rdsum, doublereal *rdscal, integer *iwork, integer
	*pq, integer *info);

/* Subroutine */ extern int dtgsyl_(char *trans, integer *ijob, integer *m, integer *
	n, doublereal *a, integer *lda, doublereal *b, integer *ldb,
	doublereal *c__, integer *ldc, doublereal *d__, integer *ldd,
	doublereal *e, integer *lde, doublereal *f, integer *ldf, doublereal *
	scale, doublereal *dif, doublereal *work, integer *lwork, integer *
	iwork, integer *info);

/* Subroutine */ extern int dtpcon_(char *norm, char *uplo, char *diag, integer *n,
	doublereal *ap, doublereal *rcond, doublereal *work, integer *iwork,
	integer *info);

/* Subroutine */ extern int dtprfs_(char *uplo, char *trans, char *diag, integer *n,
	integer *nrhs, doublereal *ap, doublereal *b, integer *ldb,
	doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr,
	doublereal *work, integer *iwork, integer *info);

/* Subroutine */ extern int dtptri_(char *uplo, char *diag, integer *n, doublereal *
	ap, integer *info);

/* Subroutine */ extern int dtptrs_(char *uplo, char *trans, char *diag, integer *n,
	integer *nrhs, doublereal *ap, doublereal *b, integer *ldb, integer *
	info);

/* Subroutine */ extern int dtrcon_(char *norm, char *uplo, char *diag, integer *n,
	doublereal *a, integer *lda, doublereal *rcond, doublereal *work,
	integer *iwork, integer *info);

/* Subroutine */ extern int dtrevc_(char *side, char *howmny, logical *select,
	integer *n, doublereal *t, integer *ldt, doublereal *vl, integer *
	ldvl, doublereal *vr, integer *ldvr, integer *mm, integer *m,
	doublereal *work, integer *info);

/* Subroutine */ extern int dtrexc_(char *compq, integer *n, doublereal *t, integer *
	ldt, doublereal *q, integer *ldq, integer *ifst, integer *ilst,
	doublereal *work, integer *info);

/* Subroutine */ extern int dtrrfs_(char *uplo, char *trans, char *diag, integer *n,
	integer *nrhs, doublereal *a, integer *lda, doublereal *b, integer *
	ldb, doublereal *x, integer *ldx, doublereal *ferr, doublereal *berr,
	doublereal *work, integer *iwork, integer *info);

/* Subroutine */ extern int dtrsen_(char *job, char *compq, logical *select, integer
	*n, doublereal *t, integer *ldt, doublereal *q, integer *ldq,
	doublereal *wr, doublereal *wi, integer *m, doublereal *s, doublereal
	*sep, doublereal *work, integer *lwork, integer *iwork, integer *
	liwork, integer *info);

/* Subroutine */ extern int dtrsna_(char *job, char *howmny, logical *select,
	integer *n, doublereal *t, integer *ldt, doublereal *vl, integer *
	ldvl, doublereal *vr, integer *ldvr, doublereal *s, doublereal *sep,
	integer *mm, integer *m, doublereal *work, integer *ldwork, integer *
	iwork, integer *info);

/* Subroutine */ extern int dtrsyl_(char *trana, char *tranb, integer *isgn, integer
	*m, integer *n, doublereal *a, integer *lda, doublereal *b, integer *
	ldb, doublereal *c__, integer *ldc, doublereal *scale, integer *info);

/* Subroutine */ extern int dtrti2_(char *uplo, char *diag, integer *n, doublereal *
	a, integer *lda, integer *info);

/* Subroutine */ extern int dtrtri_(char *uplo, char *diag, integer *n, doublereal *
	a, integer *lda, integer *info);

/* Subroutine */ extern int dtrtrs_(char *uplo, char *trans, char *diag, integer *n,
	integer *nrhs, doublereal *a, integer *lda, doublereal *b, integer *
	ldb, integer *info);

/* Subroutine */ extern int dtzrqf_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, integer *info);

/* Subroutine */ extern int dtzrzf_(integer *m, integer *n, doublereal *a, integer *
	lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);

integer ieeeck_(integer *ispec, real *zero, real *one);

integer ilaenv_(integer *ispec, char *name__, char *opts, integer *n1,
	integer *n2, integer *n3, integer *n4, ftnlen name_len, ftnlen
	opts_len);

/* Subroutine */ extern int xerbla_(char *srname, integer *info);

# endif
