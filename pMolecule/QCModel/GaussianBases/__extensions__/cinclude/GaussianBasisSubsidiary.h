# ifndef _GAUSSIANBASISSUBSIDIARY
# define _GAUSSIANBASISSUBSIDIARY

# include "Boolean.h"
# include "Integer.h"
# include "NumericalMacros.h"
# include "Real.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . 1-D subsidiary integral macros. */
/* . Overlap. */
# define _1Overlap( o, aI, nMaximum ) \
{ \
    auto Integer n ; \
    o[0] = PI12 / sqrt ( aI ) ; \
    for ( n = 1 ; n <= nMaximum ; n++ ) { if ( IsEven ( n ) ) o[n] = ( Real ) ( n - 1 ) * o[n-2] / ( 2.0e+00 * aI ) ; else o[n] = 0.0e+00 ; } \
}
/* . Derivatives. Higher multipole integrals, m1, are calculated in terms of m0 and r = function center - multipole center. */
# define _1Derivative( m1, m0, r, nMaximum ) \
{ \
    auto Integer n ; \
    for ( n = 0 ; n <= nMaximum ; n++ ) { m1[n] = m0[n+1] + r * m0[n] ; } \
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void GaussianBasisSubsidiary_f1Ag1  ( const Integer  nI        ,
                                             const Integer  nJ        ,
                                             const Integer  gStrideI  ,
                                             const Real    *Gx        ,
                                             const Real    *Gy        ,
                                             const Real    *Gz        ,
                                             const Real     xIJ       ,
                                             const Real     yIJ       ,
                                             const Real     zIJ       ,
                                             const Integer  hStrideI  ,
                                                   Real    *Hx        ,
                                                   Real    *Hy        ,
                                                   Real    *Hz        ) ;
extern void GaussianBasisSubsidiary_f1Cg1  ( const Integer  nI        ,
                                             const Integer  nJ        ,
                                             const Real     b00       ,
                                             const Real     b10       ,
                                             const Real     bp01      ,
                                             const Real     f00       ,
                                             const Real     xc00      ,
                                             const Real     xcp00     ,
                                             const Real     yc00      ,
                                             const Real     ycp00     ,
                                             const Real     zc00      ,
                                             const Real     zcp00     ,
                                             const Integer  strideI   ,
                                                   Real    *Gx        ,
                                                   Real    *Gy        ,
                                                   Real    *Gz        ) ;
extern void GaussianBasisSubsidiary_f1Dg1  (       Real    *x         ,
                                                   Real    *y         ,
                                                   Real    *z         ,
                                             const Real    aa         ,
                                             const Real    *r0        ,
                                             const Real    *ri        ,
                                             const Real    *rj        ,
                                             const Real    *center    ,
                                             const Integer  ni        ,
                                             const Integer  nj        ) ;
extern void GaussianBasisSubsidiary_f1Kg1  ( const Real    *x         ,
                                             const Real    *y         ,
                                             const Real    *z         ,
                                                   Real    *xt        ,
                                                   Real    *yt        ,
                                                   Real    *zt        ,
                                             const Real     aj        ,
                                             const Integer  ni        ,
                                             const Integer  nj        ,
                                             const Integer  jdimo     ,
                                             const Integer  jdimt     ) ;
extern void GaussianBasisSubsidiary_f1Og1  (       Real    *x         ,
                                                   Real    *y         ,
                                                   Real    *z         ,
                                             const Real    aa         ,
                                             const Real    *r0        ,
                                             const Real    *ri        ,
                                             const Real    *rj        ,
                                             const Integer  ni        ,
                                             const Integer  nj        ) ;
extern void GaussianBasisSubsidiary_f1Og2  (       Real    *x         ,
                                                   Real    *y         ,
                                                   Real    *z         ,
                                             const Real    aa         ,
                                             const Real    *r0        ,
                                             const Real    *ri        ,
                                             const Real    *rj        ,
                                             const Real    *rk        ,
                                             const Integer  ni        ,
                                             const Integer  nj        ,
                                             const Integer  nk        ) ;
extern void GaussianBasisSubsidiary_f1Qg1  (       Real    *x         ,
                                                   Real    *y         ,
                                                   Real    *z         ,
                                             const Real    aa         ,
                                             const Real    *r0        ,
                                             const Real    *ri        ,
                                             const Real    *rj        ,
                                             const Real    *center    ,
                                             const Integer  ni        ,
                                             const Integer  nj        ) ;
extern void GaussianBasisSubsidiary_f1Xg1r ( const Real    *Gx        ,
                                             const Real    *Gy        ,
                                             const Real    *Gz        ,
                                             const Real     aI        ,
                                             const Integer  nI        ,
                                             const Integer  nJ        ,
                                             const Integer  gStrideI  ,
                                             const Integer  dStrideI  ,
                                                   Real    *GxD       ,
                                                   Real    *GyD       ,
                                                   Real    *GzD       ) ;
extern void GaussianBasisSubsidiary_f1Xg2i ( const Integer  nI        ,
                                             const Integer  nJ        ,
                                             const Integer  nK        ,
                                             const Integer  gStrideIJ ,
                                             const Integer  gStrideK  ,
                                             const Real    *Gx        ,
                                             const Real    *Gy        ,
                                             const Real    *Gz        ,
                                             const Real     xIJ       ,
                                             const Real     yIJ       ,
                                             const Real     zIJ       ,
                                             const Integer  tStrideI  ,
                                             const Integer  tStrideJ  ,
                                             const Integer  tStrideK  ,
                                                   Real    *Tx        ,
                                                   Real    *Ty        ,
                                                   Real    *Tz        ) ;
extern void GaussianBasisSubsidiary_f1Xg2r ( const Real    *x         ,
                                             const Real    *y         ,
                                             const Real    *z         ,
                                                   Real    *xg        ,
                                                   Real    *yg        ,
                                                   Real    *zg        ,
                                                   Real    *xh        ,
                                                   Real    *yh        ,
                                                   Real    *zh        ,
                                             const Real     ag        ,
                                             const Real     ah        ,
                                             const Integer  ni        ,
                                             const Integer  nj        ,
                                             const Integer  nf        ,
                                             const Integer  dim1      ,
                                             const Integer  dim2      ,
                                             const Integer  ddim1     ,
                                             const Integer  ddim2     ) ;
extern void GaussianBasisSubsidiary_f2Xg2r ( const Integer  nI        , 
                                             const Integer  nJ        , 
                                             const Integer  nK        , 
                                             const Integer  nL        , 
                                             const Integer  strideI   ,
                                             const Integer  strideJ   ,
                                             const Integer  strideK   ,
                                             const Integer  strideL   ,
                                             const Integer  dStrideI  ,
                                             const Integer  dStrideJ  ,
                                             const Integer  dStrideK  ,
                                             const Integer  dStrideL  ,
                                             const Real     aI        , 
                                             const Real     aJ        , 
                                             const Real     aK        ,
                                             const Real    *x         ,
                                             const Real    *y         ,
                                             const Real    *z         ,
                                                   Real    *dXi       ,
                                                   Real    *dYi       ,
                                                   Real    *dZi       ,
                                                   Real    *dXj       ,
                                                   Real    *dYj       ,
                                                   Real    *dZj       ,
                                                   Real    *dXk       ,
                                                   Real    *dYk       ,
                                                   Real    *dZk       ) ;
# endif
