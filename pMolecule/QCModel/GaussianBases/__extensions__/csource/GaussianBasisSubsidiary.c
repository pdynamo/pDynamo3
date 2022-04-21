/*==============================================================================
! . Subsidiary integral procedures.
!=============================================================================*/

# include <math.h>
# include <stdarg.h>
# include <stdlib.h>

# include "GaussianBasis.h"
# include "GaussianBasisSubsidiary.h"

/* # define CHECKGHPOINTS */
# ifdef CHECKGHPOINTS
# include <stdio.h>
# endif

/*----------------------------------------------------------------------------------------------------------------------------------
!
! . Notes:
!
!   All rectangular arrays used by these methods are compact. Thus, for a 4-D case with nI * nJ * nK * nL elements one has:
!
!   A[i,j,k,l] = A[ i * sI + j * sJ + k * sK + l * sL ]
!
!   with sL = 1, sK = nL * sL ; sJ = nK * sK ; sI = nJ * sJ.
!
!   Consecutive dimensions can be contracted. For example:
!
!   A[ij,kl] = A[ ij * sIJ + kl * sKL ]
!
!   with nIJ = nI * nJ, nKL = nK * nL, sIJ = sJ, sKL = 1.
!
!   In general, all methods assume the stride of the last dimension has a stride of 1. Where this is not the case, the stride is
!   explicitly specified.
!
!---------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Produce H anti-Coulomb integrals from G Coulomb integrals.
! . nI and nJ refer to H, those for G are nI+2 and nJ+2.
! . The dimensions for G and H are gStrideI and gStrideJ, and hStrideI and hStrideJ, respectively.
! . ( xIJ, yIJ, zIJ ) is the distance vector between the reference centers on I and J, not the
! . centers of I and J themselves.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _gStrideJ_ 1
# define _hStrideJ_ 1
void GaussianBasisSubsidiary_f1Ag1 ( const Integer  nI       ,
                                     const Integer  nJ       ,
                                     const Integer  gStrideI ,
                                     const Real    *Gx       ,
                                     const Real    *Gy       ,
                                     const Real    *Gz       ,
                                     const Real     xIJ      ,
                                     const Real     yIJ      ,
                                     const Real     zIJ      ,
                                     const Integer  hStrideI ,
                                           Real    *Hx       ,
                                           Real    *Hy       ,
                                           Real    *Hz       )
{
    Integer gI, gIJ, hI, hIJ, i, j ;
    Real    xIJ2, yIJ2, zIJ2 ;
    xIJ2 = xIJ * xIJ ;
    yIJ2 = yIJ * yIJ ;
    zIJ2 = zIJ * zIJ ;
    for ( i = 0 ; i <= nI ; i++ )
    {
        gI = i*gStrideI ;
        hI = i*hStrideI ;
        for ( j = 0 ; j <= nJ ; j++ )
        {
            gIJ = gI + j*_gStrideJ_ ;
            hIJ = hI + j*_hStrideJ_ ;
            Hx[hIJ] = xIJ2 * Gx[gIJ] + 2.0e+00 * xIJ * ( Gx[gIJ+gStrideI] - Gx[gIJ+_gStrideJ_] ) + Gx[gIJ+2*gStrideI] + Gx[gIJ+2*_gStrideJ_] - 2.0e+00 * Gx[gIJ+gStrideI+_gStrideJ_] ;
            Hy[hIJ] = yIJ2 * Gy[gIJ] + 2.0e+00 * yIJ * ( Gy[gIJ+gStrideI] - Gy[gIJ+_gStrideJ_] ) + Gy[gIJ+2*gStrideI] + Gy[gIJ+2*_gStrideJ_] - 2.0e+00 * Gy[gIJ+gStrideI+_gStrideJ_] ;
            Hz[hIJ] = zIJ2 * Gz[gIJ] + 2.0e+00 * zIJ * ( Gz[gIJ+gStrideI] - Gz[gIJ+_gStrideJ_] ) + Gz[gIJ+2*gStrideI] + Gz[gIJ+2*_gStrideJ_] - 2.0e+00 * Gz[gIJ+gStrideI+_gStrideJ_] ;
        }
    }
}
# undef _gStrideJ_
# undef _hStrideJ_

/*----------------------------------------------------------------------------------------------------------------------------------
! . Determine the 2-center nuclear subsidiary integrals.
! . The dimensions are strideI and strideJ = 1.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _strideJ_ 1
void GaussianBasisSubsidiary_f1Cg1 ( const Integer  nI      ,
                                     const Integer  nJ      ,
                                     const Real     b00     ,
                                     const Real     b10     ,
                                     const Real     bp01    ,
                                     const Real     f00     ,
                                     const Real     xc00    ,
                                     const Real     xcp00   ,
                                     const Real     yc00    ,
                                     const Real     ycp00   ,
                                     const Real     zc00    ,
                                     const Real     zcp00   ,
                                     const Integer  strideI ,
                                           Real    *Gx      ,
                                           Real    *Gy      ,
                                           Real    *Gz      )
{
    Integer i, j ;
    Real    cp01, cp10 = 0.0e+00, c01, c10 ;
    /* . (0,0). */
    Gx[0] = 1.0e+00 ;
    Gy[0] = 1.0e+00 ;
    Gz[0] = f00 ;
    /* . (1,0). */
    if ( nI > 0 )
    {
        Gx[strideI] = xc00 ;
        Gy[strideI] = yc00 ;
        Gz[strideI] = zc00 * f00 ;
    }
    /* . (0,1). */
    if ( nJ > 0 )
    {
        Gx[_strideJ_] = xcp00 ;
        Gy[_strideJ_] = ycp00 ;
        Gz[_strideJ_] = zcp00 * f00 ;
        /* . (1,1). */
        if ( nI > 0 )
        {
            cp10 = b00 ;
            Gx[strideI+_strideJ_] = xcp00 * Gx[strideI] + cp10 ;
            Gy[strideI+_strideJ_] = ycp00 * Gy[strideI] + cp10 ;
            Gz[strideI+_strideJ_] = zcp00 * Gz[strideI] + cp10 * f00 ;
        }
    }
    /* . nI > 1. */
    if ( nI > 1 )
    {
        c10 = 0.0e+00 ;
        for ( i = 2 ; i <= nI ; i++ )
        {
            /* . (i,0). */
            c10 += b10 ;
            Gx[i*strideI] = c10 * Gx[(i-2)*strideI] + xc00 * Gx[(i-1)*strideI] ;
            Gy[i*strideI] = c10 * Gy[(i-2)*strideI] + yc00 * Gy[(i-1)*strideI] ;
            Gz[i*strideI] = c10 * Gz[(i-2)*strideI] + zc00 * Gz[(i-1)*strideI] ;
            /* . (i,1). */
            if ( nJ > 0 )
            {
                cp10 += b00 ;
                Gx[i*strideI+_strideJ_] = xcp00 * Gx[i*strideI] + cp10 * Gx[(i-1)*strideI] ;
                Gy[i*strideI+_strideJ_] = ycp00 * Gy[i*strideI] + cp10 * Gy[(i-1)*strideI] ;
                Gz[i*strideI+_strideJ_] = zcp00 * Gz[i*strideI] + cp10 * Gz[(i-1)*strideI] ;
            }
        }
    }
    /* . nJ > 1. */
    if ( nJ > 1 )
    {
        cp01 = 0.0e+00 ;
        c01  = b00 ;
        for ( j = 2 ; j <= nJ ; j++ )
        {
            /* . (0,j). */
            cp01 += bp01 ;
            Gx[j*_strideJ_] = cp01 * Gx[(j-2)*_strideJ_] + xcp00 * Gx[(j-1)*_strideJ_] ;
            Gy[j*_strideJ_] = cp01 * Gy[(j-2)*_strideJ_] + ycp00 * Gy[(j-1)*_strideJ_] ;
            Gz[j*_strideJ_] = cp01 * Gz[(j-2)*_strideJ_] + zcp00 * Gz[(j-1)*_strideJ_] ;
            /* . (1,j). */
            if ( nI > 0 )
            {
                c01 += b00 ;
                Gx[j*_strideJ_+strideI] = xc00 * Gx[j*_strideJ_] + c01 * Gx[(j-1)*_strideJ_] ;
                Gy[j*_strideJ_+strideI] = yc00 * Gy[j*_strideJ_] + c01 * Gy[(j-1)*_strideJ_] ;
                Gz[j*_strideJ_+strideI] = zc00 * Gz[j*_strideJ_] + c01 * Gz[(j-1)*_strideJ_] ;
            }
        }
    }
    /* . nI and nJ > 1. */
    if ( ( nI > 1 ) && ( nJ > 1 ) )
    {
        /* . (j,i). */
        c01 = b00 ;
        for ( j = 2 ; j <= nJ ; j++ )
        {
            c01 += b00 ;
            c10  = b10 ;
            for ( i = 2 ; i <= nI ; i++ )
            {
               Gx[j*_strideJ_+i*strideI] = c10 * Gx[j*_strideJ_+(i-2)*strideI] + xc00 * Gx[j*_strideJ_+(i-1)*strideI] + c01 * Gx[(j-1)*_strideJ_+(i-1)*strideI] ;
               Gy[j*_strideJ_+i*strideI] = c10 * Gy[j*_strideJ_+(i-2)*strideI] + yc00 * Gy[j*_strideJ_+(i-1)*strideI] + c01 * Gy[(j-1)*_strideJ_+(i-1)*strideI] ;
               Gz[j*_strideJ_+i*strideI] = c10 * Gz[j*_strideJ_+(i-2)*strideI] + zc00 * Gz[j*_strideJ_+(i-1)*strideI] + c01 * Gz[(j-1)*_strideJ_+(i-1)*strideI] ;
               c10 += b10 ;
            }
        }
    }
}
# undef _strideJ_

/*----------------------------------------------------------------------------------------------------------------------------------
! . Dipole integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisSubsidiary_f1Dg1 (       Real    *x      ,
                                           Real    *y      ,
                                           Real    *z      ,
                                     const Real    aa      ,
                                     const Real    *r0     ,
                                     const Real    *ri     ,
                                     const Real    *rj     ,
                                     const Real    *center ,
                                     const Integer  ni     ,
                                     const Integer  nj     )
{
   Integer i, j, k, n, npts, p, q ;
   Real    ar, br, ptr, rint[3], t, temp, tinv ;
   t    = sqrt ( aa ) ;
   tinv = 1.0e+00 / t ;
   for ( i = 0, n = 0 ; i <= ni ; i++ )
   {
      for ( j = 0 ; j <= nj ; j++, n++ )
      {
	 for ( k = 0 ; k < 3 ; k++ ) rint[k] = 0.0e+00 ;
	 npts = ( i + j + 1 ) / 2 ;
# ifdef CHECKGHPOINTS
if ( npts > GHMAXPT ) printf ( "\nINVALID NUMBER OF POINTS IN GAUSS-HERMITE QUADRATURE = %d\n", npts ) ;
# endif
	 for ( p = GHINDEX[npts] ; p < GHINDEX[npts+1] ; p++ )
         {
	    for ( k = 0 ; k < 3 ; k++ )
            {
	       ptr  = GHABSCISSAE[p] / t + r0[k] ;
               ar   = ptr - ri[k] ;
               br   = ptr - rj[k] ;
               temp = GHWEIGHTS[p] * ( ptr - center[k] ) ;
	       for ( q = 1 ; q <= i ; q++ ) temp *= ar ;
	       for ( q = 1 ; q <= j ; q++ ) temp *= br ;
	       rint[k] += temp ;
            }
	 }
	 x[n] = rint[0] * tinv ;
	 y[n] = rint[1] * tinv ;
	 z[n] = rint[2] * tinv ;
      }
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Determine the kinetic energy subsidiary integrals.
! . The jth dimension of x, y and z is assumed to be nj+3.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisSubsidiary_f1Kg1 ( const Real    *x     ,
                                     const Real    *y     ,
                                     const Real    *z     ,
                                           Real    *xt    ,
                                           Real    *yt    ,
                                           Real    *zt    ,
                                     const Real     aj    ,
                                     const Integer  ni    ,
                                     const Integer  nj    ,
                                     const Integer  jdimo ,
                                     const Integer  jdimt )
{
    Integer i, io, it, j ;
    Real    a2 ;
    a2 = aj + aj ;
    for ( i = 0, io = 0, it = 0 ; i <= ni ; i++, io += jdimo, it += jdimt )
    {
        xt[it] = ( x[io] - x[io+2] * a2 ) * aj ;
        yt[it] = ( y[io] - y[io+2] * a2 ) * aj ;
        zt[it] = ( z[io] - z[io+2] * a2 ) * aj ;
    }
    if ( nj > 0 )
    {
        for ( i = 0, io = 0, it = 0 ; i <= ni ; i++, io += jdimo, it += jdimt )
        {
            xt[it+1] = ( x[io+1] * 3.0e+00 - x[io+3] * a2 ) * aj ;
            yt[it+1] = ( y[io+1] * 3.0e+00 - y[io+3] * a2 ) * aj ;
            zt[it+1] = ( z[io+1] * 3.0e+00 - z[io+3] * a2 ) * aj ;
        }
        if ( nj > 1 )
        {
            for ( j = 2 ; j <= nj ; j++ )
            {
                for ( i = 0, io = 0, it = 0 ; i <= ni ; i++, io += jdimo, it += jdimt )
                {
	            xt[it+j] = ( x[io+j] * ( Real ) ( 2*j+1 ) - x[io+j+2] * a2 ) * aj - x[io+j-2] * ( Real ) ( (j*(j-1))/2 ) ;
	            yt[it+j] = ( y[io+j] * ( Real ) ( 2*j+1 ) - y[io+j+2] * a2 ) * aj - y[io+j-2] * ( Real ) ( (j*(j-1))/2 ) ;
	            zt[it+j] = ( z[io+j] * ( Real ) ( 2*j+1 ) - z[io+j+2] * a2 ) * aj - z[io+j-2] * ( Real ) ( (j*(j-1))/2 ) ;
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Double overlap integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisSubsidiary_f1Og1 (       Real    *x  ,
                                           Real    *y  ,
                                           Real    *z  ,
                                     const Real     aa ,
                                     const Real    *r0 ,
                                     const Real    *ri ,
                                     const Real    *rj ,
                                     const Integer  ni ,
                                     const Integer  nj )
{
   Integer i, j, k, n, npts, p, q ;
   Real    ar, br, ptr, rint[3], t, temp, tinv ;
   t    = sqrt ( aa ) ;
   tinv = 1.0e+00 / t ;
   for ( i = 0, n = 0 ; i <= ni ; i++ )
   {
      for ( j = 0 ; j <= nj ; j++, n++ )
      {
	 for ( k = 0 ; k < 3 ; k++ ) rint[k] = 0.0e+00 ;
	 npts = ( i + j ) / 2 ;
# ifdef CHECKGHPOINTS
if ( npts > GHMAXPT ) printf ( "\nINVALID NUMBER OF POINTS IN GAUSS-HERMITE QUADRATURE = %d\n", npts ) ;
# endif
	 for ( p = GHINDEX[npts] ; p < GHINDEX[npts+1] ; p++ )
         {
	    for ( k = 0 ; k < 3 ; k++ )
            {
	       ptr  = GHABSCISSAE[p] / t + r0[k] ;
               ar   = ptr - ri[k] ;
               br   = ptr - rj[k] ;
               temp = GHWEIGHTS[p] ;
	       for ( q = 1 ; q <= i ; q++ ) temp *= ar ;
	       for ( q = 1 ; q <= j ; q++ ) temp *= br ;
	       rint[k] += temp ;
            }
	 }
	 x[n] = rint[0] * tinv ;
	 y[n] = rint[1] * tinv ;
	 z[n] = rint[2] * tinv ;
      }
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Triple overlap integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisSubsidiary_f1Og2 (       Real    *x  ,
                                           Real    *y  ,
                                           Real    *z  ,
                                     const Real     aa ,
                                     const Real    *r0 ,
                                     const Real    *ri ,
                                     const Real    *rj ,
                                     const Real    *rk ,
                                     const Integer  ni ,
                                     const Integer  nj ,
                                     const Integer  nk )
{
    Integer i, j, k, m, n, npts, p, q ;
    Real    ar, br, cr, ptr, rint[3], t, temp, tinv ;
    t    = sqrt ( aa ) ;
    tinv = 1.0e+00 / t ;
    for ( i = 0, n = 0 ; i <= ni ; i++ )
    {
        for ( j = 0 ; j <= nj ; j++ )
        {
	    for ( k = 0 ; k <= nk ; k++, n++ )
            {
                for ( m = 0 ; m < 3 ; m++ ) rint[m] = 0.0e+00 ;
	        npts = ( i + j + k ) / 2 ;
# ifdef CHECKGHPOINTS
if ( npts > GHMAXPT ) printf ( "\nINVALID NUMBER OF POINTS IN GAUSS-HERMITE QUADRATURE = %d\n", npts ) ;
# endif
	        for ( p = GHINDEX[npts] ; p < GHINDEX[npts+1] ; p++ )
                {
	            for ( m = 0 ; m < 3 ; m++ )
                    {
	                ptr  = GHABSCISSAE[p] / t + r0[m] ;
                        ar   = ptr - ri[m] ;
                        br   = ptr - rj[m] ;
                        cr   = ptr - rk[m] ;
                        temp = GHWEIGHTS[p] ;
	                for ( q = 1 ; q <= i ; q++ ) temp *= ar ;
	                for ( q = 1 ; q <= j ; q++ ) temp *= br ;
	                for ( q = 1 ; q <= k ; q++ ) temp *= cr ;
	                rint[m] += temp ;
                    }
	        }
	        x[n] = rint[0] * tinv ;
	        y[n] = rint[1] * tinv ;
	        z[n] = rint[2] * tinv ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Quadrupole integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisSubsidiary_f1Qg1 (       Real    *x      ,
                                           Real    *y      ,
                                           Real    *z      ,
                                     const Real    aa      ,
                                     const Real    *r0     ,
                                     const Real    *ri     ,
                                     const Real    *rj     ,
                                     const Real    *center ,
                                     const Integer  ni     ,
                                     const Integer  nj     )
{
   Integer i, j, k, n, npts, p, q ;
   Real    ar, br, ptr, rint[3], t, temp, tinv ;
   t    = sqrt ( aa ) ;
   tinv = 1.0e+00 / t ;
   for ( i = 0, n = 0 ; i <= ni ; i++ )
   {
      for ( j = 0 ; j <= nj ; j++, n++ )
      {
	 for ( k = 0 ; k < 3 ; k++ ) rint[k] = 0.0e+00 ;
	 npts = ( i + j + 2 ) / 2 ;
# ifdef CHECKGHPOINTS
if ( npts > GHMAXPT ) printf ( "\nINVALID NUMBER OF POINTS IN GAUSS-HERMITE QUADRATURE = %d\n", npts ) ;
# endif
	 for ( p = GHINDEX[npts] ; p < GHINDEX[npts+1] ; p++ )
         {
	    for ( k = 0 ; k < 3 ; k++ )
            {
	       ptr  = GHABSCISSAE[p] / t + r0[k] ;
               ar   = ptr - ri[k] ;
               br   = ptr - rj[k] ;
               temp = GHWEIGHTS[p] * ( ptr - center[k] ) * ( ptr - center[k] ) ;
	       for ( q = 1 ; q <= i ; q++ ) temp *= ar ;
	       for ( q = 1 ; q <= j ; q++ ) temp *= br ;
	       rint[k] += temp ;
            }
	 }
	 x[n] = rint[0] * tinv ;
	 y[n] = rint[1] * tinv ;
	 z[n] = rint[2] * tinv ;
      }
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Determine derivative integrals from input integrals.
! . jdim is the jth-dimension of both x, y, z and xd, yd, zd.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisSubsidiary_f1Xg1r ( const Real    *Gx       ,
                                      const Real    *Gy       ,
                                      const Real    *Gz       ,
                                      const Real     aI       ,
                                      const Integer  nI       ,
                                      const Integer  nJ       ,
                                      const Integer  gStrideI ,
                                      const Integer  dStrideI ,
                                            Real    *GxD      ,
                                            Real    *GyD      ,
                                            Real    *GzD      )
{
   Integer i, ijD, ijG, j ;
   Real    a2 ;
   a2 = aI + aI ;
   for ( j = 0 ; j <= nJ ; j++ )
   {
      /* . (0,j). */
      GxD[j] = a2 * Gx[gStrideI+j] ;
      GyD[j] = a2 * Gy[gStrideI+j] ;
      GzD[j] = a2 * Gz[gStrideI+j] ;
      /* . (i,j). */
      for ( i = 1 ; i <= nI ; i++ )
      {
         ijD      = i * dStrideI + j;
         ijG      = i * gStrideI + j;
	 GxD[ijD] = a2 * Gx[ijG+gStrideI] - ( ( Real ) i ) * Gx[ijG-gStrideI] ;
	 GyD[ijD] = a2 * Gy[ijG+gStrideI] - ( ( Real ) i ) * Gy[ijG-gStrideI] ;
	 GzD[ijD] = a2 * Gz[ijG+gStrideI] - ( ( Real ) i ) * Gz[ijG-gStrideI] ;
      }
   }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Generate g2 integrals from g1 integrals by applying the transfer relations.
! . The input integrals are G[NI*NJ,NK] whereas the output integrals are T[NI,NJ,NK].
! . The Ts need to be initialized on entry.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisSubsidiary_f1Xg2i ( const Integer  nI        ,
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
                                            Real    *Tz        )
{
    Integer b, g, i, j, k, t ;
    Real    c, cX, cY, cZ ;
    for ( j = 0 ; j <= nJ ; j++ )
    {
        for ( b = 0 ; b <= j ; b++ )
        {
            if ( ( b == 0 ) || ( j == 0 ) )
            {
                cX = cY = cZ = 1.0 ;
            }
            else
            {
                c   = ( Real ) ( j - b + 1 ) / ( Real ) b ;
                cX *= ( xIJ * c ) ;
                cY *= ( yIJ * c ) ;
                cZ *= ( zIJ * c ) ; 
            }
            for ( i = 0 ; i <= nI ; i++ )
            {
                g = (i+j-b)*gStrideIJ ;
                t = i*tStrideI+j*tStrideJ ;
                for ( k = 0 ; k <= nK ; k++ )
                {
                    Tx[t] += cX * Gx[g] ;
                    Ty[t] += cY * Gy[g] ;
                    Tz[t] += cZ * Gz[g] ;
                    g     += gStrideK ;
                    t     += tStrideK ;
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Determine derivative integrals from input integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisSubsidiary_f1Xg2r ( const Real    *x     ,
                                      const Real    *y     ,
                                      const Real    *z     ,
                                            Real    *xg    ,
                                            Real    *yg    ,
                                            Real    *zg    ,
                                            Real    *xh    ,
                                            Real    *yh    ,
                                            Real    *zh    ,
                                      const Real     ag    ,
                                      const Real     ah    ,
                                      const Integer  ni    ,
                                      const Integer  nj    ,
                                      const Integer  nf    ,
                                      const Integer  dim1  ,
                                      const Integer  dim2  ,
                                      const Integer  ddim1 ,
                                      const Integer  ddim2 )
{
    Integer f, i, j ;
    Real    ag2, ah2 ;
    /* . Initialization. */
    ag2 = ag + ag ;
    ah2 = ah + ah ;
    /* . Loop over the values of nf. */
    for ( f = 0 ; f <= nf ; f++ )
    {
        /* . Loop over the values of nj. */
        for ( j = 0 ; j <= nj ; j++ )
        {
            /* . i = 0. */
            xg[f+j*ddim1] = ag2 * x[f+j*dim1+dim2] ;
            yg[f+j*ddim1] = ag2 * y[f+j*dim1+dim2] ;
            zg[f+j*ddim1] = ag2 * z[f+j*dim1+dim2] ;
            /* . Loop over the remaining values of ni. */
            for ( i = 1 ; i <= ni ; i++ )
            {
	        xg[f+j*ddim1+i*ddim2] = ag2 * x[f+j*dim1+(i+1)*dim2] - ( Real ) i * x[f+j*dim1+(i-1)*dim2] ;
	        yg[f+j*ddim1+i*ddim2] = ag2 * y[f+j*dim1+(i+1)*dim2] - ( Real ) i * y[f+j*dim1+(i-1)*dim2] ;
	        zg[f+j*ddim1+i*ddim2] = ag2 * z[f+j*dim1+(i+1)*dim2] - ( Real ) i * z[f+j*dim1+(i-1)*dim2] ;
            }
        }
        /* . Loop over the values of ni. */
        for ( i = 0 ; i <= ni ; i++ )
        {
            /* . j = 0.*/
            xh[f+i*ddim2] = ah2 * x[f+dim1+i*dim2] ;
            yh[f+i*ddim2] = ah2 * y[f+dim1+i*dim2] ;
            zh[f+i*ddim2] = ah2 * z[f+dim1+i*dim2] ;
            /* . Loop over the remaining values of nj. */
            for ( j = 1 ; j <= nj ; j++ )
            {
	        xh[f+j*ddim1+i*ddim2] = ah2 * x[f+(j+1)*dim1+i*dim2] - ( Real ) j * x[f+(j-1)*dim1+i*dim2] ;
	        yh[f+j*ddim1+i*ddim2] = ah2 * y[f+(j+1)*dim1+i*dim2] - ( Real ) j * y[f+(j-1)*dim1+i*dim2] ;
	        zh[f+j*ddim1+i*ddim2] = ah2 * z[f+(j+1)*dim1+i*dim2] - ( Real ) j * z[f+(j-1)*dim1+i*dim2] ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Determine derivative integrals from input integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisSubsidiary_f2Xg2r ( const Integer  nI       , 
                                      const Integer  nJ       , 
                                      const Integer  nK       , 
                                      const Integer  nL       , 
                                      const Integer  strideI  ,
                                      const Integer  strideJ  ,
                                      const Integer  strideK  ,
                                      const Integer  strideL  ,
                                      const Integer  dStrideI ,
                                      const Integer  dStrideJ ,
                                      const Integer  dStrideK ,
                                      const Integer  dStrideL ,
                                      const Real     aI       ,
                                      const Real     aJ       ,
                                      const Real     aK       ,
                                      const Real    *x        ,
                                      const Real    *y        ,
                                      const Real    *z        ,
                                            Real    *dXi      ,
                                            Real    *dYi      ,
                                            Real    *dZi      ,
                                            Real    *dXj      ,
                                            Real    *dYj      ,
                                            Real    *dZj      ,
                                            Real    *dXk      ,
                                            Real    *dYk      ,
                                            Real    *dZk      )
{
    Integer i, j, k, l ;
    Real    aI2, aJ2, aK2 ;
    /* . Initialization. */
    aI2 = aI + aI ;
    aJ2 = aJ + aJ ;
    aK2 = aK + aK ;
    /* . Loop over all angular momentum values, n. */
    /* . All n = 0 cases are special. */
    for ( l = 0 ; l <= nL ; l++ )
    {
        for ( k = 0 ; k <= nK ; k++ )
        {
            for ( j = 0 ; j <= nJ ; j++ )
            {
                dXi[l*dStrideL+k*dStrideK+j*dStrideJ] = aI2 * x[l*strideL+k*strideK+j*strideJ+strideI] ;
                dYi[l*dStrideL+k*dStrideK+j*dStrideJ] = aI2 * y[l*strideL+k*strideK+j*strideJ+strideI] ;
                dZi[l*dStrideL+k*dStrideK+j*dStrideJ] = aI2 * z[l*strideL+k*strideK+j*strideJ+strideI] ;
                for ( i = 1 ; i <= nI ; i++ )
                {
                    dXi[l*dStrideL+k*dStrideK+j*dStrideJ+i*dStrideI] = aI2 * x[l*strideL+k*strideK+j*strideJ+(i+1)*strideI] - ( Real ) i * x[l*strideL+k*strideK+j*strideJ+(i-1)*strideI] ;
                    dYi[l*dStrideL+k*dStrideK+j*dStrideJ+i*dStrideI] = aI2 * y[l*strideL+k*strideK+j*strideJ+(i+1)*strideI] - ( Real ) i * y[l*strideL+k*strideK+j*strideJ+(i-1)*strideI] ;
                    dZi[l*dStrideL+k*dStrideK+j*dStrideJ+i*dStrideI] = aI2 * z[l*strideL+k*strideK+j*strideJ+(i+1)*strideI] - ( Real ) i * z[l*strideL+k*strideK+j*strideJ+(i-1)*strideI] ;
                }
            }
        }
        for ( i = 0 ; i <= nI ; i++ )
        {
            for ( k = 0 ; k <= nK ; k++ )
            {
                dXj[l*dStrideL+k*dStrideK+i*dStrideI] = aJ2 * x[l*strideL+k*strideK+strideJ+i*strideI] ;
                dYj[l*dStrideL+k*dStrideK+i*dStrideI] = aJ2 * y[l*strideL+k*strideK+strideJ+i*strideI] ;
                dZj[l*dStrideL+k*dStrideK+i*dStrideI] = aJ2 * z[l*strideL+k*strideK+strideJ+i*strideI] ;
                for ( j = 1 ; j <= nJ ; j++ )
                {
                    dXj[l*dStrideL+k*dStrideK+j*dStrideJ+i*dStrideI] = aJ2 * x[l*strideL+k*strideK+(j+1)*strideJ+i*strideI] - ( Real ) j * x[l*strideL+k*strideK+(j-1)*strideJ+i*strideI] ;
                    dYj[l*dStrideL+k*dStrideK+j*dStrideJ+i*dStrideI] = aJ2 * y[l*strideL+k*strideK+(j+1)*strideJ+i*strideI] - ( Real ) j * y[l*strideL+k*strideK+(j-1)*strideJ+i*strideI] ;
                    dZj[l*dStrideL+k*dStrideK+j*dStrideJ+i*dStrideI] = aJ2 * z[l*strideL+k*strideK+(j+1)*strideJ+i*strideI] - ( Real ) j * z[l*strideL+k*strideK+(j-1)*strideJ+i*strideI] ;
                }
            }
        }
        for ( j = 0 ; j <= nJ ; j++ )
        {
            for ( i = 0 ; i <= nI ; i++ )
            {
                dXk[l*dStrideL+j*dStrideJ+i*dStrideI] = aK2 * x[l*strideL+strideK+j*strideJ+i*strideI] ;
                dYk[l*dStrideL+j*dStrideJ+i*dStrideI] = aK2 * y[l*strideL+strideK+j*strideJ+i*strideI] ;
                dZk[l*dStrideL+j*dStrideJ+i*dStrideI] = aK2 * z[l*strideL+strideK+j*strideJ+i*strideI] ;
                for ( k = 1 ; k <= nK ; k++ )
                {
                    dXk[l*dStrideL+k*dStrideK+j*dStrideJ+i*dStrideI] = aK2 * x[l*strideL+(k+1)*strideK+j*strideJ+i*strideI] - ( Real ) k * x[l*strideL+(k-1)*strideK+j*strideJ+i*strideI] ;
                    dYk[l*dStrideL+k*dStrideK+j*dStrideJ+i*dStrideI] = aK2 * y[l*strideL+(k+1)*strideK+j*strideJ+i*strideI] - ( Real ) k * y[l*strideL+(k-1)*strideK+j*strideJ+i*strideI] ;
                    dZk[l*dStrideL+k*dStrideK+j*dStrideJ+i*dStrideI] = aK2 * z[l*strideL+(k+1)*strideK+j*strideJ+i*strideI] - ( Real ) k * z[l*strideL+(k-1)*strideK+j*strideJ+i*strideI] ;
                }
            }
        }
    }
}
