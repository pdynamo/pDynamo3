/*==================================================================================================================================
! . Integrals - 2 bases, 1 electron, 0 nuclei/points.
!=================================================================================================================================*/

# include <math.h>
# include <stdarg.h>
# include <stdlib.h>

# include "Boolean.h"
# include "GaussianBasisIntegrals_b2e1n0.h"
# include "GaussianBasisSubsidiary.h"
# include "Integer.h"
# include "RysQuadrature.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Coulomb integrals.
! . integrals is overwritten by this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisIntegrals_2Coulomb ( const GaussianBasis *iBasis    ,
                                       const Real          *rI        ,
                                       const GaussianBasis *jBasis    ,
                                       const Real          *rJ        ,
                                             RealArray2D   *integrals )
{
    Boolean       iIsJ ;
    Integer       i, iammax, icbfind, ip, iShell, ix, iy, iz, j, jammax, jcbfind, jdim, jdimm,
                  jp, jShell, jUpper, jxix, jyiy, jziz, m, n, ncfunci, ncfuncj, nroots ;
    Real          aandb, ab, ai, aj, bp01, b00, b10, c1x, c1y, c1z, c3x, c3y, c3z, dfi, dfIJ,
                  fac, fac2, f00, rho, rIJ2, ti, tIJ, u2, xc00, xcp00, xIJ, yc00, ycp00, yIJ, zc00, zcp00, zIJ ;
    Real          g[MAXCBF*MAXCBF],
                  xint[MAXAMP1*MAXAMP1*_MAXRYS], yint[MAXAMP1*MAXAMP1*_MAXRYS], zint[MAXAMP1*MAXAMP1*_MAXRYS] ;
    RysQuadrature roots ;
    /* . Initialization. */
    iIsJ = ( iBasis == jBasis ) && ( rI == rJ ) ;
    xIJ  = rI[0] - rJ[0] ;
    yIJ  = rI[1] - rJ[1] ;
    zIJ  = rI[2] - rJ[2] ;
    rIJ2 = xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ;
    RealArray2D_Set ( integrals, 0.0e+00 ) ;
    /* . Outer loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nshells ; iShell++ )
    {
        iammax  = iBasis->shells[iShell].type->angularmomentum_high ;
        icbfind = iBasis->shells[iShell].type->cbfindex ;
        ncfunci = iBasis->shells[iShell].type->ncbf     ;
        /* . Inner loop over shells. */
        if ( iIsJ ) jUpper = iShell + 1 ;
        else        jUpper = jBasis->nshells ;
        for ( jShell = 0 ; jShell < jUpper ; jShell++ )
        {
            jammax  = jBasis->shells[jShell].type->angularmomentum_high ;
            jdim    = jammax + 1 ;
            jdimm   = ( iammax + 1 ) * ( jammax + 1 ) ;
            jcbfind = jBasis->shells[jShell].type->cbfindex ;
            ncfuncj = jBasis->shells[jShell].type->ncbf     ;
            /* . Get the number of roots. */
            nroots = ( iammax + jammax ) / 2 + 1 ;
            /* . Initialize the integral blocks. */
            for ( i = 0 ; i < ( ncfunci * ncfuncj ) ; i++ ) g[i] = 0.0e+00 ;
            /* . Outer loop over primitives. */
            for ( ip = 0 ; ip < iBasis->shells[iShell].nprimitives ; ip++ )
            {
                /* . Get some information for the primitive. */
	        ai  = iBasis->shells[iShell].primitives[ip].exponent ;
                dfi = PI252 / ai ;
                /* . Inner loop over primitives. */
                for ( jp = 0 ; jp < jBasis->shells[jShell].nprimitives ; jp++ )
                {
                    /* . Get some information for the primitive. */
	            aj = jBasis->shells[jShell].primitives[jp].exponent ;
                    /* . Calculate some factors. */
                    ab    = ai * aj ;
                    aandb = ai + aj ;
                    rho   = ab / aandb ;
                    dfIJ  = dfi / ( aj * sqrt ( aandb ) ) ;
                    /* . Calculate some displacements. */
                    c1x  = ai * ( rI[0] - rJ[0] ) ;
                    c1y  = ai * ( rI[1] - rJ[1] ) ;
                    c1z  = ai * ( rI[2] - rJ[2] ) ;
                    c3x  = aj * ( rJ[0] - rI[0] ) ;
                    c3y  = aj * ( rJ[1] - rI[1] ) ;
                    c3z  = aj * ( rJ[2] - rI[2] ) ;
                    /* . Calculate the rys polynomial roots. */
                    RysQuadrature_Roots ( &roots, nroots, ( rho * rIJ2 ) ) ;
                    /* . Loop over the roots and construct the subsidiary integrals. */
                    for ( m = 0 ; m < nroots ; m++ )
                    {
                        u2    = roots.roots[m] * rho ;
                        f00   = roots.weights[m] ;
                        fac   = 1.0e+00 / ( ab + u2 * aandb ) ;
                        fac2  = 0.5e+00 * fac ;
                        bp01  = ( ai + u2 ) * fac2 ;
                        b00   =        u2   * fac2 ;
                        b10   = ( aj + u2 ) * fac2 ;
                        xcp00 = u2 * c1x * fac ;
                        xc00  = u2 * c3x * fac ;
                        ycp00 = u2 * c1y * fac ;
                        yc00  = u2 * c3y * fac ;
                        zcp00 = u2 * c1z * fac ;
                        zc00  = u2 * c3z * fac ;
                        Subsidiary_Integral_Nuclear2C ( iammax         ,
                                                        jammax         ,
                                                        b00            ,
                                                        b10            ,
                                                        bp01           ,
                                                        f00            ,
                                                        xc00           ,
                                                        xcp00          ,
                                                        yc00           ,
                                                        ycp00          ,
                                                        zc00           ,
                                                        zcp00          ,
                                                        jdim           ,
                                                        &xint[m*jdimm] ,
                                                        &yint[m*jdimm] ,
                                                        &zint[m*jdimm] ) ;
                    }
                    /* . Assemble the integrals. */
                    for ( i = 0, n = 0 ; i < ncfunci ; i++ )
                    {
   	                ix = CBFPOWX[i+icbfind] * jdim ;
	                iy = CBFPOWY[i+icbfind] * jdim ;
	                iz = CBFPOWZ[i+icbfind] * jdim ;
                        ti = dfIJ * iBasis->shells[iShell].primitives[ip].ccbf[i] ;
                        for ( j = 0 ; j < ncfuncj ; j++ )
                        {
	                    jxix = CBFPOWX[j+jcbfind] + ix ;
	                    jyiy = CBFPOWY[j+jcbfind] + iy ;
	                    jziz = CBFPOWZ[j+jcbfind] + iz ;
                            for ( m = 0, fac = 0.0e+00 ; m < nroots ; m++ ) fac += xint[jxix+m*jdimm] * yint[jyiy+m*jdimm] * zint[jziz+m*jdimm] ;
                            tIJ = ti * jBasis->shells[jShell].primitives[jp].ccbf[j] ;
                            g[n] += tIJ * fac ;
                            n++ ;
                        }
                    }
                } /* . jP. */
            } /* . iP. */
            /* . Put the integrals in the proper place. */
            for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nbasisw ; i++ )
            {
                for ( j = 0 ; j < jBasis->shells[jShell].nbasisw ; j++, n++ )
                {
                    Array2D_Item ( integrals, i+iBasis->shells[iShell].nstartw, j+jBasis->shells[jShell].nstartw ) = g[n] ;
                }
            }
        } /* . jShell. */
    } /* . iShell. */
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Coulomb derivatives.
! . sX, sY and sZ are overwritten by this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisIntegrals_2CoulombD ( const GaussianBasis *iBasis ,
                                        const Real          *rI     ,
                                        const GaussianBasis *jBasis ,
                                        const Real          *rJ     ,
                                              RealArray2D   *sX     ,
                                              RealArray2D   *sY     ,
                                              RealArray2D   *sZ     )
{
    Integer       i, iammax, icbfind, ip, iShell, ix, iy, iz,
                  j, jammax, jcbfind, jdim, jdimd,  jdimm, jp, jShell,
                  jxix, jyiy, jziz, m, n, ncfunci, ncfuncj, nroots ;
    Real          aandb, ab, ai, aj, bp01, b00, b10, c1x, c1y, c1z, c3x, c3y, c3z, dfi, dfIJ,
                  fac, facx, facy, facz, fac2, f00, rho, rIJ2, ti, tIJ, u2,
                  xc00, xcp00, xIJ, yc00, ycp00, yIJ, zc00, zcp00, zIJ ;
    Real          gx[MAXCBF*MAXCBF], gy[MAXCBF*MAXCBF], gz[MAXCBF*MAXCBF],
                  xind[MAXAMP1*MAXAMP1*_MAXRYS], yind[MAXAMP1*MAXAMP1*_MAXRYS], zind[MAXAMP1*MAXAMP1*_MAXRYS],
                  xint[MAXAMP1*MAXAMP2*_MAXRYS], yint[MAXAMP1*MAXAMP2*_MAXRYS], zint[MAXAMP1*MAXAMP2*_MAXRYS] ;
    RysQuadrature roots ;
    /* . Initialization. */
    xIJ  = rI[0] - rJ[0] ;
    yIJ  = rI[1] - rJ[1] ;
    zIJ  = rI[2] - rJ[2] ;
    rIJ2 = xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ;
    RealArray2D_Set ( sX, 0.0e+00 ) ;
    RealArray2D_Set ( sY, 0.0e+00 ) ;
    RealArray2D_Set ( sZ, 0.0e+00 ) ;
    /* . Outer loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nshells ; iShell++ )
    {
        iammax  = iBasis->shells[iShell].type->angularmomentum_high ;
        icbfind = iBasis->shells[iShell].type->cbfindex ;
        ncfunci = iBasis->shells[iShell].type->ncbf     ;
        /* . Inner loop over shells. */
        for ( jShell = 0 ; jShell < jBasis->nshells ; jShell++ )
        {

            jammax  = jBasis->shells[jShell].type->angularmomentum_high ;
            jdim    = jammax + 1 ;
            jdimd   = ( iammax + 1 ) * ( jammax + 1 ) ;
            jdimm   = ( iammax + 2 ) * ( jammax + 1 ) ;
            jcbfind = jBasis->shells[jShell].type->cbfindex ;
            ncfuncj = jBasis->shells[jShell].type->ncbf     ;
            /* . Get the number of roots. */
            nroots = ( iammax + jammax + 1 ) / 2 + 1 ;
            /* . Initialize the integral blocks. */
            for ( i = 0 ; i < ( ncfunci * ncfuncj ) ; i++ ) gx[i] = gy[i] = gz[i] = 0.0e+00 ;
            /* . Outer loop over primitives. */
            for ( ip = 0 ; ip < iBasis->shells[iShell].nprimitives ; ip++ )
            {
                /* . Get some information for the primitive. */
	        ai  = iBasis->shells[iShell].primitives[ip].exponent ;
                dfi = PI252 / ai ;
                /* . Inner loop over primitives. */
                for ( jp = 0 ; jp < jBasis->shells[jShell].nprimitives ; jp++ )
                {
                    /* . Get some information for the primitive. */
	            aj = jBasis->shells[jShell].primitives[jp].exponent ;
                    /* . Calculate some factors. */
                    ab    = ai * aj ;
                    aandb = ai + aj ;
                    rho   = ab / aandb ;
                    dfIJ  = dfi / ( aj * sqrt ( aandb ) ) ;
                    /* . Calculate some displacements. */
                    c1x  = ai * ( rI[0] - rJ[0] ) ;
                    c1y  = ai * ( rI[1] - rJ[1] ) ;
                    c1z  = ai * ( rI[2] - rJ[2] ) ;
                    c3x  = aj * ( rJ[0] - rI[0] ) ;
                    c3y  = aj * ( rJ[1] - rI[1] ) ;
                    c3z  = aj * ( rJ[2] - rI[2] ) ;
                    /* . Calculate the rys polynomial roots. */
                    RysQuadrature_Roots ( &roots, nroots, ( rho * rIJ2 ) ) ;
                    /* . Loop over the roots and construct the subsidiary integrals. */
                    for ( m = 0 ; m < nroots ; m++ )
                    {
                        u2    = roots.roots[m] * rho ;
                        f00   = roots.weights[m] ;
                        fac   = 1.0e+00 / ( ab + u2 * aandb ) ;
                        fac2  = 0.5e+00 * fac ;
                        bp01  = ( ai + u2 ) * fac2 ;
                        b00   =        u2   * fac2 ;
                        b10   = ( aj + u2 ) * fac2 ;
                        xcp00 = u2 * c1x * fac ;
                        xc00  = u2 * c3x * fac ;
                        ycp00 = u2 * c1y * fac ;
                        yc00  = u2 * c3y * fac ;
                        zcp00 = u2 * c1z * fac ;
                        zc00  = u2 * c3z * fac ;
                        Subsidiary_Integral_Nuclear2C ( iammax+1, jammax, b00, b10, bp01, f00, xc00, xcp00, yc00, ycp00, zc00, zcp00,
                                                        jdim, &xint[m*jdimm], &yint[m*jdimm], &zint[m*jdimm] ) ;
                        Subsidiary_Integral_Derivative2 ( &xint[m*jdimm], &yint[m*jdimm], &zint[m*jdimm], ai, iammax, jammax,
                                                          jdim, &xind[m*jdimd], &yind[m*jdimd], &zind[m*jdimd] ) ;
                    }
                    /* . Assemble the integrals. */
                    for ( i = 0, n = 0 ; i < ncfunci ; i++ )
                    {
   	                ix = CBFPOWX[i+icbfind] * jdim ;
	                iy = CBFPOWY[i+icbfind] * jdim ;
	                iz = CBFPOWZ[i+icbfind] * jdim ;
                        ti = dfIJ * iBasis->shells[iShell].primitives[ip].ccbf[i] ;
                        for ( j = 0 ; j < ncfuncj ; j++ )
                        {
	                    jxix = CBFPOWX[j+jcbfind] + ix ;
	                    jyiy = CBFPOWY[j+jcbfind] + iy ;
	                    jziz = CBFPOWZ[j+jcbfind] + iz ;
                            for ( m = 0, facx = facy = facz = 0.0e+00 ; m < nroots ; m++ )
                            {
                                facx += xind[jxix+m*jdimd] * yint[jyiy+m*jdimm] * zint[jziz+m*jdimm] ;
                                facy += xint[jxix+m*jdimm] * yind[jyiy+m*jdimd] * zint[jziz+m*jdimm] ;
                                facz += xint[jxix+m*jdimm] * yint[jyiy+m*jdimm] * zind[jziz+m*jdimd] ;
                            }
                            tIJ = ti * jBasis->shells[jShell].primitives[jp].ccbf[j] ;
                            gx[n] += tIJ * facx ;
                            gy[n] += tIJ * facy ;
                            gz[n] += tIJ * facz ;
                            n++ ;
                        }
                    }
                } /* . jP. */
            } /* . iP. */
            /* . Put the integrals in the proper place. */
            for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nbasisw ; i++ )
            {
                for ( j = 0 ; j < jBasis->shells[jShell].nbasisw ; j++, n++ )
                {
                    Array2D_Item ( sX, i+iBasis->shells[iShell].nstartw, j+jBasis->shells[jShell].nstartw ) = gx[n] ;
                    Array2D_Item ( sY, i+iBasis->shells[iShell].nstartw, j+jBasis->shells[jShell].nstartw ) = gy[n] ;
                    Array2D_Item ( sZ, i+iBasis->shells[iShell].nstartw, j+jBasis->shells[jShell].nstartw ) = gz[n] ;
                }
            }
        } /* . jShell. */
    } /* . iShell. */
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Overlap integrals.
! . integrals is overwritten by this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisIntegrals_2Overlap ( const GaussianBasis *iBasis    ,
                                       const Real          *rI        ,
                                       const GaussianBasis *jBasis    ,
                                       const Real          *rJ        ,
                                             RealArray2D   *integrals )
{
    Boolean  iIsJ ;
    Integer  i, iammax, icbfind, ip, iShell, ix, iy, iz,
             j, jammax, jcbfind, jdim, jp, jShell, jUpper, jxix, jyiy, jziz,
             n, ncfunci, ncfuncj ;
    Real     aa, aainv, ai, aj, arri, expfac, fac, rIJ2, xIJ, yIJ, zIJ ;
    Real     ar[3], arI[3], s[MAXCBF*MAXCBF],
             xo[MAXAMP1*MAXAMP1], yo[MAXAMP1*MAXAMP1], zo[MAXAMP1*MAXAMP1] ;
    /* . Initialization. */
    iIsJ = ( iBasis == jBasis ) && ( rI == rJ ) ;
    xIJ  = rI[0] - rJ[0] ;
    yIJ  = rI[1] - rJ[1] ;
    zIJ  = rI[2] - rJ[2] ;
    rIJ2 = xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ;
    RealArray2D_Set ( integrals, 0.0e+00 ) ;
    /* . Outer loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nshells ; iShell++ )
    {
        iammax  = iBasis->shells[iShell].type->angularmomentum_high ;
        icbfind = iBasis->shells[iShell].type->cbfindex ;
        ncfunci = iBasis->shells[iShell].type->ncbf     ;
        /* . Inner loop over shells. */
        if ( iIsJ ) jUpper = iShell + 1 ;
        else        jUpper = jBasis->nshells ;
        for ( jShell = 0 ; jShell < jUpper ; jShell++ )
        {
            jammax  = jBasis->shells[jShell].type->angularmomentum_high ;
            jdim    = jammax + 1 ;
            jcbfind = jBasis->shells[jShell].type->cbfindex ;
            ncfuncj = jBasis->shells[jShell].type->ncbf     ;
            /* . Initialize the integral block. */
            for ( i = 0 ; i < ( ncfunci * ncfuncj ) ; i++ ) s[i] = 0.0e+00 ;
            /* . Outer loop over primitives. */
            for ( ip = 0 ; ip < iBasis->shells[iShell].nprimitives ; ip++ )
            {
                /* . Get some information for the primitive. */
	        ai   = iBasis->shells[iShell].primitives[ip].exponent ;
	        arri = ai * rIJ2 ;
	        for ( i = 0 ; i < 3 ; i++ ) arI[i] = ai * rI[i] ;
                /* . Inner loop over primitives. */
                for ( jp = 0 ; jp < jBasis->shells[jShell].nprimitives ; jp++ )
                {
                    /* . Get some information for the primitive. */
	            aj    = jBasis->shells[jShell].primitives[jp].exponent ;
	            aa    = ai + aj ;
	            aainv = 1.0e+00 / aa ;
	            fac   = aj * arri * aainv ;
	            if ( fac > PRIMITIVE_OVERLAP_TOLERANCE ) continue ;
                    expfac = exp ( - fac ) ;
	            for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( arI[i] + aj * rJ[i] ) * aainv ;
                    /* . Calculate the overlap integrals. */
                    Subsidiary_Integral_Overlap2 ( xo, yo, zo, aa, ar, rI, rJ, iammax, jammax ) ;
                    /* . Add in the contributions to the full integrals. */
                    for ( i = n = 0 ; i < ncfunci ; i++ )
                    {
   	                ix = CBFPOWX[i+icbfind] * jdim ;
	                iy = CBFPOWY[i+icbfind] * jdim ;
	                iz = CBFPOWZ[i+icbfind] * jdim ;
                        for ( j = 0 ; j < ncfuncj ; j++ )
                        {
	                    jxix = CBFPOWX[j+jcbfind] + ix ;
	                    jyiy = CBFPOWY[j+jcbfind] + iy ;
	                    jziz = CBFPOWZ[j+jcbfind] + iz ;
    		            s[n] += expfac * iBasis->shells[iShell].primitives[ip].ccbf[i] *
                                             jBasis->shells[jShell].primitives[jp].ccbf[j] * xo[jxix] * yo[jyiy] * zo[jziz] ;
                            n++ ;
                        }
                    }
                } /* . jP. */
            } /* . iP. */
            /* . Put the integrals in the proper place. */
            for ( i = n = 0 ; i < iBasis->shells[iShell].nbasisw ; i++ )
            {
                for ( j = 0 ; j < jBasis->shells[jShell].nbasisw ; j++, n++ )
                {
                    Array2D_Item ( integrals, i+iBasis->shells[iShell].nstartw, j+jBasis->shells[jShell].nstartw ) = s[n] ;
                }
            }
        } /* . jShell. */
    } /* . iShell. */
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Overlap derivatives.
! . overlapX, overlapY and overlapZ are overwritten by this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisIntegrals_2OverlapD ( const GaussianBasis *iBasis   ,
                                        const Real          *rI       ,
                                        const GaussianBasis *jBasis   ,
                                        const Real          *rJ       ,
                                              RealArray2D   *overlapX ,
                                              RealArray2D   *overlapY ,
                                              RealArray2D   *overlapZ )
{
    Integer  i, iammax, icbfind, ip, iShell, ix, iy, iz,
             j, jammax, jcbfind, jdim, jp, jShell, jxix, jyiy, jziz, n, ncfunci, ncfuncj ;
    Real     aa, aainv, ai, aj, arri, denfac, expfac, fac, rIJ2, xIJ, yIJ, zIJ ;
    Real     ar[3], arI[3], sx[MAXCBF*MAXCBF],       sy[MAXCBF*MAXCBF],       sz[MAXCBF*MAXCBF],
                            xd[MAXAMP1*MAXAMP1],     yd[MAXAMP1*MAXAMP1],     zd[MAXAMP1*MAXAMP1],
                            xo[MAXAMP1*(MAXAMP1+1)], yo[MAXAMP1*(MAXAMP1+1)], zo[MAXAMP1*(MAXAMP1+1)] ;
    /* . Initialization. */
    xIJ  = rI[0] - rJ[0] ;
    yIJ  = rI[1] - rJ[1] ;
    zIJ  = rI[2] - rJ[2] ;
    rIJ2 = xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ;
    RealArray2D_Set ( overlapX, 0.0e+00 ) ;
    RealArray2D_Set ( overlapY, 0.0e+00 ) ;
    RealArray2D_Set ( overlapZ, 0.0e+00 ) ;
    /* . Outer loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nshells ; iShell++ )
    {
        iammax  = iBasis->shells[iShell].type->angularmomentum_high ;
        icbfind = iBasis->shells[iShell].type->cbfindex ;
        ncfunci = iBasis->shells[iShell].type->ncbf     ;
        /* . Inner loop over shells. */
        for ( jShell = 0 ; jShell < jBasis->nshells ; jShell++ )
        {
            jammax  = jBasis->shells[jShell].type->angularmomentum_high ;
            jdim    = jammax + 1 ;
            jcbfind = jBasis->shells[jShell].type->cbfindex ;
            ncfuncj = jBasis->shells[jShell].type->ncbf     ;
            /* . Initialize the integral block. */
            for ( i = 0 ; i < ( ncfunci * ncfuncj ) ; i++ ) { sx[i] = 0.0e+00 ; sy[i] = 0.0e+00 ; sz[i] = 0.0e+00 ; }
            /* . Outer loop over primitives. */
            for ( ip = 0 ; ip < iBasis->shells[iShell].nprimitives ; ip++ )
            {
                /* . Get some information for the primitive. */
	        ai   = iBasis->shells[iShell].primitives[ip].exponent ;
	        arri = ai * rIJ2 ;
	        for ( i = 0 ; i < 3 ; i++ ) arI[i] = ai * rI[i] ;
                /* . Inner loop over primitives. */
                for ( jp = 0 ; jp < jBasis->shells[jShell].nprimitives ; jp++ )
                {
                    /* . Get some information for the primitive. */
	            aj    = jBasis->shells[jShell].primitives[jp].exponent ;
	            aa    = ai + aj ;
	            aainv = 1.0e+00 / aa ;
	            fac   = aj * arri * aainv ;
	            if ( fac > PRIMITIVE_OVERLAP_TOLERANCE ) continue ;
                    expfac = exp ( - fac ) ;
	            for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( arI[i] + aj * rJ[i] ) * aainv ;
                    /* . Calculate the overlap integrals and derivatives. */
                    Subsidiary_Integral_Overlap2    ( xo, yo, zo, aa, ar, rI, rJ, iammax+1, jammax ) ;
                    Subsidiary_Integral_Derivative2 ( xo, yo, zo, ai, iammax, jammax, jdim, xd, yd, zd ) ;
                    /* . Add in the contributions to the full integrals. */
                    for ( i = 0, n = 0 ; i < ncfunci ; i++ )
                    {
   	                ix = CBFPOWX[i+icbfind] * jdim ;
	                iy = CBFPOWY[i+icbfind] * jdim ;
	                iz = CBFPOWZ[i+icbfind] * jdim ;
                        for ( j = 0 ; j < ncfuncj ; j++, n++ )
                        {
	                    jxix = CBFPOWX[j+jcbfind] + ix ;
	                    jyiy = CBFPOWY[j+jcbfind] + iy ;
	                    jziz = CBFPOWZ[j+jcbfind] + iz ;
                            denfac = expfac * iBasis->shells[iShell].primitives[ip].ccbf[i] * jBasis->shells[jShell].primitives[jp].ccbf[j] ;
    		            sx[n] += denfac * xd[jxix] * yo[jyiy] * zo[jziz] ;
    		            sy[n] += denfac * xo[jxix] * yd[jyiy] * zo[jziz] ;
    		            sz[n] += denfac * xo[jxix] * yo[jyiy] * zd[jziz] ;
                        }
                    }
                } /* . jP. */
            } /* . iP. */
            /* . Put the integrals in the proper place. */
            for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nbasisw ; i++ )
            {
                for ( j = 0 ; j < jBasis->shells[jShell].nbasisw ; j++, n++ )
                {
                    Array2D_Item ( overlapX, i+iBasis->shells[iShell].nstartw, j+jBasis->shells[jShell].nstartw ) = sx[n] ;
                    Array2D_Item ( overlapY, i+iBasis->shells[iShell].nstartw, j+jBasis->shells[jShell].nstartw ) = sy[n] ;
                    Array2D_Item ( overlapZ, i+iBasis->shells[iShell].nstartw, j+jBasis->shells[jShell].nstartw ) = sz[n] ;
                }
            }
        } /* . jShell. */
    } /* . iShell. */
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Dipole integrals.
! . dipoleX, dipoleY and dipoleZ are overwritten by this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisIntegrals_Dipole ( const GaussianBasis *iBasis  ,
                                     const Real          *rI      ,
                                     const GaussianBasis *jBasis  ,
                                     const Real          *rJ      ,
                                     const Real          *center  ,
                                           RealArray2D   *dipoleX ,
                                           RealArray2D   *dipoleY ,
                                           RealArray2D   *dipoleZ )
{
    Boolean  iIsJ ;
    Integer  i, iammax, icbfind,  ip, iShell,  ix, iy, iz,
             j, jammax, jcbfind, jdim, jp,  jShell,  jUpper, jxix, jyiy, jziz, n, ncfunci, ncfuncj ;
    Real     aa, aainv, ai, aj, arri, expfac, fac, rIJ2, ti, tIJ, xIJ, yIJ, zIJ ;
    Real     ar[3], arI[3], sx[MAXCBF*MAXCBF], sy[MAXCBF*MAXCBF], sz[MAXCBF*MAXCBF],
             xo[MAXAMP1*MAXAMP1], yo[MAXAMP1*MAXAMP1], zo[MAXAMP1*MAXAMP1] ,
             xd[MAXAMP1*MAXAMP1], yd[MAXAMP1*MAXAMP1], zd[MAXAMP1*MAXAMP1] ;
    /* . Initialization. */
    iIsJ = ( iBasis == jBasis ) && ( rI == rJ ) ;
    xIJ  = rI[0] - rJ[0] ;
    yIJ  = rI[1] - rJ[1] ;
    zIJ  = rI[2] - rJ[2] ;
    rIJ2 = xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ;
    RealArray2D_Set ( dipoleX, 0.0e+00 ) ;
    RealArray2D_Set ( dipoleY, 0.0e+00 ) ;
    RealArray2D_Set ( dipoleZ, 0.0e+00 ) ;
    /* . Outer loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nshells ; iShell++ )
    {
        iammax  = iBasis->shells[iShell].type->angularmomentum_high ;
        icbfind = iBasis->shells[iShell].type->cbfindex ;
        ncfunci = iBasis->shells[iShell].type->ncbf     ;
        /* . Inner loop over shells. */
        if ( iIsJ ) jUpper = iShell + 1 ;
        else        jUpper = jBasis->nshells ;
        for ( jShell = 0 ; jShell < jUpper ; jShell++ )
        {
            jammax  = jBasis->shells[jShell].type->angularmomentum_high ;
            jdim    = jammax + 1 ;
            jcbfind = jBasis->shells[jShell].type->cbfindex ;
            ncfuncj = jBasis->shells[jShell].type->ncbf     ;
            /* . Initialize the integral blocks. */
            for ( i = 0 ; i < ( ncfunci * ncfuncj ) ; i++ )
            {
                sx[i] = 0.0e+00 ;
                sy[i] = 0.0e+00 ;
                sz[i] = 0.0e+00 ;
            }
            /* . Outer loop over primitives. */
            for ( ip = 0 ; ip < iBasis->shells[iShell].nprimitives ; ip++ )
            {
                /* . Get some information for the primitive. */
	        ai   = iBasis->shells[iShell].primitives[ip].exponent ;
	        arri = ai * rIJ2 ;
	        for ( i = 0 ; i < 3 ; i++ ) arI[i] = ai * rI[i] ;
                /* . Inner loop over primitives. */
                for ( jp = 0 ; jp < jBasis->shells[jShell].nprimitives ; jp++ )
                {
                    /* . Get some information for the primitive. */
	            aj    = jBasis->shells[jShell].primitives[jp].exponent ;
	            aa    = ai + aj ;
	            aainv = 1.0e+00 / aa ;
	            fac   = aj * arri * aainv ;
	            if ( fac > PRIMITIVE_OVERLAP_TOLERANCE ) continue ;
                    expfac = exp ( - fac ) ;
	            for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( arI[i] + aj * rJ[i] ) * aainv ;
                    /* . Calculate the subsidiary integrals. */
                    Subsidiary_Integral_Overlap2 ( xo, yo, zo, aa, ar, rI, rJ,         iammax, jammax ) ;
                    Subsidiary_Integral_Dipole   ( xd, yd, zd, aa, ar, rI, rJ, center, iammax, jammax ) ;
                    /* . Add in the contributions to the full integrals. */
                    for ( i = 0, n = 0 ; i < ncfunci ; i++ )
                    {
   	                ix = CBFPOWX[i+icbfind] * jdim ;
	                iy = CBFPOWY[i+icbfind] * jdim ;
	                iz = CBFPOWZ[i+icbfind] * jdim ;
                        ti = expfac * iBasis->shells[iShell].primitives[ip].ccbf[i] ;
                        for ( j = 0 ; j < ncfuncj ; j++, n++ )
                        {
	                    jxix = CBFPOWX[j+jcbfind] + ix ;
	                    jyiy = CBFPOWY[j+jcbfind] + iy ;
	                    jziz = CBFPOWZ[j+jcbfind] + iz ;
                            tIJ  = ti * jBasis->shells[jShell].primitives[jp].ccbf[j] ;
    		            sx[n] += tIJ * xd[jxix] * yo[jyiy] * zo[jziz] ;
    		            sy[n] += tIJ * xo[jxix] * yd[jyiy] * zo[jziz] ;
    		            sz[n] += tIJ * xo[jxix] * yo[jyiy] * zd[jziz] ;
                        }
                    }
                } /* . jP. */
            } /* . iP. */
            /* . Put the integrals in the proper place. */
            for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nbasisw ; i++ )
            {
                for ( j = 0 ; j < jBasis->shells[jShell].nbasisw ; j++, n++ )
                {
                    Array2D_Item ( dipoleX, i+iBasis->shells[iShell].nstartw, j+jBasis->shells[jShell].nstartw ) = sx[n] ;
                    Array2D_Item ( dipoleY, i+iBasis->shells[iShell].nstartw, j+jBasis->shells[jShell].nstartw ) = sy[n] ;
                    Array2D_Item ( dipoleZ, i+iBasis->shells[iShell].nstartw, j+jBasis->shells[jShell].nstartw ) = sz[n] ;
                }
            }
        } /* . jShell. */
    } /* . iShell. */
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Kinetic energy and overlap integrals.
! . Kinetic and overlap are overwritten by this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisIntegrals_Kinetic2Overlap ( const GaussianBasis *iBasis  ,
                                              const Real          *rI      ,
                                              const GaussianBasis *jBasis  ,
                                              const Real          *rJ      ,
                                                    RealArray2D   *overlap ,
                                                    RealArray2D   *kinetic )
{
    Boolean  iIsJ ;
    Integer  i, iammax, icbfind,  ip, iShell,  ixo, ixt, iyo, iyt, izo, izt,
             j, jammax, jcbfind, jdimo, jdimt, jp,  jShell,  jUpper,
             jxixo, jyiyo, jzizo, jxixt, jyiyt, jzizt, n, ncfunci, ncfuncj ;
    Real     aa, aainv, ai, aj, arri, expfac, fac, rIJ2, ti, tIJ, xIJ, yIJ, zIJ ;
    Real     ar[3], arI[3], s[MAXCBF*MAXCBF], t[MAXCBF*MAXCBF],
             xo[MAXAMP1*MAXAMP3], yo[MAXAMP1*MAXAMP3], zo[MAXAMP1*MAXAMP3] ,
             xt[MAXAMP1*MAXAMP1], yt[MAXAMP1*MAXAMP1], zt[MAXAMP1*MAXAMP1] ;
    /* . Initialization. */
    iIsJ = ( iBasis == jBasis ) && ( rI == rJ ) ;
    xIJ  = rI[0] - rJ[0] ;
    yIJ  = rI[1] - rJ[1] ;
    zIJ  = rI[2] - rJ[2] ;
    rIJ2 = xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ;
    RealArray2D_Set ( kinetic, 0.0e+00 ) ;
    RealArray2D_Set ( overlap, 0.0e+00 ) ;
    /* . Outer loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nshells ; iShell++ )
    {
        iammax  = iBasis->shells[iShell].type->angularmomentum_high ;
        icbfind = iBasis->shells[iShell].type->cbfindex ;
        ncfunci = iBasis->shells[iShell].type->ncbf     ;
        /* . Inner loop over shells. */
        if ( iIsJ ) jUpper = iShell + 1 ;
        else        jUpper = jBasis->nshells ;
        for ( jShell = 0 ; jShell < jUpper ; jShell++ )
        {
            jammax  = jBasis->shells[jShell].type->angularmomentum_high ;
            jdimo   = jammax + 3 ;
            jdimt   = jammax + 1 ;
            jcbfind = jBasis->shells[jShell].type->cbfindex ;
            ncfuncj = jBasis->shells[jShell].type->ncbf     ;
            /* . Initialize the integral blocks. */
            for ( i = 0 ; i < ( ncfunci * ncfuncj ) ; i++ )
            {
                s[i] = 0.0e+00 ;
                t[i] = 0.0e+00 ;
            }
            /* . Outer loop over primitives. */
            for ( ip = 0 ; ip < iBasis->shells[iShell].nprimitives ; ip++ )
            {
                /* . Get some information for the primitive. */
	        ai   = iBasis->shells[iShell].primitives[ip].exponent ;
	        arri = ai * rIJ2 ;
	        for ( i = 0 ; i < 3 ; i++ ) arI[i] = ai * rI[i] ;

                /* . Inner loop over primitives. */
                for ( jp = 0 ; jp < jBasis->shells[jShell].nprimitives ; jp++ )
                {
                    /* . Get some information for the primitive. */
	            aj    = jBasis->shells[jShell].primitives[jp].exponent ;
	            aa    = ai + aj ;
	            aainv = 1.0e+00 / aa ;
	            fac   = aj * arri * aainv ;
	            if ( fac > PRIMITIVE_OVERLAP_TOLERANCE ) continue ;
                    expfac = exp ( - fac ) ;
	            for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( arI[i] + aj * rJ[i] ) * aainv ;

                    /* . Calculate the subsidiary integrals. */
                    Subsidiary_Integral_Overlap2 ( xo, yo, zo, aa, ar, rI, rJ, iammax, jammax + 2 ) ;
                    Subsidiary_Integral_Kinetic  ( xo, yo, zo, xt, yt, zt, aj, iammax, jammax, jdimo, jdimt ) ;

                    /* . Add in the contributions to the full integrals. */
                    for ( i = 0, n = 0 ; i < ncfunci ; i++ )
                    {
   	                ixo = CBFPOWX[i+icbfind] * jdimo ;
	                iyo = CBFPOWY[i+icbfind] * jdimo ;
	                izo = CBFPOWZ[i+icbfind] * jdimo ;
   	                ixt = CBFPOWX[i+icbfind] * jdimt ;
	                iyt = CBFPOWY[i+icbfind] * jdimt ;
	                izt = CBFPOWZ[i+icbfind] * jdimt ;
                        ti = expfac * iBasis->shells[iShell].primitives[ip].ccbf[i] ;
                        for ( j = 0 ; j < ncfuncj ; j++, n++ )
                        {
	                    jxixo = CBFPOWX[j+jcbfind] + ixo ;
	                    jyiyo = CBFPOWY[j+jcbfind] + iyo ;
	                    jzizo = CBFPOWZ[j+jcbfind] + izo ;
	                    jxixt = CBFPOWX[j+jcbfind] + ixt ;
	                    jyiyt = CBFPOWY[j+jcbfind] + iyt ;
	                    jzizt = CBFPOWZ[j+jcbfind] + izt ;
                            tIJ  = ti * jBasis->shells[jShell].primitives[jp].ccbf[j] ;
    		            s[n] += tIJ * xo[jxixo] * yo[jyiyo] * zo[jzizo] ;
                            t[n] += tIJ * ( xt[jxixt] * yo[jyiyo] * zo[jzizo] +
                                            xo[jxixo] * yt[jyiyt] * zo[jzizo] +
                                            xo[jxixo] * yo[jyiyo] * zt[jzizt] ) ;
                        }
                    }
                } /* . jP. */
            } /* . iP. */
            /* . Put the integrals in the proper place. */
            for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nbasisw ; i++ )
            {
                for ( j = 0 ; j < jBasis->shells[jShell].nbasisw ; j++, n++ )
                {
                    Array2D_Item ( overlap, i+iBasis->shells[iShell].nstartw, j+jBasis->shells[jShell].nstartw ) = s[n] ;
                    Array2D_Item ( kinetic, i+iBasis->shells[iShell].nstartw, j+jBasis->shells[jShell].nstartw ) = t[n] ;
                }
            }
        } /* . jShell. */
    } /* . iShell. */
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Kinetic energy and overlap derivatives.
! . kineticX, kineticY, kineticZ, overlapX, overlapY and overlapZ are overwritten by this function.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisIntegrals_Kinetic2OverlapD ( const GaussianBasis *iBasis   ,
                                               const Real          *rI       ,
                                               const GaussianBasis *jBasis   ,
                                               const Real          *rJ       ,
                                                     RealArray2D   *overlapX ,
                                                     RealArray2D   *overlapY ,
                                                     RealArray2D   *overlapZ ,
                                                     RealArray2D   *kineticX ,
                                                     RealArray2D   *kineticY ,
                                                     RealArray2D   *kineticZ )
{
     Integer  i, iammax, icbfind, ip, iShell,  ixo, ixt, iyo, iyt, izo, izt,
              j, jammax, jcbfind, jdimo, jdimt, jp,  jShell,  jxixo, jyiyo, jzizo,
              jxixt, jyiyt, jzizt, n, ncfunci, ncfuncj ;
     Real     aa, aainv, ai, aj, arri, expfac, fac, rIJ2, ti, tIJ, xIJ, yIJ, zIJ ;
     Real     ar[3], arI[3], sx[MAXCBF*MAXCBF], sy[MAXCBF*MAXCBF], sz[MAXCBF*MAXCBF],
                             tx[MAXCBF*MAXCBF], ty[MAXCBF*MAXCBF], tz[MAXCBF*MAXCBF],
                             xo[MAXAMP2*MAXAMP3], yo[MAXAMP2*MAXAMP3], zo[MAXAMP2*MAXAMP3],
                             xt[MAXAMP1*MAXAMP2], yt[MAXAMP1*MAXAMP2], zt[MAXAMP1*MAXAMP2],
                             xod[MAXAMP1*MAXAMP3], yod[MAXAMP1*MAXAMP3], zod[MAXAMP1*MAXAMP3],
                             xtd[MAXAMP1*MAXAMP1], ytd[MAXAMP1*MAXAMP1], ztd[MAXAMP1*MAXAMP1] ;
    /* . Initialization. */
    xIJ  = rI[0] - rJ[0] ;
    yIJ  = rI[1] - rJ[1] ;
    zIJ  = rI[2] - rJ[2] ;
    rIJ2 = xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ;
    RealArray2D_Set ( kineticX, 0.0e+00 ) ;
    RealArray2D_Set ( kineticY, 0.0e+00 ) ;
    RealArray2D_Set ( kineticZ, 0.0e+00 ) ;
    RealArray2D_Set ( overlapX, 0.0e+00 ) ;
    RealArray2D_Set ( overlapY, 0.0e+00 ) ;
    RealArray2D_Set ( overlapZ, 0.0e+00 ) ;
    /* . Outer loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nshells ; iShell++ )
    {
        iammax  = iBasis->shells[iShell].type->angularmomentum_high ;
        icbfind = iBasis->shells[iShell].type->cbfindex ;
        ncfunci = iBasis->shells[iShell].type->ncbf     ;
        /* . Inner loop over shells. */
        for ( jShell = 0 ; jShell < jBasis->nshells ; jShell++ )
        {
            jammax  = jBasis->shells[jShell].type->angularmomentum_high ;
            jdimo   = jammax + 3 ;
            jdimt   = jammax + 1 ;
            jcbfind = jBasis->shells[jShell].type->cbfindex ;
            ncfuncj = jBasis->shells[jShell].type->ncbf     ;
            /* . Initialize the integral blocks. */
            for ( i = 0 ; i < ( ncfunci * ncfuncj ) ; i++ )
            {
                sx[i] = sy[i] = sz[i] = 0.0e+00 ;
                tx[i] = ty[i] = tz[i] = 0.0e+00 ;
            }
            /* . Outer loop over primitives. */
            for ( ip = 0 ; ip < iBasis->shells[iShell].nprimitives ; ip++ )
            {
                /* . Get some information for the primitive. */
	        ai   = iBasis->shells[iShell].primitives[ip].exponent ;
	        arri = ai * rIJ2 ;
	        for ( i = 0 ; i < 3 ; i++ ) arI[i] = ai * rI[i] ;
                /* . Inner loop over primitives. */
                for ( jp = 0 ; jp < jBasis->shells[jShell].nprimitives ; jp++ )
                {
                    /* . Get some information for the primitive. */
	            aj    = jBasis->shells[jShell].primitives[jp].exponent ;
	            aa    = ai + aj ;
	            aainv = 1.0e+00 / aa ;
	            fac   = aj * arri * aainv ;
	            if ( fac > PRIMITIVE_OVERLAP_TOLERANCE ) continue ;
                    expfac = exp ( - fac ) ;
	            for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( arI[i] + aj * rJ[i] ) * aainv ;
                    /* . Calculate the subsidiary integrals. */
                    Subsidiary_Integral_Overlap2 ( xo, yo, zo, aa, ar, rI, rJ, iammax + 1, jammax + 2 ) ;
                    Subsidiary_Integral_Kinetic  ( xo, yo, zo, xt, yt, zt, aj, iammax + 1, jammax, jdimo, jdimt ) ;
                    Subsidiary_Integral_Derivative2 ( xo, yo, zo, ai, iammax, jammax, jdimo, xod, yod, zod ) ;
                    Subsidiary_Integral_Derivative2 ( xt, yt, zt, ai, iammax, jammax, jdimt, xtd, ytd, ztd ) ;
                    /* . Add in the contributions to the full integrals. */
                    for ( i = 0, n = 0 ; i < ncfunci ; i++ )
                    {
   	                ixo = CBFPOWX[i+icbfind] * jdimo ;
	                iyo = CBFPOWY[i+icbfind] * jdimo ;
	                izo = CBFPOWZ[i+icbfind] * jdimo ;
   	                ixt = CBFPOWX[i+icbfind] * jdimt ;
	                iyt = CBFPOWY[i+icbfind] * jdimt ;
	                izt = CBFPOWZ[i+icbfind] * jdimt ;
                        ti  = expfac * iBasis->shells[iShell].primitives[ip].ccbf[i] ;
                        for ( j = 0 ; j < ncfuncj ; j++, n++ )
                        {
	                    jxixo = CBFPOWX[j+jcbfind] + ixo ;
	                    jyiyo = CBFPOWY[j+jcbfind] + iyo ;
	                    jzizo = CBFPOWZ[j+jcbfind] + izo ;
	                    jxixt = CBFPOWX[j+jcbfind] + ixt ;
	                    jyiyt = CBFPOWY[j+jcbfind] + iyt ;
	                    jzizt = CBFPOWZ[j+jcbfind] + izt ;
                            tIJ   = ti * jBasis->shells[jShell].primitives[jp].ccbf[j] ;
    		            sx[n] += tIJ * xod[jxixo] * yo[jyiyo] * zo[jzizo] ;
    		            sy[n] += tIJ * xo[jxixo] * yod[jyiyo] * zo[jzizo] ;
    		            sz[n] += tIJ * xo[jxixo] * yo[jyiyo] * zod[jzizo] ;
                            tx[n] += tIJ * ( xtd[jxixt] * yo[jyiyo] * zo[jzizo] +
                                             xod[jxixo] * yt[jyiyt] * zo[jzizo] +
                                             xod[jxixo] * yo[jyiyo] * zt[jzizt] ) ;
                            ty[n] += tIJ * ( xt[jxixt] * yod[jyiyo] * zo[jzizo] +
                                             xo[jxixo] * ytd[jyiyt] * zo[jzizo] +
                                             xo[jxixo] * yod[jyiyo] * zt[jzizt] ) ;
                            tz[n] += tIJ * ( xt[jxixt] * yo[jyiyo] * zod[jzizo] +
                                             xo[jxixo] * yt[jyiyt] * zod[jzizo] +
                                             xo[jxixo] * yo[jyiyo] * ztd[jzizt] ) ;
                        }
                    }
                } /* . jP. */
            } /* . iP. */
            /* . Put the integrals in the proper place. */
            for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nbasisw ; i++ )
            {
                for ( j = 0 ; j < jBasis->shells[jShell].nbasisw ; j++, n++ )
                {
                    Array2D_Item ( overlapX, i+iBasis->shells[iShell].nstartw, j+jBasis->shells[jShell].nstartw ) = sx[n] ;
                    Array2D_Item ( overlapY, i+iBasis->shells[iShell].nstartw, j+jBasis->shells[jShell].nstartw ) = sy[n] ;
                    Array2D_Item ( overlapZ, i+iBasis->shells[iShell].nstartw, j+jBasis->shells[jShell].nstartw ) = sz[n] ;
                    Array2D_Item ( kineticX, i+iBasis->shells[iShell].nstartw, j+jBasis->shells[jShell].nstartw ) = tx[n] ;
                    Array2D_Item ( kineticY, i+iBasis->shells[iShell].nstartw, j+jBasis->shells[jShell].nstartw ) = ty[n] ;
                    Array2D_Item ( kineticZ, i+iBasis->shells[iShell].nstartw, j+jBasis->shells[jShell].nstartw ) = tz[n] ;
                }
            }
        } /* . jShell. */
    } /* . iShell. */
}
