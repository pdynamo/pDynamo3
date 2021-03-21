/*==================================================================================================================================
! . Integrals - 3 basis, 2 electrons, 0 nuclei/points.
!=================================================================================================================================*/

# include <math.h>
# include <stdlib.h>

# include "GaussianBasisIntegrals_b3e2n0.h"
# include "GaussianBasisSubsidiary.h"
# include "RysQuadrature.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Electron-fit integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisIntegrals_ElectronFit ( const GaussianBasis *iBasis ,
                                          const Real          *rI     ,
                                          const GaussianBasis *jBasis ,
                                          const Real          *rJ     ,
                                          const Real          *rIJ    ,
                                          const Real           rIJ2   ,
                                          const GaussianBasis *fBasis ,
                                          const Real          *rF     ,
                                                Block         *block  )
{
    Boolean       iIsJ, isDiagonal, QF0, QF1, QIJ0, QIJ1 ;
    Integer       dim1, dim2, dim3,
                  f, fammax,          fcbfind, fijx, fijy, fijz, fp, fShell,
                  i, iammax, iammaxt, icbfind, ip, iShell, ix, iy, iz,
                  j, jammax, jammaxt, jcbfind, jix, jiy, jiz, jp, jShell,
                  jUpper, m, n, ncfuncf, ncfunci, ncfuncj, nroots ;
    Real          aa, aandb, aainv, ab, ai, aj, arri, axac, ayac, azac, bp01, b00, b10,
                  c1x, c1y, c1z, c3x, c3y, c3z, c4x, c4y, c4z, dnuc, dxIJt, dyIJt, dzIJt,
                  expfac, expf, fac, fac2, f00, rho, ti, tij, tijf, u2,
                  xc00, xcp00, yc00, ycp00, zc00, zcp00 ;
    Real          ar[3], arI[3], g[MAXCBF*MAXCBF*MAXCBF],
                  xint[MAXAMP1*MAXAMP1*MAXAMP1*_MAXRYS], yint[MAXAMP1*MAXAMP1*MAXAMP1*_MAXRYS], zint[MAXAMP1*MAXAMP1*MAXAMP1*_MAXRYS] ;
    const Real   *rC ;
    RysQuadrature roots ;
    /* . Initialization. */
    block->count = 0 ;
    iIsJ         = ( iBasis == jBasis ) && ( rI == rJ ) ;
    /* . Triple loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nshells ; iShell++ )
    {
        /* . Get information about the shell. */
        iammax  = iBasis->shells[iShell].type->angularmomentum_high ;
        icbfind = iBasis->shells[iShell].type->cbfindex ;
        ncfunci = iBasis->shells[iShell].type->ncbf     ;
        /* . Set the Upper limit for the JSHELL loops. */
        if ( iIsJ ) jUpper = iShell + 1 ;
        else        jUpper = jBasis->nshells ;
        /* . Inner loop over shells. */
        for ( jShell = 0 ; jShell < jUpper ; jShell++ )
        {
            /* . Get information about the shell. */
            jammax  = jBasis->shells[jShell].type->angularmomentum_high ;
            jcbfind = jBasis->shells[jShell].type->cbfindex ;
            ncfuncj = jBasis->shells[jShell].type->ncbf     ;
            /* . Set the diagonal block flag. */
            isDiagonal = iIsJ && ( iShell == jShell ) ;
            /* . Set some flags. */
            QIJ0 = ( iammax + jammax == 0 ) ;
            QIJ1 = ( iammax + jammax <= 1 ) ;
            /* . Select the expansion center for the recurrence relations. */
            if ( iammax >= jammax )
            {
               iammaxt = iammax ;
               jammaxt = jammax ;
               dxIJt   = rIJ[0] ;
               dyIJt   = rIJ[1] ;
               dzIJt   = rIJ[2] ;
               rC      = rI ;
            }
            else
            {
               iammaxt = jammax ;
               jammaxt = iammax ;
               dxIJt   = - rIJ[0] ;
               dyIJt   = - rIJ[1] ;
               dzIJt   = - rIJ[2] ;
               rC      = rJ ;
            }
            /* . Loop over fitting shells. */
            for ( fShell = 0 ; fShell < fBasis->nshells ; fShell++ )
            {
                /* . Get information about the shell. */
                fammax  = fBasis->shells[fShell].type->angularmomentum_high ;
                fcbfind = fBasis->shells[fShell].type->cbfindex ;
                ncfuncf = fBasis->shells[fShell].type->ncbf     ;
                /* . Set some flags. */
                QF0 = ( fammax == 0 ) ;
                QF1 = ( fammax <= 1 ) ;
                /* . Get the number of roots. */
                nroots = ( fammax + iammax + jammax ) / 2 + 1 ;
                /* . Initialize the integral blocks. */
                for ( i = 0 ; i < ( ncfuncf * ncfunci * ncfuncj ) ; i++ ) g[i] = 0.0e+00 ;
                /* . Set the array dimensions. */
                dim1 = fammax + 1 ;
                if ( iammax >= jammax ) dim2 = dim1 * ( jammax + 1 ) ;
                else                    dim2 = dim1 * ( iammax + 1 ) ;
                dim3 = dim1 * ( iammax + 1 ) * ( jammax + 1 ) ;
                /* . Triple loop over primitives. */
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
                        expfac = exp ( - fac ) * PI252 * aainv ;
	                for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( arI[i] + aj * rJ[i] ) * aainv ;
                        /* . Loop over fitting primitives. */
                        for ( fp = 0 ; fp < fBasis->shells[fShell].nprimitives ; fp++ )
                        {
                            /* . Get some information for the primitive. */
	                    expf = fBasis->shells[fShell].primitives[fp].exponent ;
                            /* . Calculate some factors. */
                            ab    = aa * expf ;
                            aandb = aa + expf ;
                            rho   = ab / aandb ;
                            dnuc  = expfac / ( expf * sqrt ( aandb ) ) ;
                            /* . Calculate the rys polynomial roots. */
                            c1x = ( ar[0] - rF[0] ) ;
                            c1y = ( ar[1] - rF[1] ) ;
                            c1z = ( ar[2] - rF[2] ) ;
                            RysQuadrature_Roots ( &roots, nroots, rho * ( c1x * c1x + c1y * c1y + c1z * c1z ) ) ;
                            /* . Calculate some displacements. */
                            axac = aa * ( ar[0] - rC[0] ) ;
                            ayac = aa * ( ar[1] - rC[1] ) ;
                            azac = aa * ( ar[2] - rC[2] ) ;
                            c1x *= aa ;
                            c1y *= aa ;
                            c1z *= aa ;
                            c3x  = expf * ( rF[0] - rC[0] ) + axac ;
                            c3y  = expf * ( rF[1] - rC[1] ) + ayac ;
                            c3z  = expf * ( rF[2] - rC[2] ) + azac ;
                            c4x  = expf * axac ;
                            c4y  = expf * ayac ;
                            c4z  = expf * azac ;
                            /* . Loop over the roots and construct the subsidiary integrals. */
                            for ( m = 0 ; m < nroots ; m++ )
                            {
                                u2    = roots.roots[m] * rho ;
                                f00   = roots.weights[m] ;
                                fac   = 1.0e+00 / ( ab + u2 * aandb ) ;
                                fac2  = 0.5e+00 * fac ;
                                bp01  = ( aa   + u2 ) * fac2 ;
                                b00   =          u2   * fac2 ;
                                b10   = ( expf + u2 ) * fac2 ;
                                xcp00 = u2 * c1x * fac ;
                                ycp00 = u2 * c1y * fac ;
                                zcp00 = u2 * c1z * fac ;
                                xc00  = ( u2 * c3x + c4x ) * fac ;
                                yc00  = ( u2 * c3y + c4y ) * fac ;
                                zc00  = ( u2 * c3z + c4z ) * fac ;
                                Subsidiary_Integral_Nuclear3C ( iammaxt, jammaxt, fammax, QIJ0, QIJ1, QF0, QF1, b00, b10, bp01, dxIJt, dyIJt, dzIJt, f00,
                                                                xc00, xcp00, yc00, ycp00, zc00, zcp00, dim1, dim2, &xint[m*dim3], &yint[m*dim3], &zint[m*dim3] ) ;
                            }
                            /* . Assemble the integrals. */
                            if ( iammax >= jammax )
                            {
                                for ( i = 0, n = 0 ; i < ncfunci ; i++ )
                                {
   	                            ix = CBFPOWX[i+icbfind] * dim2 ;
	                            iy = CBFPOWY[i+icbfind] * dim2 ;
	                            iz = CBFPOWZ[i+icbfind] * dim2 ;
                                    ti = dnuc * iBasis->shells[iShell].primitives[ip].ccbf[i] ;
                                    for ( j = 0 ; j < ncfuncj ; j++ )
                                    {
	                                jix = CBFPOWX[j+jcbfind] * dim1 + ix ;
	                                jiy = CBFPOWY[j+jcbfind] * dim1 + iy ;
	                                jiz = CBFPOWZ[j+jcbfind] * dim1 + iz ;
                                        tij = ti * jBasis->shells[jShell].primitives[jp].ccbf[j] ;
                                        for ( f = 0 ; f < ncfuncf ; f++, n++ )
                                        {
                                            fijx = CBFPOWX[f+fcbfind] + jix ;
                                            fijy = CBFPOWY[f+fcbfind] + jiy ;
                                            fijz = CBFPOWZ[f+fcbfind] + jiz ;
                                            for ( m = 0, fac = 0.0e+00 ; m < nroots ; m++ ) fac += xint[fijx+m*dim3] * yint[fijy+m*dim3] * zint[fijz+m*dim3] ;
                                            tijf = tij * fBasis->shells[fShell].primitives[fp].ccbf[f] ;
                                            g[n] += tijf * fac ;
                                        }
                                    }
                                }
                            }
                            else
                            {
                                for ( j = 0 ; j < ncfuncj ; j++ )
                                {
   	                            ix = CBFPOWX[j+jcbfind] * dim2 ;
	                            iy = CBFPOWY[j+jcbfind] * dim2 ;
	                            iz = CBFPOWZ[j+jcbfind] * dim2 ;
                                    ti = dnuc * jBasis->shells[jShell].primitives[jp].ccbf[j] ;
                                    for ( i = 0 ; i < ncfunci ; i++ )
                                    {
                                        n   = j * ncfuncf + i * ( ncfuncf * ncfuncj ) ;
	                                jix = CBFPOWX[i+icbfind] * dim1 + ix ;
	                                jiy = CBFPOWY[i+icbfind] * dim1 + iy ;
	                                jiz = CBFPOWZ[i+icbfind] * dim1 + iz ;
                                        tij = ti * iBasis->shells[iShell].primitives[ip].ccbf[i] ;
                                        for ( f = 0 ; f < ncfuncf ; f++, n++ )
                                        {
                                            fijx = CBFPOWX[f+fcbfind] + jix ;
                                            fijy = CBFPOWY[f+fcbfind] + jiy ;
                                            fijz = CBFPOWZ[f+fcbfind] + jiz ;
                                            for ( m = 0, fac = 0.0e+00 ; m < nroots ; m++ ) fac += xint[fijx+m*dim3] * yint[fijy+m*dim3] * zint[fijz+m*dim3] ;
                                            tijf = tij * fBasis->shells[fShell].primitives[fp].ccbf[f] ;
                                            g[n] += tijf * fac ;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                /* . Save the integrals. */
                {
                    auto Boolean     skip ;
                    auto Integer     ii, jj, m, n ;
                    auto Cardinal16 *indices16 = block->indices16 ;
                    auto Real       *integrals = block->data      ;
                    m = block->count ;
                    for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nbasisw ; i++ )
                    {
                        ii = iBasis->shells[iShell].nstartw + i ;
                        for ( j = 0 ; j < jBasis->shells[jShell].nbasisw ; j++ )
                        {
                            skip = isDiagonal && ( j > i ) ;
                            jj   = jBasis->shells[jShell].nstartw + j ;
                            for ( f = 0 ; f < fBasis->shells[fShell].nbasisw ; f++, n++ )
                            {
                                if ( ! skip )
                                {
                                    auto Integer  m3 = 3 * m ;
                                    indices16[m3  ] = ii ;
                                    indices16[m3+1] = jj ;
                                    indices16[m3+2] = fBasis->shells[fShell].nstartw + f ;
                                    integrals[m]    = g[n] ;
                                    m++ ;
                                }
                            }
                        }
                    }
                    block->count = m ;
                }
            } /* . fShell. */
        } /* . jShell. */
    } /* . iShell. */
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Electron-fit integral derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisIntegrals_ElectronFitD ( const GaussianBasis *iBasis ,
                                           const Real          *rI     ,
                                           const GaussianBasis *jBasis ,
                                           const Real          *rJ     ,
                                           const Real          *rIJ    ,
                                           const Real           rIJ2   ,
                                           const GaussianBasis *fBasis ,
                                           const Real          *rF     ,
                                                 Block         *block  )
{
    Boolean       iIsJ, isDiagonal, QF0, QF1 ;
    Integer       dim1, ddim2, dim2, ddim3, dim3,
                  f, fammax,          fcbfind, fijx, fijxd, fijy, fijyd, fijz, fijzd, fp, fShell,
                  i, iammax, iammaxt, icbfind, ip, iShell, ix, ixd, iy, iyd, iz, izd,
                  j, jammax, jammaxt, jcbfind, jix, jixd, jiy, jiyd, jiz, jizd, jp, jShell,
                  jUpper, m, n, ncfuncf, ncfunci, ncfuncj, nroots ;
    Real          aa, aandb, aainv, ab, ag, ah, ai, aj, arri, axac, ayac, azac, bp01, b00, b10,
                  c1x, c1y, c1z, c3x, c3y, c3z, c4x, c4y, c4z, dnuc,
                  dxIJt, dyIJt, dzIJt, expfac, expf, fac, facgx, facgy, facgz, fachx, fachy, fachz,
                  fac2, f00, rho, ti, tij, tijf, u2,
                  xc00, xcp00, yc00, ycp00, zc00, zcp00 ;
    Real          ar[3], arI[3],
                  gx[MAXCBF*MAXCBF*MAXCBF], gy[MAXCBF*MAXCBF*MAXCBF], gz[MAXCBF*MAXCBF*MAXCBF],
                  hx[MAXCBF*MAXCBF*MAXCBF], hy[MAXCBF*MAXCBF*MAXCBF], hz[MAXCBF*MAXCBF*MAXCBF],
                  xidg[MAXAMP1*MAXAMP1*MAXAMP1*_MAXRYS], yidg[MAXAMP1*MAXAMP1*MAXAMP1*_MAXRYS], zidg[MAXAMP1*MAXAMP1*MAXAMP1*_MAXRYS] ,
                  xidh[MAXAMP1*MAXAMP1*MAXAMP1*_MAXRYS], yidh[MAXAMP1*MAXAMP1*MAXAMP1*_MAXRYS], zidh[MAXAMP1*MAXAMP1*MAXAMP1*_MAXRYS] ,
                  xint[MAXAMP1*MAXAMP2*MAXAMP2*_MAXRYS], yint[MAXAMP1*MAXAMP2*MAXAMP2*_MAXRYS], zint[MAXAMP1*MAXAMP2*MAXAMP2*_MAXRYS] ;
    const Real   *rC ;
    RysQuadrature roots ;
    /* . Initialization. */
    block->count = 0 ;
    iIsJ         = ( iBasis == jBasis ) && ( rI == rJ ) ;
    /* . Triple loop over shells. */
    for ( iShell = 0 ; iShell < iBasis->nshells ; iShell++ )
    {
        /* . Get information about the shell. */
        iammax  = iBasis->shells[iShell].type->angularmomentum_high ;
        icbfind = iBasis->shells[iShell].type->cbfindex ;
        ncfunci = iBasis->shells[iShell].type->ncbf     ;
        /* . Set the Upper limit for the JSHELL loops. */
        if ( iIsJ ) jUpper = iShell + 1 ;
        else        jUpper = jBasis->nshells ;
        /* . Inner loop over shells. */
        for ( jShell = 0 ; jShell < jUpper ; jShell++ )
        {
            /* . Get information about the shell. */
            jammax  = jBasis->shells[jShell].type->angularmomentum_high ;
            jcbfind = jBasis->shells[jShell].type->cbfindex ;
            ncfuncj = jBasis->shells[jShell].type->ncbf     ;
            /* . Set the diagonal block flag. */
            isDiagonal = iIsJ && ( iShell == jShell ) ;
            /* . Select the expansion center for the recurrence relations. */
            if ( iammax >= jammax )
            {
               iammaxt = iammax ;
               jammaxt = jammax ;
               dxIJt   = rIJ[0] ;
               dyIJt   = rIJ[1] ;
               dzIJt   = rIJ[2] ;
               rC      = rI ;
            }
            else
            {
               iammaxt = jammax ;
               jammaxt = iammax ;
               dxIJt   = - rIJ[0] ;
               dyIJt   = - rIJ[1] ;
               dzIJt   = - rIJ[2] ;
               rC      = rJ ;
            }
            /* . Loop over fitting shells. */
            for ( fShell = 0 ; fShell < fBasis->nshells ; fShell++ )
            {
                /* . Get information about the shell. */
                fammax  = fBasis->shells[fShell].type->angularmomentum_high ;
                fcbfind = fBasis->shells[fShell].type->cbfindex ;
                ncfuncf = fBasis->shells[fShell].type->ncbf     ;
                /* . Set some flags. */
                QF0 = ( fammax == 0 ) ;
                QF1 = ( fammax <= 1 ) ;
                /* . Get the number of roots. */
                nroots = ( fammax + iammax + jammax + 2 ) / 2 + 1 ;
                /* . Initialize the integral blocks. */
                for ( i = 0 ; i < ( ncfuncf * ncfunci * ncfuncj ) ; i++ ) gx[i] = gy[i] = gz[i] = hx[i] = hy[i] = hz[i] = 0.0e+00 ;
                /* . Set the array dimensions. */
                dim1 = fammax + 1 ;
                if ( iammax >= jammax )
                {
                    ddim2 = dim1 * ( jammax + 1 ) ;
                    dim2  = dim1 * ( jammax + 2 ) ;
                }
                else
                {
                    ddim2 = dim1 * ( iammax + 1 ) ;
                    dim2  = dim1 * ( iammax + 2 ) ;
                }
                ddim3 = dim1 * ( iammax + 1 ) * ( jammax + 1 ) ;
                dim3  = dim1 * ( iammax + 2 ) * ( jammax + 2 ) ;
                /* . Triple loop over primitives. */
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
                        expfac = exp ( - fac ) * PI252 * aainv ;
	                for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( arI[i] + aj * rJ[i] ) * aainv ;
                        /* . Set the primitive exponents. */
                        if ( iammax >= jammax )
                        {
                           ag = ai ;
                           ah = aj ;
                        }
                        else
                        {
                           ag = aj ;
                           ah = ai ;
                        }
                        /* . Loop over fitting primitives. */
                        for ( fp = 0 ; fp < fBasis->shells[fShell].nprimitives ; fp++ )
                        {
                            /* . Get some information for the primitive. */
	                    expf = fBasis->shells[fShell].primitives[fp].exponent ;
                            /* . Calculate some factors. */
                            ab    = aa * expf ;
                            aandb = aa + expf ;
                            rho   = ab / aandb ;
                            dnuc  = expfac / ( expf * sqrt ( aandb ) ) ;
                            /* . Calculate the rys polynomial roots. */
                            c1x = ( ar[0] - rF[0] ) ;
                            c1y = ( ar[1] - rF[1] ) ;
                            c1z = ( ar[2] - rF[2] ) ;
                            RysQuadrature_Roots ( &roots, nroots, rho * ( c1x * c1x + c1y * c1y + c1z * c1z ) ) ;
                            /* . Calculate some displacements. */
                            axac = aa * ( ar[0] - rC[0] ) ;
                            ayac = aa * ( ar[1] - rC[1] ) ;
                            azac = aa * ( ar[2] - rC[2] ) ;
                            c1x *= aa ;
                            c1y *= aa ;
                            c1z *= aa ;
                            c3x  = expf * ( rF[0] - rC[0] ) + axac ;
                            c3y  = expf * ( rF[1] - rC[1] ) + ayac ;
                            c3z  = expf * ( rF[2] - rC[2] ) + azac ;
                            c4x  = expf * axac ;
                            c4y  = expf * ayac ;
                            c4z  = expf * azac ;
                            /* . Loop over the roots and construct the subsidiary integrals. */
                            for ( m = 0 ; m < nroots ; m++ )
                            {
                                u2    = roots.roots[m] * rho ;
                                f00   = roots.weights[m] ;
                                fac   = 1.0e+00 / ( ab + u2 * aandb ) ;
                                fac2  = 0.5e+00 * fac ;
                                bp01  = ( aa   + u2 ) * fac2 ;
                                b00   =          u2   * fac2 ;
                                b10   = ( expf + u2 ) * fac2 ;
                                xcp00 = u2 * c1x * fac ;
                                ycp00 = u2 * c1y * fac ;
                                zcp00 = u2 * c1z * fac ;
                                xc00  = ( u2 * c3x + c4x ) * fac ;
                                yc00  = ( u2 * c3y + c4y ) * fac ;
                                zc00  = ( u2 * c3z + c4z ) * fac ;
                                Subsidiary_Integral_Nuclear3C ( iammaxt+1, jammaxt+1, fammax, False, False, QF0, QF1, b00, b10, bp01, dxIJt, dyIJt, dzIJt, f00,
                                                                xc00, xcp00, yc00, ycp00, zc00, zcp00, dim1, dim2, &xint[m*dim3], &yint[m*dim3], &zint[m*dim3] ) ;
                                Subsidiary_Integral_Derivative3 ( &xint[m*dim3], &yint[m*dim3], &zint[m*dim3],
                                                                  &xidg[m*ddim3], &yidg[m*ddim3], &zidg[m*ddim3],
                                                                  &xidh[m*ddim3], &yidh[m*ddim3], &zidh[m*ddim3],
                                                                  ag, ah, iammaxt, jammaxt, fammax, dim1, dim2, dim1, ddim2 ) ;
                            }
                            /* . Assemble the integrals. */
                            if ( iammax >= jammax )
                            {
                                for ( i = 0, n = 0 ; i < ncfunci ; i++ )
                                {
   	                            ix  = CBFPOWX[i+icbfind] * dim2  ;
	                            iy  = CBFPOWY[i+icbfind] * dim2  ;
	                            iz  = CBFPOWZ[i+icbfind] * dim2  ;
   	                            ixd = CBFPOWX[i+icbfind] * ddim2 ;
	                            iyd = CBFPOWY[i+icbfind] * ddim2 ;
	                            izd = CBFPOWZ[i+icbfind] * ddim2 ;
                                    ti  = dnuc * iBasis->shells[iShell].primitives[ip].ccbf[i] ;
                                    for ( j = 0 ; j < ncfuncj ; j++ )
                                    {
	                                jix  = CBFPOWX[j+jcbfind] * dim1 + ix  ;
	                                jiy  = CBFPOWY[j+jcbfind] * dim1 + iy  ;
	                                jiz  = CBFPOWZ[j+jcbfind] * dim1 + iz  ;
	                                jixd = CBFPOWX[j+jcbfind] * dim1 + ixd ;
	                                jiyd = CBFPOWY[j+jcbfind] * dim1 + iyd ;
	                                jizd = CBFPOWZ[j+jcbfind] * dim1 + izd ;
                                        tij  = ti * jBasis->shells[jShell].primitives[jp].ccbf[j] ;
                                        for ( f = 0 ; f < ncfuncf ; f++, n++ )
                                        {
                                            fijx  = CBFPOWX[f+fcbfind] + jix  ;
                                            fijy  = CBFPOWY[f+fcbfind] + jiy  ;
                                            fijz  = CBFPOWZ[f+fcbfind] + jiz  ;
                                            fijxd = CBFPOWX[f+fcbfind] + jixd ;
                                            fijyd = CBFPOWY[f+fcbfind] + jiyd ;
                                            fijzd = CBFPOWZ[f+fcbfind] + jizd ;
                                            for ( m = 0, facgx = facgy = facgz = fachx = fachy = fachz = 0.0e+00 ; m < nroots ; m++ )
                                            {
                                                facgx += xidg[fijxd+m*ddim3] * yint[fijy+m*dim3] * zint[fijz+m*dim3] ;
                                                facgy += xint[fijx+m*dim3] * yidg[fijyd+m*ddim3] * zint[fijz+m*dim3] ;
                                                facgz += xint[fijx+m*dim3] * yint[fijy+m*dim3] * zidg[fijzd+m*ddim3] ;
                                                fachx += xidh[fijxd+m*ddim3] * yint[fijy+m*dim3] * zint[fijz+m*dim3] ;
                                                fachy += xint[fijx+m*dim3] * yidh[fijyd+m*ddim3] * zint[fijz+m*dim3] ;
                                                fachz += xint[fijx+m*dim3] * yint[fijy+m*dim3] * zidh[fijzd+m*ddim3] ;
                                            }
                                            tijf = tij * fBasis->shells[fShell].primitives[fp].ccbf[f] ;
                                            gx[n] += tijf * facgx ;
                                            gy[n] += tijf * facgy ;
                                            gz[n] += tijf * facgz ;
                                            hx[n] += tijf * fachx ;
                                            hy[n] += tijf * fachy ;
                                            hz[n] += tijf * fachz ;
                                        }
                                    }
                                }
                            }
                            else
                            {
                                for ( j = 0 ; j < ncfuncj ; j++ )
                                {
   	                            ix  = CBFPOWX[j+jcbfind] * dim2  ;
	                            iy  = CBFPOWY[j+jcbfind] * dim2  ;
	                            iz  = CBFPOWZ[j+jcbfind] * dim2  ;
    	                            ixd = CBFPOWX[j+jcbfind] * ddim2 ;
	                            iyd = CBFPOWY[j+jcbfind] * ddim2 ;
	                            izd = CBFPOWZ[j+jcbfind] * ddim2 ;
                                    ti  = dnuc * jBasis->shells[jShell].primitives[jp].ccbf[j] ;
                                    for ( i = 0 ; i < ncfunci ; i++ )
                                    {
                                        n    = j * ncfuncf + i * ( ncfuncf * ncfuncj ) ;
	                                jix  = CBFPOWX[i+icbfind] * dim1 + ix  ;
	                                jiy  = CBFPOWY[i+icbfind] * dim1 + iy  ;
	                                jiz  = CBFPOWZ[i+icbfind] * dim1 + iz  ;
	                                jixd = CBFPOWX[i+icbfind] * dim1 + ixd ;
	                                jiyd = CBFPOWY[i+icbfind] * dim1 + iyd ;
	                                jizd = CBFPOWZ[i+icbfind] * dim1 + izd ;
                                        tij  = ti * iBasis->shells[iShell].primitives[ip].ccbf[i] ;
                                        for ( f = 0 ; f < ncfuncf ; f++, n++ )
                                        {
                                            fijx  = CBFPOWX[f+fcbfind] + jix  ;
                                            fijy  = CBFPOWY[f+fcbfind] + jiy  ;
                                            fijz  = CBFPOWZ[f+fcbfind] + jiz  ;
                                            fijxd = CBFPOWX[f+fcbfind] + jixd ;
                                            fijyd = CBFPOWY[f+fcbfind] + jiyd ;
                                            fijzd = CBFPOWZ[f+fcbfind] + jizd ;
                                            for ( m = 0, facgx = facgy = facgz = fachx = fachy = fachz = 0.0e+00 ; m < nroots ; m++ )
                                            {
                                                facgx += xidh[fijxd+m*ddim3] * yint[fijy+m*dim3] * zint[fijz+m*dim3] ;
                                                facgy += xint[fijx+m*dim3] * yidh[fijyd+m*ddim3] * zint[fijz+m*dim3] ;
                                                facgz += xint[fijx+m*dim3] * yint[fijy+m*dim3] * zidh[fijzd+m*ddim3] ;
                                                fachx += xidg[fijxd+m*ddim3] * yint[fijy+m*dim3] * zint[fijz+m*dim3] ;
                                                fachy += xint[fijx+m*dim3] * yidg[fijyd+m*ddim3] * zint[fijz+m*dim3] ;
                                                fachz += xint[fijx+m*dim3] * yint[fijy+m*dim3] * zidg[fijzd+m*ddim3] ;
                                            }
                                            tijf = tij * fBasis->shells[fShell].primitives[fp].ccbf[f] ;
                                            gx[n] += tijf * facgx ;
                                            gy[n] += tijf * facgy ;
                                            gz[n] += tijf * facgz ;
                                            hx[n] += tijf * fachx ;
                                            hy[n] += tijf * fachy ;
                                            hz[n] += tijf * fachz ;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                /* . Save the integrals. */
                {
/*                    auto Boolean     skip ;*/
                    auto Integer     ii, jj, m, m3, m6, n ;
                    auto Cardinal16 *indices16 = block->indices16   ;
                    auto Real       *integrals = block->data, scale ;
                    if ( isDiagonal ) scale = 1.0e+00 ;
                    else              scale = 2.0e+00 ;
                    m = block->count ;
                    for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nbasisw ; i++ )
                    {
                        ii = iBasis->shells[iShell].nstartw + i ;
                        for ( j = 0 ; j < jBasis->shells[jShell].nbasisw ; j++ )
                        {
/*                            skip = isDiagonal && ( j > i ) ;*/
                            jj   = jBasis->shells[jShell].nstartw + j ;
                            for ( f = 0 ; f < fBasis->shells[fShell].nbasisw ; f++, m++, n++ )
                            {
/*
                                if ( ! skip )
                                {
*/
                                    m3              = 3 * m ;
                                    m6              = 6 * m ;
                                    indices16[m3  ] = ii ;
                                    indices16[m3+1] = jj ;
                                    indices16[m3+2] = fBasis->shells[fShell].nstartw + f ;
                                    integrals[m6  ] = scale * gx[n] ;
                                    integrals[m6+1] = scale * gy[n] ;
                                    integrals[m6+2] = scale * gz[n] ;
                                    integrals[m6+3] = scale * hx[n] ;
                                    integrals[m6+4] = scale * hy[n] ;
                                    integrals[m6+5] = scale * hz[n] ;
                                    /*m++ ;*/
/*
                                }
*/
                            }
                        }
                    }
                    block->count = m ;
                }
            } /* . fShell. */
        } /* . jShell. */
    } /* . iShell. */
}
