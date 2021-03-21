/*==================================================================================================================================
! . Integrals - 2 basis, 2 electrons, 1 nucleus/point.
!=================================================================================================================================*/

# include <math.h>
# include <stdlib.h>

# include "GaussianBasisIntegrals_b2e1n1.h"
# include "GaussianNucleus.h"
# include "GaussianBasisSubsidiary.h"
# include "RysQuadrature.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _Selected( selection, i ) ( (selection) == NULL ? True : Block_Item ( selection->flags, i ) )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Electron-nuclear/point integrals.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisIntegrals_ElectronNuclear ( const GaussianBasis *iBasis     ,
                                              const Real          *rI         ,
                                              const GaussianBasis *jBasis     ,
                                              const Real          *rJ         ,
                                              const RealArray1D   *charges    ,
                                              const RealArray1D   *widthsE    ,
                                              const RealArray1D   *widthsN    ,
                                              const Coordinates3  *rNP        ,
                                              const Selection     *selectionN ,
                                                    RealArray2D   *integrals  )
{
    Boolean       iIsJ, QIJ0, QIJ1 ;
    Integer       i, iammax, iammaxt, icbfind, ii, ip, iShell, ix, iy, iz,
                  j, jammax, jammaxt, jcbfind, jdim, jdimm, jj, jp, jShell,
                  jUpper, jxix, jyiy, jziz, k, m, n, ncfunci, ncfuncj, nroots ;
    Real          aa, aandb, aainv, ab, ai, aj, arri, axac, ayac, azac, bp01, b00, b10,
                  c1x, c1y, c1z, c3x, c3y, c3z, c4x, c4y, c4z, dnuc, dxIJt, dyIJt, dzIJt,
                  expfac, expN, fac, facN, fac2, f00, qN, rho, rIJ2, ti, tij, u2,
                  xc00, xcp00, xIJ, yc00, ycp00, yIJ, zc00, zcp00, zIJ ;
    Real          ar[3], ari[3], g[MAXCBF*MAXCBF],
                  xint[MAXAMP1*MAXAMP1*_MAXRYS], yint[MAXAMP1*MAXAMP1*_MAXRYS], zint[MAXAMP1*MAXAMP1*_MAXRYS] ;
    const Real   *rC ;
    Real         *rN ;
    RysQuadrature roots ;
    /* . Initialization. */
    iIsJ = ( iBasis == jBasis ) && ( rI == rJ ) ;
    xIJ  = rI[0] - rJ[0] ;
    yIJ  = rI[1] - rJ[1] ;
    zIJ  = rI[2] - rJ[2] ;
    rIJ2 = xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ;
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
            /* . Get information about the shell. */
            jammax  = jBasis->shells[jShell].type->angularmomentum_high ;
            jdimm   = ( iammax + 1 ) * ( jammax + 1 ) ;
            jcbfind = jBasis->shells[jShell].type->cbfindex ;
            ncfuncj = jBasis->shells[jShell].type->ncbf     ;
            /* . Get the number of roots. */
            nroots = ( iammax + jammax ) / 2 + 1 ;
            /* . Set some flags. */
            QIJ0 = ( iammax + jammax == 0 ) ;
            QIJ1 = ( iammax + jammax <= 1 ) ;
            /* . Initialize the integral blocks. */
            for ( i = 0 ; i < ( ncfunci * ncfuncj ) ; i++ ) g[i] = 0.0e+00 ;
            /* . Select the expansion center for the recurrence relations. */
            if ( iammax >= jammax )
            {
               iammaxt = iammax ;
               jammaxt = jammax ;
               jdim    = jammax + 1 ;
               dxIJt   = xIJ ;
               dyIJt   = yIJ ;
               dzIJt   = zIJ ;
               rC      = rI ;
            }
            else
            {
               iammaxt = jammax ;
               jammaxt = iammax ;
               jdim    = iammax + 1 ;
               dxIJt   = - xIJ ;
               dyIJt   = - yIJ ;
               dzIJt   = - zIJ ;
               rC      = rJ ;
            }
            /* . Outer loop over primitives. */
            for ( ip = 0 ; ip < iBasis->shells[iShell].nprimitives ; ip++ )
            {
                /* . Get some information for the primitive. */
	        ai   = iBasis->shells[iShell].primitives[ip].exponent ;
	        arri = ai * rIJ2 ;
	        for ( i = 0 ; i < 3 ; i++ ) ari[i] = ai * rI[i] ;
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
	            for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( ari[i] + aj * rJ[i] ) * aainv ;
                    /* . Loop over the nuclear densities. */
                    for ( k = 0 ; k < Coordinates3_Rows ( rNP ) ; k++ )
                    {
                        if ( _Selected ( selectionN, k ) )
                        {
                            expN = _GetWidthE ( widthsE, k ) ;
                            facN = _GetWidthN ( widthsN, k ) ;
                            qN   = - Array1D_Item ( charges, k ) ;
                            rN   = Coordinates3_RowPointer ( rNP, k ) ;
                            /* . Calculate some factors. */
                            ab     = aa * expN ;
                            aandb  = aa + expN ;
                            rho    = ab / aandb ;
                            dnuc   = expfac * ( facN * qN ) / ( expN * sqrt ( aandb ) ) ;
                            /* . Calculate the rys polynomial roots. */
                            c1x = ( ar[0] - rN[0] ) ;
                            c1y = ( ar[1] - rN[1] ) ;
                            c1z = ( ar[2] - rN[2] ) ;
                            RysQuadrature_Roots ( &roots, nroots, rho * ( c1x * c1x + c1y * c1y + c1z * c1z ) ) ;
                            /* . Calculate some displacements. */
                            axac = aa * ( ar[0] - rC[0] ) ;
                            ayac = aa * ( ar[1] - rC[1] ) ;
                            azac = aa * ( ar[2] - rC[2] ) ;
                            c1x *= aa ;
                            c1y *= aa ;
                            c1z *= aa ;
                            c3x  = expN * ( rN[0] - rC[0] ) + axac ;
                            c3y  = expN * ( rN[1] - rC[1] ) + ayac ;
                            c3z  = expN * ( rN[2] - rC[2] ) + azac ;
                            c4x  = expN * axac ;
                            c4y  = expN * ayac ;
                            c4z  = expN * azac ;
                            /* . Loop over the roots and construct the subsidiary integrals. */
                            for ( m = 0 ; m < nroots ; m++ )
                            {
                                u2    = roots.roots[m] * rho ;
                                f00   = roots.weights[m] ;
                                fac   = 1.0e+00 / ( ab + u2 * aandb ) ;
                                fac2  = 0.5e+00 * fac ;
                                bp01  = ( aa   + u2 ) * fac2 ;
                                b00   =          u2   * fac2 ;
                                b10   = ( expN + u2 ) * fac2 ;
                                xcp00 = u2 * c1x * fac ;
                                ycp00 = u2 * c1y * fac ;
                                zcp00 = u2 * c1z * fac ;
                                xc00  = ( u2 * c3x + c4x ) * fac ;
                                yc00  = ( u2 * c3y + c4y ) * fac ;
                                zc00  = ( u2 * c3z + c4z ) * fac ;
                                Subsidiary_Integral_Nuclear3C ( iammaxt, jammaxt, 0, QIJ0, QIJ1, True, True, b00, b10, bp01, dxIJt, dyIJt,dzIJt, f00,
                                                                xc00, xcp00, yc00, ycp00, zc00, zcp00, 1, jdim, &xint[m*jdimm], &yint[m*jdimm], &zint[m*jdimm] ) ;
                            }
                            /* . Assemble the integrals. */
                            if ( iammax >= jammax )
                            {
                                for ( i = 0, n = 0 ; i < ncfunci ; i++ )
                                {
   	                            ix = CBFPOWX[i+icbfind] * jdim ;
	                            iy = CBFPOWY[i+icbfind] * jdim ;
	                            iz = CBFPOWZ[i+icbfind] * jdim ;
                                    ti = dnuc * iBasis->shells[iShell].primitives[ip].ccbf[i] ;
                                    for ( j = 0 ; j < ncfuncj ; j++, n++ )
                                    {
	                                jxix = CBFPOWX[j+jcbfind] + ix ;
	                                jyiy = CBFPOWY[j+jcbfind] + iy ;
	                                jziz = CBFPOWZ[j+jcbfind] + iz ;
                                        for ( m = 0, fac = 0.0e+00 ; m < nroots ; m++ ) fac += xint[jxix+m*jdimm] * yint[jyiy+m*jdimm] * zint[jziz+m*jdimm] ;
                                        tij = ti * jBasis->shells[jShell].primitives[jp].ccbf[j] ;
                                        g[n] += tij * fac ;
                                    }
                                }
                            }
                            else
                            {
                                for ( j = 0 ; j < ncfuncj ; j++ )
                                {
   	                            ix = CBFPOWX[j+jcbfind] * jdim ;
	                            iy = CBFPOWY[j+jcbfind] * jdim ;
	                            iz = CBFPOWZ[j+jcbfind] * jdim ;
                                    ti = dnuc * jBasis->shells[jShell].primitives[jp].ccbf[j] ;
                                    for ( i = 0, n = j ; i < ncfunci ; i++, n+= ncfuncj )
                                    {
	                                jxix = CBFPOWX[i+icbfind] + ix ;
	                                jyiy = CBFPOWY[i+icbfind] + iy ;
	                                jziz = CBFPOWZ[i+icbfind] + iz ;
                                        for ( m = 0, fac = 0.0e+00 ; m < nroots ; m++ ) fac += xint[jxix+m*jdimm] * yint[jyiy+m*jdimm] * zint[jziz+m*jdimm] ;
                                        tij = ti * iBasis->shells[iShell].primitives[ip].ccbf[i] ;
                                        g[n] += tij * fac ;
                                    }
                                }
                            }
                        } /* . Selected. */
                    } /* . k. */
                } /* . jp. */
            } /* . ip. */
            /* . Save the integrals. */
            for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nbasisw ; i++ )
            {
                ii = iBasis->shells[iShell].nstartw + i ;
                for ( j = 0 ; j < jBasis->shells[jShell].nbasisw ; j++, n++ )
                {
                    jj = jBasis->shells[jShell].nstartw + j ;
                    Array2D_Item ( integrals, ii, jj ) = g[n] ;
	        }
            }
        } /* . jShell. */
    } /* . iShell. */
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Electron-nuclear/point derivatives.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void GaussianBasisIntegrals_ElectronNuclearD ( const GaussianBasis *iBasis     ,
                                                      const Real          *rI         ,
                                                      const GaussianBasis *jBasis     ,
                                                      const Real          *rJ         ,
                                                      const RealArray1D   *charges    ,
                                                      const RealArray1D   *widthsE    ,
                                                      const RealArray1D   *widthsN    ,
                                                      const Coordinates3  *rNP        ,
                                                      const Selection     *selectionN ,
                                                      const RealArray2D   *dOneIJ     ,
                                                            Real          *dRi        ,
                                                            Real          *dRj        ,
                                                            Coordinates3  *gN         )

{
    Boolean       iIsJ, isDiagonal ;
    Integer       ddim1, ddim2,
                  i, iammax, iammaxt, icbfind, ii, ip, iShell, ix, ixd, iy, iyd, iz, izd,
                  j, jammax, jammaxt, jcbfind, jdim, jdimm, jj, jp, jShell, jUpper,
                  jxix, jxixd, jyiy, jyiyd, jziz, jzizd, k, m, n, ncfunci, ncfuncj, nroots ;
    Real          aa, aandb, aainv, ab, ag, ah, ai, aj, arri, axac, ayac, azac, bp01, b00, b10,
                  c1x, c1y, c1z, c3x, c3y, c3z, c4x, c4y, c4z, dGx, dGy, dGz, dHx, dHy, dHz, dnuc,
                  dxIJt, dyIJt, dzIJt, expfac, expN, fac, facgx, facgy, facgz, fachx, fachy, fachz,
                  facN, fac2, f00, qN, rho, rIJ2, scale, ti, tij, u2,
                  xc00, xcp00, xIJ, yc00, ycp00, yIJ, zc00, zcp00, zIJ ;
    Real          ar[3], ari[3],
                  gx[MAXCBF*MAXCBF], gy[MAXCBF*MAXCBF], gz[MAXCBF*MAXCBF],
                  hx[MAXCBF*MAXCBF], hy[MAXCBF*MAXCBF], hz[MAXCBF*MAXCBF],
                  xidG[MAXAMP1*MAXAMP1*_MAXRYS], yidG[MAXAMP1*MAXAMP1*_MAXRYS], zidG[MAXAMP1*MAXAMP1*_MAXRYS] ,
                  xidH[MAXAMP1*MAXAMP1*_MAXRYS], yidH[MAXAMP1*MAXAMP1*_MAXRYS], zidH[MAXAMP1*MAXAMP1*_MAXRYS] ,
                  xint[MAXAMP2*MAXAMP2*_MAXRYS], yint[MAXAMP2*MAXAMP2*_MAXRYS], zint[MAXAMP2*MAXAMP2*_MAXRYS] ;
    const Real   *rC ;
    Real         *rN ;
    RysQuadrature roots ;
    /* . Initialization. */
    iIsJ = ( iBasis == jBasis ) && ( rI == rJ ) ;
    xIJ  = rI[0] - rJ[0] ;
    yIJ  = rI[1] - rJ[1] ;
    zIJ  = rI[2] - rJ[2] ;
    rIJ2 = xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ;
    for ( i = 0 ; i < 3 ; i++ ) { dRi[i] = dRj[i] = 0.0e+00 ; }
    /* . Loop over the nuclear densities. */
    for ( k = 0 ; k < Coordinates3_Rows ( rNP ) ; k++ )
    {
        if ( _Selected ( selectionN, k ) )
        {
            expN = _GetWidthE ( widthsE, k ) ;
            facN = _GetWidthN ( widthsN, k ) ;
            qN   = - Array1D_Item ( charges, k ) ;
            rN   = Coordinates3_RowPointer ( rNP, k ) ;
            /* . Initialize some accumulators. */
            dGx = dGy = dGz = dHx = dHy = dHz = 0.0e+00 ;
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
                    ddim2   = ( iammax + 1 ) * ( jammax + 1 ) ;
                    jdimm   = ( iammax + 2 ) * ( jammax + 2 ) ;
                    jcbfind = jBasis->shells[jShell].type->cbfindex ;
                    ncfuncj = jBasis->shells[jShell].type->ncbf     ;
                    /* . Set the diagonal block flag. */
                    isDiagonal = iIsJ && ( iShell == jShell ) ;
                    /* . Get the number of roots. */
                    nroots = ( iammax + jammax + 2 ) / 2 + 1 ;
                    /* . Initialize the integral blocks. */
                    for ( i = 0 ; i < ( ncfunci * ncfuncj ) ; i++ ) gx[i] = gy[i] = gz[i] = hx[i] = hy[i] = hz[i] = 0.0e+00 ;
                    /* . Select the expansion center for the recurrence relations. */
                    if ( iammax >= jammax )
                    {
                       iammaxt = iammax ;
                       jammaxt = jammax ;
                       ddim1   = jammax + 1 ;
                       jdim    = jammax + 2 ;
                       dxIJt   = xIJ ;
                       dyIJt   = yIJ ;
                       dzIJt   = zIJ ;
                       rC      = rI ;
                    }
                    else
                    {
                       iammaxt = jammax ;
                       jammaxt = iammax ;
                       ddim1   = iammax + 1 ;
                       jdim    = iammax + 2 ;
                       dxIJt   = - xIJ ;
                       dyIJt   = - yIJ ;
                       dzIJt   = - zIJ ;
                       rC      = rJ ;
                    }
                    /* . Outer loop over primitives. */
                    for ( ip = 0 ; ip < iBasis->shells[iShell].nprimitives ; ip++ )
                    {
                        /* . Get some information for the primitive. */
	                ai   = iBasis->shells[iShell].primitives[ip].exponent ;
	                arri = ai * rIJ2 ;
	                for ( i = 0 ; i < 3 ; i++ ) ari[i] = ai * rI[i] ;
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
	                    for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( ari[i] + aj * rJ[i] ) * aainv ;
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
                            /* . Calculate some factors. */
                            ab    = aa * expN ;
                            aandb = aa + expN ;
                            rho   = ab / aandb ;
                            dnuc  = expfac * ( facN * qN ) / ( expN * sqrt ( aandb ) ) ;
                            /* . Calculate the rys polynomial roots. */
                            c1x = ( ar[0] - rN[0] ) ;
                            c1y = ( ar[1] - rN[1] ) ;
                            c1z = ( ar[2] - rN[2] ) ;
                            RysQuadrature_Roots ( &roots, nroots, rho * ( c1x * c1x + c1y * c1y + c1z * c1z ) ) ;
                            /* . Calculate some displacements. */
                            axac = aa * ( ar[0] - rC[0] ) ;
                            ayac = aa * ( ar[1] - rC[1] ) ;
                            azac = aa * ( ar[2] - rC[2] ) ;
                            c1x *= aa ;
                            c1y *= aa ;
                            c1z *= aa ;
                            c3x  = expN * ( rN[0] - rC[0] ) + axac ;
                            c3y  = expN * ( rN[1] - rC[1] ) + ayac ;
                            c3z  = expN * ( rN[2] - rC[2] ) + azac ;
                            c4x  = expN * axac ;
                            c4y  = expN * ayac ;
                            c4z  = expN * azac ;
                            /* . Loop over the roots and construct the subsidiary integrals. */
                            for ( m = 0 ; m < nroots ; m++ )
                            {
                                u2    = roots.roots[m] * rho ;
                                f00   = roots.weights[m] ;
                                fac   = 1.0e+00 / ( ab + u2 * aandb ) ;
                                fac2  = 0.5e+00 * fac ;
                                bp01  = ( aa   + u2 ) * fac2 ;
                                b00   =          u2   * fac2 ;
                                b10   = ( expN + u2 ) * fac2 ;
                                xcp00 = u2 * c1x * fac ;
                                ycp00 = u2 * c1y * fac ;
                                zcp00 = u2 * c1z * fac ;
                                xc00  = ( u2 * c3x + c4x ) * fac ;
                                yc00  = ( u2 * c3y + c4y ) * fac ;
                                zc00  = ( u2 * c3z + c4z ) * fac ;
                                Subsidiary_Integral_Nuclear3C ( iammaxt+1, jammaxt+1, 0, False, False, True, True, b00, b10, bp01, dxIJt, dyIJt, dzIJt, f00,
                                                                xc00, xcp00, yc00, ycp00, zc00, zcp00, 1, jdim, &xint[m*jdimm], &yint[m*jdimm], &zint[m*jdimm] ) ;
                                Subsidiary_Integral_Derivative3 ( &xint[m*jdimm], &yint[m*jdimm], &zint[m*jdimm],
                                                                  &xidG[m*ddim2], &yidG[m*ddim2], &zidG[m*ddim2],
                                                                  &xidH[m*ddim2], &yidH[m*ddim2], &zidH[m*ddim2],
                                                                  ag, ah, iammaxt, jammaxt, 0, 1, jdim, 1, ddim1 ) ;
                            }
                            /* . Assemble the integrals. */
                            if ( iammax >= jammax )
                            {
                                for ( i = 0, n = 0 ; i < ncfunci ; i++ )
                                {
   	                            ix  = CBFPOWX[i+icbfind] * jdim  ;
	                            iy  = CBFPOWY[i+icbfind] * jdim  ;
	                            iz  = CBFPOWZ[i+icbfind] * jdim  ;
   	                            ixd = CBFPOWX[i+icbfind] * ddim1 ;
	                            iyd = CBFPOWY[i+icbfind] * ddim1 ;
	                            izd = CBFPOWZ[i+icbfind] * ddim1 ;
                                    ti  = dnuc * iBasis->shells[iShell].primitives[ip].ccbf[i] ;
                                    for ( j = 0 ; j < ncfuncj ; j++, n++ )
                                    {
	                                jxix  = CBFPOWX[j+jcbfind] + ix  ;
	                                jyiy  = CBFPOWY[j+jcbfind] + iy  ;
	                                jziz  = CBFPOWZ[j+jcbfind] + iz  ;
	                                jxixd = CBFPOWX[j+jcbfind] + ixd ;
	                                jyiyd = CBFPOWY[j+jcbfind] + iyd ;
	                                jzizd = CBFPOWZ[j+jcbfind] + izd ;
                                        for ( m = 0, facgx = facgy = facgz = fachx = fachy = fachz = 0.0e+00 ; m < nroots ; m++ )
                                        {
                                            facgx += xidG[jxixd+m*ddim2] * yint[jyiy+m*jdimm] * zint[jziz+m*jdimm] ;
                                            facgy += xint[jxix+m*jdimm] * yidG[jyiyd+m*ddim2] * zint[jziz+m*jdimm] ;
                                            facgz += xint[jxix+m*jdimm] * yint[jyiy+m*jdimm] * zidG[jzizd+m*ddim2] ;
                                            fachx += xidH[jxixd+m*ddim2] * yint[jyiy+m*jdimm] * zint[jziz+m*jdimm] ;
                                            fachy += xint[jxix+m*jdimm] * yidH[jyiyd+m*ddim2] * zint[jziz+m*jdimm] ;
                                            fachz += xint[jxix+m*jdimm] * yint[jyiy+m*jdimm] * zidH[jzizd+m*ddim2] ;
                                        }
                                        tij = ti * jBasis->shells[jShell].primitives[jp].ccbf[j] ;
                                        gx[n] += tij * facgx ;
                                        gy[n] += tij * facgy ;
                                        gz[n] += tij * facgz ;
                                        hx[n] += tij * fachx ;
                                        hy[n] += tij * fachy ;
                                        hz[n] += tij * fachz ;
                                    }
                                }
                            }
                            else
                            {
                                for ( j = 0 ; j < ncfuncj ; j++ )
                                {
   	                            ix  = CBFPOWX[j+jcbfind] * jdim  ;
	                            iy  = CBFPOWY[j+jcbfind] * jdim  ;
	                            iz  = CBFPOWZ[j+jcbfind] * jdim  ;
   	                            ixd = CBFPOWX[j+jcbfind] * ddim1 ;
	                            iyd = CBFPOWY[j+jcbfind] * ddim1 ;
	                            izd = CBFPOWZ[j+jcbfind] * ddim1 ;
                                    ti  = dnuc * jBasis->shells[jShell].primitives[jp].ccbf[j] ;
                                    for ( i = 0, n = j ; i < ncfunci ; i++, n+= ncfuncj )
                                    {
	                                jxix  = CBFPOWX[i+icbfind] + ix  ;
	                                jyiy  = CBFPOWY[i+icbfind] + iy  ;
	                                jziz  = CBFPOWZ[i+icbfind] + iz  ;
	                                jxixd = CBFPOWX[i+icbfind] + ixd ;
	                                jyiyd = CBFPOWY[i+icbfind] + iyd ;
	                                jzizd = CBFPOWZ[i+icbfind] + izd ;
                                        for ( m = 0, facgx = facgy = facgz = fachx = fachy = fachz = 0.0e+00 ; m < nroots ; m++ )
                                        {
                                            facgx += xidH[jxixd+m*ddim2] * yint[jyiy+m*jdimm] * zint[jziz+m*jdimm] ;
                                            facgy += xint[jxix+m*jdimm] * yidH[jyiyd+m*ddim2] * zint[jziz+m*jdimm] ;
                                            facgz += xint[jxix+m*jdimm] * yint[jyiy+m*jdimm] * zidH[jzizd+m*ddim2] ;
                                            fachx += xidG[jxixd+m*ddim2] * yint[jyiy+m*jdimm] * zint[jziz+m*jdimm] ;
                                            fachy += xint[jxix+m*jdimm] * yidG[jyiyd+m*ddim2] * zint[jziz+m*jdimm] ;
                                            fachz += xint[jxix+m*jdimm] * yint[jyiy+m*jdimm] * zidG[jzizd+m*ddim2] ;
                                        }
                                        tij = ti * iBasis->shells[iShell].primitives[ip].ccbf[i] ;
                                        gx[n] += tij * facgx ;
                                        gy[n] += tij * facgy ;
                                        gz[n] += tij * facgz ;
                                        hx[n] += tij * fachx ;
                                        hy[n] += tij * fachy ;
                                        hz[n] += tij * fachz ;
                                    }
                                }
                            }
                        } /* . jP. */
                    } /* . iP. */
                    /* . Get the scale factor. */
                    if ( isDiagonal ) scale = 1.0e+00 ;
                    else              scale = 2.0e+00 ;
                    /* .  Add in the blocks of integrals to the derivatives. */
                    for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nbasisw ; i++ )
                    {
                        ii = iBasis->shells[iShell].nstartw + i ;
                        for ( j = 0 ; j < jBasis->shells[jShell].nbasisw ; j++, n++ )
                        {
                            jj   = jBasis->shells[jShell].nstartw + j ;
                            fac  = scale * Array2D_Item ( dOneIJ, ii, jj ) ;
                            dGx += fac * gx[n] ;
                            dGy += fac * gy[n] ;
                            dGz += fac * gz[n] ;
                            dHx += fac * hx[n] ;
                            dHy += fac * hy[n] ;
                            dHz += fac * hz[n] ;
	                }
                    }
                } /* . jShell. */
            } /* . iShell. */
            /* . Sum in the contributions to the gradients. */
            dRi[0] += dGx ; dRi[1] += dGy ; dRi[2] += dGz ;
            dRj[0] += dHx ; dRj[1] += dHy ; dRj[2] += dHz ;
            Coordinates3_DecrementRow ( gN, k, dGx, dGy, dGz ) ;
            Coordinates3_DecrementRow ( gN, k, dHx, dHy, dHz ) ;
        } /* . Selected. */
    } /* . k. */
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Electron-nuclear/point potentials.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisIntegrals_ElectronNuclearPotentials ( const GaussianBasis *iBasis     ,
                                                        const Real          *rI         ,
                                                        const GaussianBasis *jBasis     ,
                                                        const Real          *rJ         ,
                                                        const RealArray1D   *widthsE    ,
                                                        const RealArray1D   *widthsN    ,
                                                        const Coordinates3  *rNP        ,
                                                        const Selection     *selectionN ,
                                                        const RealArray2D   *dOneIJ     ,
                                                              RealArray1D   *potentials )
{
    Boolean       iIsJ, isDiagonal, QIJ0, QIJ1 ;
    Integer       i, iammax, iammaxt, icbfind, ii, ip, iShell, ix, iy, iz,
                  j, jammax, jammaxt, jcbfind, jdim, jdimm, jj, jp, jShell,
                  jUpper, jxix, jyiy, jziz, k, m, n, ncfunci, ncfuncj, nroots ;
    Real          aa, aandb, aainv, ab, ai, aj, arri, axac, ayac, azac, bp01, b00, b10, c1x, c1y, c1z, c3x, c3y, c3z, c4x, c4y, c4z,
                  dnuc, dxIJt, dyIJt, dzIJt, expfac, expN, fac, facN, fac2, f00, pot, rho, rIJ2, scale, ti, tij, u2,
                  xc00, xcp00, xIJ, yc00, ycp00, yIJ, zc00, zcp00, zIJ ;
    Real          ar[3], ari[3], g[MAXCBF*MAXCBF],
                  xint[MAXAMP1*MAXAMP1*_MAXRYS], yint[MAXAMP1*MAXAMP1*_MAXRYS], zint[MAXAMP1*MAXAMP1*_MAXRYS] ;
    const Real   *rC ;
    Real         *rN ;
    RysQuadrature roots ;
    /* . Initialization. */
    iIsJ = ( iBasis == jBasis ) && ( rI == rJ ) ;
    xIJ  = rI[0] - rJ[0] ;
    yIJ  = rI[1] - rJ[1] ;
    zIJ  = rI[2] - rJ[2] ;
    rIJ2 = xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ;
    /* . Loop over the points. */
    for ( k = 0 ; k < Coordinates3_Rows ( rNP ) ; k++ )
    {
        if ( _Selected ( selectionN, k ) )
        {
            expN = _GetWidthE ( widthsE, k ) ;
            facN = _GetWidthN ( widthsN, k ) ;
            rN   = Coordinates3_RowPointer ( rNP, k ) ;
            pot  = 0.0e+00 ;
            /* . Outer loop over shells. */
            for ( iShell = 0 ; iShell < iBasis->nshells ; iShell++ )
            {
                iammax  = iBasis->shells[iShell].type->angularmomentum_high ;
                icbfind = iBasis->shells[iShell].type->cbfindex ;
                ncfunci = iBasis->shells[iShell].type->ncbf     ;
                /* . Inner loop over shells. */
                if ( iIsJ ) jUpper = iShell + 1 ;
                else           jUpper = jBasis->nshells ;
                for ( jShell = 0 ; jShell < jUpper ; jShell++ )
                {
                    jammax  = jBasis->shells[jShell].type->angularmomentum_high ;
                    jdimm   = ( iammax + 1 ) * ( jammax + 1 ) ;
                    jcbfind = jBasis->shells[jShell].type->cbfindex ;
                    ncfuncj = jBasis->shells[jShell].type->ncbf     ;
                    /* . Set the diagonal block flag. */
                    isDiagonal = iIsJ && ( iShell == jShell ) ;
                    /* . Get the number of roots. */
                    nroots = ( iammax + jammax ) / 2 + 1 ;
                    /* . Set some flags. */
                    QIJ0 = ( iammax + jammax == 0 ) ;
                    QIJ1 = ( iammax + jammax <= 1 ) ;
                    /* . Initialize the integral blocks. */
                    for ( i = 0 ; i < ( ncfunci * ncfuncj ) ; i++ ) g[i] = 0.0e+00 ;
                    /* . Select the expansion center for the recurrence relations. */
                    if ( iammax >= jammax )
                    {
                       iammaxt = iammax ;
                       jammaxt = jammax ;
                       jdim    = jammax + 1 ;
                       dxIJt   = xIJ ;
                       dyIJt   = yIJ ;
                       dzIJt   = zIJ ;
                       rC      = rI ;
                    }
                    else
                    {
                       iammaxt = jammax ;
                       jammaxt = iammax ;
                       jdim    = iammax + 1 ;
                       dxIJt   = - xIJ ;
                       dyIJt   = - yIJ ;
                       dzIJt   = - zIJ ;
                       rC      = rJ ;
                    }
                    /* . Outer loop over primitives. */
                    for ( ip = 0 ; ip < iBasis->shells[iShell].nprimitives ; ip++ )
                    {
                        /* . Get some information for the primitive. */
	                ai   = iBasis->shells[iShell].primitives[ip].exponent ;
	                arri = ai * rIJ2 ;
	                for ( i = 0 ; i < 3 ; i++ ) ari[i] = ai * rI[i] ;
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
	                    for ( i = 0 ; i < 3 ; i++ ) ar[i] = ( ari[i] + aj * rJ[i] ) * aainv ;
                            /* . Start of point-specific code. */
                            /* . Calculate some factors. */
                            ab     = aa * expN ;
                            aandb  = aa + expN ;
                            rho    = ab / aandb ;
                            dnuc   = expfac * facN / ( expN * sqrt ( aandb ) ) ;
                            /* . Calculate the rys polynomial roots. */
                            c1x = ( ar[0] - rN[0] ) ;
                            c1y = ( ar[1] - rN[1] ) ;
                            c1z = ( ar[2] - rN[2] ) ;
                            RysQuadrature_Roots ( &roots, nroots, rho * ( c1x * c1x + c1y * c1y + c1z * c1z ) ) ;
                            /* . Calculate some displacements. */
                            axac = aa * ( ar[0] - rC[0] ) ;
                            ayac = aa * ( ar[1] - rC[1] ) ;
                            azac = aa * ( ar[2] - rC[2] ) ;
                            c1x *= aa ;
                            c1y *= aa ;
                            c1z *= aa ;
                            c3x  = expN * ( rN[0] - rC[0] ) + axac ;
                            c3y  = expN * ( rN[1] - rC[1] ) + ayac ;
                            c3z  = expN * ( rN[2] - rC[2] ) + azac ;
                            c4x  = expN * axac ;
                            c4y  = expN * ayac ;
                            c4z  = expN * azac ;
                            /* . Loop over the roots and construct the subsidiary integrals. */
                            for ( m = 0 ; m < nroots ; m++ )
                            {
                                u2    = roots.roots[m] * rho ;
                                f00   = roots.weights[m] ;
                                fac   = 1.0e+00 / ( ab + u2 * aandb ) ;
                                fac2  = 0.5e+00 * fac ;
                                bp01  = ( aa   + u2 ) * fac2 ;
                                b00   =          u2   * fac2 ;
                                b10   = ( expN + u2 ) * fac2 ;
                                xcp00 = u2 * c1x * fac ;
                                ycp00 = u2 * c1y * fac ;
                                zcp00 = u2 * c1z * fac ;
                                xc00  = ( u2 * c3x + c4x ) * fac ;
                                yc00  = ( u2 * c3y + c4y ) * fac ;
                                zc00  = ( u2 * c3z + c4z ) * fac ;
                                Subsidiary_Integral_Nuclear3C ( iammaxt, jammaxt, 0, QIJ0, QIJ1, True, True, b00, b10, bp01, dxIJt, dyIJt,dzIJt, f00,
                                                                xc00, xcp00, yc00, ycp00, zc00, zcp00, 1, jdim, &xint[m*jdimm], &yint[m*jdimm], &zint[m*jdimm] ) ;
                            }
                            /* . Assemble the integrals. */
                            if ( iammax >= jammax )
                            {
                                for ( i = 0, n = 0 ; i < ncfunci ; i++ )
                                {
   	                            ix = CBFPOWX[i+icbfind] * jdim ;
	                            iy = CBFPOWY[i+icbfind] * jdim ;
	                            iz = CBFPOWZ[i+icbfind] * jdim ;
                                    ti = dnuc * iBasis->shells[iShell].primitives[ip].ccbf[i] ;
                                    for ( j = 0 ; j < ncfuncj ; j++, n++ )
                                    {
	                                jxix = CBFPOWX[j+jcbfind] + ix ;
	                                jyiy = CBFPOWY[j+jcbfind] + iy ;
	                                jziz = CBFPOWZ[j+jcbfind] + iz ;
                                        for ( m = 0, fac = 0.0e+00 ; m < nroots ; m++ ) fac += xint[jxix+m*jdimm] * yint[jyiy+m*jdimm] * zint[jziz+m*jdimm] ;
                                        tij = ti * jBasis->shells[jShell].primitives[jp].ccbf[j] ;
                                        g[n] += tij * fac ;
                                    }
                                }
                            }
                            else
                            {
                                for ( j = 0 ; j < ncfuncj ; j++ )
                                {
   	                            ix = CBFPOWX[j+jcbfind] * jdim ;
	                            iy = CBFPOWY[j+jcbfind] * jdim ;
	                            iz = CBFPOWZ[j+jcbfind] * jdim ;
                                    ti = dnuc * jBasis->shells[jShell].primitives[jp].ccbf[j] ;
                                    for ( i = 0, n = j ; i < ncfunci ; i++, n+= ncfuncj )
                                    {
	                                jxix = CBFPOWX[i+icbfind] + ix ;
	                                jyiy = CBFPOWY[i+icbfind] + iy ;
	                                jziz = CBFPOWZ[i+icbfind] + iz ;
                                        for ( m = 0, fac = 0.0e+00 ; m < nroots ; m++ ) fac += xint[jxix+m*jdimm] * yint[jyiy+m*jdimm] * zint[jziz+m*jdimm] ;
                                        tij = ti * iBasis->shells[iShell].primitives[ip].ccbf[i] ;
                                        g[n] += tij * fac ;
                                    }
                                }
                            }
                        } /* . jP. */
                    } /* . iP. */
                    /* . Determine scaling. */
                    if ( isDiagonal ) scale = 1.0e+00 ;
                    else              scale = 2.0e+00 ;
                    /* .  Add in the block of integrals to the potential - i is usually greater than j. */
                    for ( i = 0, n = 0 ; i < iBasis->shells[iShell].nbasisw ; i++ )
                    {
                        ii = iBasis->shells[iShell].nstartw + i ;
                        for ( j = 0 ; j < jBasis->shells[jShell].nbasisw ; j++, n++ )
                        {
                            jj = jBasis->shells[jShell].nstartw + j ;
                            pot += ( Array2D_Item ( dOneIJ, ii, jj ) * g[n] * scale ) ;
	                }
                    }
                } /* . jShell. */
            } /* . iShell. */
            /* . Save the potential. */
            Array1D_Item ( potentials, k ) -= pot ;
        } /* . Selected. */
    } /* . k. */
}

# undef _Selected
