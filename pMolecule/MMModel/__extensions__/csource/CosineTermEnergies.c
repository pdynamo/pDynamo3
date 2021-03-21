/*==================================================================================================================================
! . Energies and gradients for MM terms expressed as cosine expansions - sum_p c_p * cos ( p x ).
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "Boolean.h"
# include "CosineTermEnergies.h"
# include "Integer.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Angle i-j-k.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real CosineTermEnergy_Angle ( const CosineTermContainer *self         ,
                              const Coordinates3        *coordinates3 ,
                                    Coordinates3        *gradients3   )
{
    Real energy = 0.0e+00 ;
    if ( ( self != NULL ) && ( coordinates3 != NULL ) )
    {
	auto Boolean doGradients ;
        auto Integer i, j, k, nt, p, t ;
        auto Real    c, cn, co, cosPhi, dF, dTxi, dTyi, dTzi, dTxj, dTyj, dTzj, dTxk, dTyk, dTzk,
                     e, rij, rkj, xij, yij, zij, xkj, ykj, zkj ;
	doGradients = ( gradients3 != NULL ) ;
	for ( nt = 0 ; nt < self->nTerms ; nt++ )
	{
	    if ( self->terms[nt].isActive )
	    {
                /* . Local data. */
	        i = self->terms[nt].indices[0] ;
	        j = self->terms[nt].indices[1] ;
	        k = self->terms[nt].indices[2] ;
	        t = self->terms[nt].type ;
	        /* . Coordinate displacements. */
	        Coordinates3_DifferenceRow ( coordinates3, i, j, xij, yij, zij ) ;
	        Coordinates3_DifferenceRow ( coordinates3, k, j, xkj, ykj, zkj ) ;
                /* . Normalize. */
                rij  = sqrt ( xij * xij + yij * yij + zij * zij ) ;
                xij /= rij ; yij /= rij ; zij /= rij ;
                rkj  = sqrt ( xkj * xkj + ykj * ykj + zkj * zkj ) ;
                xkj /= rkj ; ykj /= rkj ; zkj /= rkj ;
                /* . Cosine of the angle. */
                cosPhi = xij * xkj + yij * ykj + zij * zkj ;
                /* . Loop over powers of the cosine. */
                cn = 1.0e+00 ;
                co = 0.0e+00 ;
                dF = 0.0e+00 ;
                e  = 0.0e+00 ;
                for ( p = 0 ; p <= self->parameters[t].nPowers ; p++ )
                {
                    c   = self->parameters[t].powerCoefficients[p] ;
                    dF += c * co * ( Real ) p ;
                    e  += c * cn ;
                    co  = cn ;
                    cn *= cosPhi ;
                }
                /* . The energy term. */
                energy += e ;
                if ( doGradients )
                {
                    /* . i terms. */
                    dTxi = dF * ( xkj - cosPhi * xij ) / rij ;
                    dTyi = dF * ( ykj - cosPhi * yij ) / rij ;
                    dTzi = dF * ( zkj - cosPhi * zij ) / rij ;
                    /* . k terms. */
                    dTxk = dF * ( xij - cosPhi * xkj ) / rkj ;
                    dTyk = dF * ( yij - cosPhi * ykj ) / rkj ;
                    dTzk = dF * ( zij - cosPhi * zkj ) / rkj ;
                    /* . j terms. */
	            dTxj  = - dTxi - dTxk ;
	            dTyj  = - dTyi - dTyk ;
	            dTzj  = - dTzi - dTzk ;
                    /* . Add in the contributions. */
	            Coordinates3_IncrementRow ( gradients3, i, dTxi, dTyi, dTzi ) ;
	            Coordinates3_IncrementRow ( gradients3, j, dTxj, dTyj, dTzj ) ;
	            Coordinates3_IncrementRow ( gradients3, k, dTxk, dTyk, dTzk ) ;
        	}
	    }
	}
    }
    return energy ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Dihedral i-j-k-l.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real CosineTermEnergy_Dihedral ( const CosineTermContainer *self         ,
                                 const Coordinates3        *coordinates3 ,
                                       Coordinates3        *gradients3   )
{
    Real energy = 0.0e+00 ;
    if ( ( self != NULL ) && ( coordinates3 != NULL ) )
    {
	auto Boolean doGradients ;
        auto Integer i, j, k, l, nt, p, t ;
        auto Real    c, cn, co, cosPhi, dF, dTxi, dTyi, dTzi, dTxj, dTyj, dTzj, dTxk, dTyk, dTzk, dTxl, dTyl, dTzl,
                     e, m, mx, my, mz, n, nx, ny, nz, sx, sy, sz, xij, yij, zij, xkj, ykj, zkj, xlk, ylk, zlk ;
	doGradients = ( gradients3 != NULL ) ;
	for ( nt = 0 ; nt < self->nTerms ; nt++ )
	{
	    if ( self->terms[nt].isActive )
	    {
                /* . Local data. */
	        i = self->terms[nt].indices[0] ;
	        j = self->terms[nt].indices[1] ;
	        k = self->terms[nt].indices[2] ;
	        l = self->terms[nt].indices[3] ;
	        t = self->terms[nt].type ;
	        /* . Coordinate displacements. */
	        Coordinates3_DifferenceRow ( coordinates3, i, j, xij, yij, zij ) ;
	        Coordinates3_DifferenceRow ( coordinates3, k, j, xkj, ykj, zkj ) ;
	        Coordinates3_DifferenceRow ( coordinates3, l, k, xlk, ylk, zlk ) ;
                /* . m and n. */
                mx = yij * zkj - zij * ykj ;
	        my = zij * xkj - xij * zkj ;
	        mz = xij * ykj - yij * xkj ;
	        nx = ylk * zkj - zlk * ykj ;
	        ny = zlk * xkj - xlk * zkj ;
	        nz = xlk * ykj - ylk * xkj ;
                /* . Normalize. */
                m = sqrt ( mx * mx + my * my + mz * mz ) ;
                mx /= m ; my /= m ; mz /= m ;
	        n = sqrt ( nx * nx + ny * ny + nz * nz ) ;
	        nx /= n ; ny /= n ; nz /= n ;
                /* . Cosine of the dihedral. */
                cosPhi = mx * nx +  my * ny +  mz * nz ;
                /* . Loop over powers of the cosine. */
                cn = 1.0e+00 ;
                co = 0.0e+00 ;
                dF = 0.0e+00 ;
                e  = 0.0e+00 ;
                for ( p = 0 ; p <= self->parameters[t].nPowers ; p++ )
                {
                    c   = self->parameters[t].powerCoefficients[p] ;
                    dF += c * co * ( Real ) p ;
                    e  += c * cn ;
                    co  = cn ;
                    cn *= cosPhi ;
                }
                /* . The energy term. */
                energy += e ;
                if ( doGradients )
                {
                    /* . rkj ^ m terms. */
                    sx = ykj * mz - zkj * my ;
	            sy = zkj * mx - xkj * mz ;
	            sz = xkj * my - ykj * mx ;
                    dTxi = - cosPhi * sx ;
                    dTyi = - cosPhi * sy ;
                    dTzi = - cosPhi * sz ;
                    dTxl = sx ;
                    dTyl = sy ;
                    dTzl = sz ;
                    /* . rkj ^ n terms. */
                    sx = ykj * nz - zkj * ny ;
	            sy = zkj * nx - xkj * nz ;
	            sz = xkj * ny - ykj * nx ;
                    dTxi += sx ;
                    dTyi += sy ;
                    dTzi += sz ;
                    dTxl -= cosPhi * sx ;
                    dTyl -= cosPhi * sy ;
                    dTzl -= cosPhi * sz ;
                    /* . Finish i and l. */
                    dTxi *= dF / m ;
                    dTyi *= dF / m ;
                    dTzi *= dF / m ;
                    dTxl *= dF / n ;
                    dTyl *= dF / n ;
                    dTzl *= dF / n ;
                    /* . Scale rij. */
                    xij /= m ; yij /= m ; zij /= m ;
                    /* . rij ^ m terms. */
                    sx = yij * mz - zij * my ;
	            sy = zij * mx - xij * mz ;
	            sz = xij * my - yij * mx ;
                    dTxk = cosPhi * sx ;
                    dTyk = cosPhi * sy ;
                    dTzk = cosPhi * sz ;
                    /* . rij ^ n terms. */
                    sx = yij * nz - zij * ny ;
	            sy = zij * nx - xij * nz ;
	            sz = xij * ny - yij * nx ;
                    dTxk -= sx ;
                    dTyk -= sy ;
                    dTzk -= sz ;
                    /* . Scale rlk. */
                    xlk /= n ; ylk /= n ; zlk /= n ;
                    /* . rlk ^ m terms. */
                    sx = ylk * mz - zlk * my ;
	            sy = zlk * mx - xlk * mz ;
	            sz = xlk * my - ylk * mx ;
                    dTxk -= sx ;
                    dTyk -= sy ;
                    dTzk -= sz ;
                    /* . rlk ^ n terms. */
                    sx = ylk * nz - zlk * ny ;
	            sy = zlk * nx - xlk * nz ;
	            sz = xlk * ny - ylk * nx ;
                    dTxk += cosPhi * sx ;
                    dTyk += cosPhi * sy ;
                    dTzk += cosPhi * sz ;
                    /* . Scale k. */
                    dTxk *= dF ; dTyk *= dF ; dTzk *= dF ;
                    /* . Finish j and k. */
	            dTxj  = - dTxk - dTxi ;
	            dTyj  = - dTyk - dTyi ;
	            dTzj  = - dTzk - dTzi ;
	            dTxk -= dTxl ;
	            dTyk -= dTyl ;
	            dTzk -= dTzl ;
                    /* . Add in the contributions. */
	            Coordinates3_IncrementRow ( gradients3, i, dTxi, dTyi, dTzi ) ;
	            Coordinates3_IncrementRow ( gradients3, j, dTxj, dTyj, dTzj ) ;
	            Coordinates3_IncrementRow ( gradients3, k, dTxk, dTyk, dTzk ) ;
	            Coordinates3_IncrementRow ( gradients3, l, dTxl, dTyl, dTzl ) ;
        	}
	    }
	}
    }
    return energy ;
}


/*----------------------------------------------------------------------------------------------------------------------------------
! . Out-of-plane i-j-(k,l).
! . The angle between ij and the normal to the plane ijk.
! . No distinction is made between positive and negative angles and so the cosine expansion should reflect this.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real CosineTermEnergy_OutOfPlane ( const CosineTermContainer *self         ,
                                   const Coordinates3        *coordinates3 ,
                                         Coordinates3        *gradients3   )
{
    Real energy = 0.0e+00 ;
    if ( ( self != NULL ) && ( coordinates3 != NULL ) )
    {
	auto Boolean doGradients ;
        auto Integer i, j, k, l, nt, p, t ;
        auto Real    c, cn, co, cosPhi, dF, dNx, dNy, dNz, dTxi, dTyi, dTzi, dTxj, dTyj, dTzj, dTxk, dTyk, dTzk, dTxl, dTyl, dTzl,
                     e, n, nx, ny, nz, rij, xij, yij, zij, xkj, ykj, zkj, xlj, ylj, zlj ;
	doGradients = ( gradients3 != NULL ) ;
	for ( nt = 0 ; nt < self->nTerms ; nt++ )
	{
	    if ( self->terms[nt].isActive )
	    {
                /* . Local data. */
	        i = self->terms[nt].indices[0] ;
	        j = self->terms[nt].indices[1] ;
	        k = self->terms[nt].indices[2] ;
	        l = self->terms[nt].indices[3] ;
	        t = self->terms[nt].type ;
	        /* . Coordinate displacements. */
	        Coordinates3_DifferenceRow ( coordinates3, i, j, xij, yij, zij ) ;
	        Coordinates3_DifferenceRow ( coordinates3, k, j, xkj, ykj, zkj ) ;
	        Coordinates3_DifferenceRow ( coordinates3, l, j, xlj, ylj, zlj ) ;
                /* . n. */
	        nx = ykj * zlj - zkj * ylj ;
	        ny = zkj * xlj - xkj * zlj ;
	        nz = xkj * ylj - ykj * xlj ;
                /* . Normalize. */
	        n   = sqrt ( nx * nx + ny * ny + nz * nz ) ;
	        nx /= n ; ny /= n ; nz /= n ;
                rij = sqrt ( xij * xij + yij * yij + zij * zij ) ;
                xij /= rij ; yij /= rij ; zij /= rij ;
                /* . Cosine of the dihedral. */
                cosPhi = nx * xij + ny * yij + nz * zij ;
                /* . Loop over powers of the cosine. */
                cn = 1.0e+00 ;
                co = 0.0e+00 ;
                dF = 0.0e+00 ;
                e  = 0.0e+00 ;
                for ( p = 0 ; p <= self->parameters[t].nPowers ; p++ )
                {
                    c   = self->parameters[t].powerCoefficients[p] ;
                    dF += c * co * ( Real ) p ;
                    e  += c * cn ;
                    co  = cn ;
                    cn *= cosPhi ;
                }
                /* . The energy term. */
                energy += e ;
                if ( doGradients )
                {
                    /* . i terms. */
                    dTxi = dF * ( nx - cosPhi * xij ) / rij ;
                    dTyi = dF * ( ny - cosPhi * yij ) / rij ;
                    dTzi = dF * ( nz - cosPhi * zij ) / rij ;
                    /* . n terms. */
                    dNx  = dF * ( xij - cosPhi * nx ) / n ; 
                    dNy  = dF * ( yij - cosPhi * ny ) / n ; 
                    dNz  = dF * ( zij - cosPhi * nz ) / n ; 
                    /* . k terms. */
                    dTxk = dNz * ylj - dNy * zlj ;
                    dTyk = dNx * zlj - dNz * xlj ;
                    dTzk = dNy * xlj - dNx * ylj ;
                    /* . l terms. */
                    dTxl = dNy * zkj - dNz * ykj ;
                    dTyl = dNz * xkj - dNx * zkj ;
                    dTzl = dNx * ykj - dNy * xkj ;
                    /* . j terms. */
	            dTxj  = - dTxi - dTxk - dTxl ;
	            dTyj  = - dTyi - dTyk - dTyl ;
	            dTzj  = - dTzi - dTzk - dTzl ;
                    /* . Add in the contributions. */
	            Coordinates3_IncrementRow ( gradients3, i, dTxi, dTyi, dTzi ) ;
	            Coordinates3_IncrementRow ( gradients3, j, dTxj, dTyj, dTzj ) ;
	            Coordinates3_IncrementRow ( gradients3, k, dTxk, dTyk, dTzk ) ;
	            Coordinates3_IncrementRow ( gradients3, l, dTxl, dTyl, dTzl ) ;
        	}
	    }
	}
    }
    return energy ;
}
