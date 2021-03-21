/*==============================================================================
!=============================================================================*/

# include <math.h>
# include <stdlib.h>

# include "BooleanUtilities.h"
# include "HarmonicImproperContainer.h"
# include "Memory.h"
# include "NumericalMacros.h"

/*------------------------------------------------------------------------------
! . Local procedures.
!-----------------------------------------------------------------------------*/
static Integer HarmonicImproperTerm_Compare ( const void *vTerm1, const void *vTerm2 ) ;

/*==============================================================================
! . Procedures.
!=============================================================================*/
/*------------------------------------------------------------------------------
! . Activate terms.
!-----------------------------------------------------------------------------*/
void HarmonicImproperContainer_ActivateTerms ( HarmonicImproperContainer *self )
{
    if ( self != NULL )
    {
        auto Integer i ;
	for ( i = 0 ; i < self->nTerms ; i++ ) self->terms[i].isActive = True ;
    }
}

/*------------------------------------------------------------------------------
! . Allocation.
!-----------------------------------------------------------------------------*/
HarmonicImproperContainer *HarmonicImproperContainer_Allocate ( const Integer nTerms, const Integer nParameters )
{
    HarmonicImproperContainer *self = NULL ;
    if ( ( nTerms != 0 ) && ( nParameters != 0 ) )
    {
        Integer i ;
	self = Memory_AllocateType ( HarmonicImproperContainer ) ;
        self->isSorted     = False       ;
	self->nTerms      = nTerms      ;
	self->nParameters = nParameters ;
	self->terms	  = Memory_AllocateArrayOfTypes ( nTerms     , HarmonicImproper          ) ;
	self->parameters  = Memory_AllocateArrayOfTypes ( nParameters, HarmonicImproperParameter ) ;
	/* . Make all terms inactive. */
	for ( i = 0 ; i < nTerms ; i++ ) self->terms[i].isActive = False ;
    }
    return self ;
}

/*------------------------------------------------------------------------------
! . Cloning.
!-----------------------------------------------------------------------------*/
HarmonicImproperContainer *HarmonicImproperContainer_Clone ( const HarmonicImproperContainer *self )
{
    HarmonicImproperContainer *new = NULL ;
    if ( self != NULL )
    {
        auto Integer i ;
        new = HarmonicImproperContainer_Allocate ( self->nTerms, self->nParameters ) ;
        for ( i = 0 ; i < self->nTerms ; i++ )
        {
            new->terms[i].isActive = self->terms[i].isActive ;
            new->terms[i].atom1   = self->terms[i].atom1   ;
            new->terms[i].atom2   = self->terms[i].atom2   ;
            new->terms[i].atom3   = self->terms[i].atom3   ;
            new->terms[i].atom4   = self->terms[i].atom4   ;
            new->terms[i].type    = self->terms[i].type    ;
        }
        for ( i = 0 ; i < self->nParameters ; i++ )
        {
            new->parameters[i].eq    = self->parameters[i].eq    ;
            new->parameters[i].fc    = self->parameters[i].fc    ;
            new->parameters[i].coseq = self->parameters[i].coseq ;
            new->parameters[i].sineq = self->parameters[i].sineq ;
        }
        new->isSorted = self->isSorted ;
    }
    return new ;
}

/*------------------------------------------------------------------------------
! . Deactivate terms that do not involve any atoms in the selection.
! . Already deactivated terms are not affected.
!-----------------------------------------------------------------------------*/
void HarmonicImproperContainer_DeactivateTerms ( HarmonicImproperContainer *self, Selection *selection )
{
    if ( ( self != NULL ) && ( selection != NULL ) )
    {
        auto BooleanBlock *flags ;
        auto Integer       i, n ;
        n     = HarmonicImproperContainer_UpperBound ( self ) ;
        flags = Selection_MakeFlags ( selection, n, NULL ) ;
	for ( i = 0 ; i < self->nTerms ; i++ )
	{
            if ( self->terms[i].isActive )
            {
                self->terms[i].isActive = ( Block_Item ( flags, self->terms[i].atom1 ) ||
                                            Block_Item ( flags, self->terms[i].atom2 ) ||
                                            Block_Item ( flags, self->terms[i].atom3 ) ||
                                            Block_Item ( flags, self->terms[i].atom4 ) ) ;
            }
	}
    }
}

/*------------------------------------------------------------------------------
! . Deallocation.
!-----------------------------------------------------------------------------*/
void HarmonicImproperContainer_Deallocate ( HarmonicImproperContainer **self )
{
    if ( (*self) != NULL )
    {
        Memory_Deallocate ( (*self)->terms      ) ;
        Memory_Deallocate ( (*self)->parameters ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*------------------------------------------------------------------------------
! . Fill the cos and sin equilibrium values of the parameter array.
!-----------------------------------------------------------------------------*/
void HarmonicImproperContainer_FillCosSinValues ( HarmonicImproperContainer *self )
{
    if ( self != NULL )
    {
        auto Integer i ;
        for ( i = 0; i < self->nParameters ; i++ )
        {
            self->parameters[i].coseq = cos ( self->parameters[i].eq ) ;
            self->parameters[i].sineq = sin ( self->parameters[i].eq ) ;
        }
    }
}

/*------------------------------------------------------------------------------
! . Energy and gradients.
! . Following Becker, Berendsen and van Gunsteren, JCC 16 p527 (1995).
!-----------------------------------------------------------------------------*/
# define LOWCOSPHI 0.1e+00
double HarmonicImproperContainer_Energy ( const HarmonicImproperContainer *self, const Coordinates3 *coordinates3, Coordinates3 *gradients3 )
{
    Real energy = 0.0e+00 ;
    if ( ( self != NULL ) && ( coordinates3 != NULL ) )
    {
	auto Boolean   QGRADIENTS ;
        auto Real cosdphi, cosPhi, df, dotij, dotlk, dphi, mn, rkj, rkj2, sindphi, sinPhi ;
        auto Real dtxi, dtyi, dtzi, dtxj, dtyj, dtzj, dtxk, dtyk, dtzk, dtxl, dtyl, dtzl, m2, mx, my, mz, n2, nx, ny, nz, sx, sy, sz,
                    xij, yij, zij, xkj, ykj, zkj, xlk, ylk, zlk ;
        auto Integer    i, j, k, l, n, t ;
	QGRADIENTS = ( gradients3 != NULL ) ;
	for ( n = 0 ; n < self->nTerms ; n++ )
	{
	    if ( self->terms[n].isActive )
	    {
                /* . Local data. */
	        i = self->terms[n].atom1 ;
	        j = self->terms[n].atom2 ;
	        k = self->terms[n].atom3 ;
	        l = self->terms[n].atom4 ;
	        t = self->terms[n].type  ;
	        /* . Coordinate displacements. */
	        Coordinates3_DifferenceRow ( coordinates3, i, j, xij, yij, zij ) ;
	        Coordinates3_DifferenceRow ( coordinates3, k, j, xkj, ykj, zkj ) ;
	        Coordinates3_DifferenceRow ( coordinates3, l, k, xlk, ylk, zlk ) ;
	        rkj2 = xkj * xkj + ykj * ykj + zkj * zkj ;
	        rkj  = sqrt ( rkj2 ) ;
                /* . m and n. */
                mx = yij * zkj - zij * ykj ;
	        my = zij * xkj - xij * zkj ;
	        mz = xij * ykj - yij * xkj ;
	        nx = ylk * zkj - zlk * ykj ;
	        ny = zlk * xkj - xlk * zkj ;
	        nz = xlk * ykj - ylk * xkj ;
                m2 = mx * mx + my * my + mz * mz ;
	        n2 = nx * nx + ny * ny + nz * nz ;
	        mn = sqrt ( m2 * n2 ) ;
                /* . Cosine and sine of the dihedral. */
                cosPhi =       (  mx * nx +  my * ny +  mz * nz ) / mn ;
                sinPhi = rkj * ( xij * nx + yij * ny + zij * nz ) / mn ;
                /* . Cosine and sine of (phi - phi0). */
                cosdphi = cosPhi * self->parameters[t].coseq + sinPhi * self->parameters[t].sineq ;
                sindphi = sinPhi * self->parameters[t].coseq - cosPhi * self->parameters[t].sineq ;
                /* . Follow CHARMM here. */
                if ( cosdphi > LOWCOSPHI ) dphi = asin ( sindphi ) ;
                else
                {
                    dphi = fabs ( acos ( Maximum ( cosdphi, -1.0e+00 ) ) ) ;
                    if ( sindphi < 0.0e+00 ) dphi *= -1.0e+00 ;
                }
                /* . The energy term. */
                df      = self->parameters[t].fc * dphi ;
                energy += df * dphi ;
                if ( QGRADIENTS )
                {
	            /* . The derivatives. */
	            df *= 2.0e+00 ;
                    /* . i and l. */
                    dtxi =   df * rkj * mx / m2 ;
                    dtyi =   df * rkj * my / m2 ;
                    dtzi =   df * rkj * mz / m2 ;
                    dtxl = - df * rkj * nx / n2 ;
                    dtyl = - df * rkj * ny / n2 ;
                    dtzl = - df * rkj * nz / n2 ;
                    /* . j and k. */
                    dotij = xij * xkj + yij * ykj + zij * zkj ;
                    dotlk = xlk * xkj + ylk * ykj + zlk * zkj ;
	            sx    = ( dotij * dtxi + dotlk * dtxl ) / rkj2 ;
	            sy    = ( dotij * dtyi + dotlk * dtyl ) / rkj2 ;
	            sz    = ( dotij * dtzi + dotlk * dtzl ) / rkj2 ;
	            dtxj  =   sx - dtxi ;
	            dtyj  =   sy - dtyi ;
	            dtzj  =   sz - dtzi ;
	            dtxk  = - sx - dtxl ;
	            dtyk  = - sy - dtyl ;
	            dtzk  = - sz - dtzl ;
                    /* . Add in the contributions. */
	            Coordinates3_IncrementRow ( gradients3, i, dtxi, dtyi, dtzi ) ;
	            Coordinates3_IncrementRow ( gradients3, j, dtxj, dtyj, dtzj ) ;
	            Coordinates3_IncrementRow ( gradients3, k, dtxk, dtyk, dtzk ) ;
	            Coordinates3_IncrementRow ( gradients3, l, dtxl, dtyl, dtzl ) ;
        	}
	    }
	}
    }
    return energy ;
}

/*------------------------------------------------------------------------------
! . Merging.
!-----------------------------------------------------------------------------*/
HarmonicImproperContainer *HarmonicImproperContainer_Merge ( const HarmonicImproperContainer *self, const HarmonicImproperContainer *other, const Integer atomincrement )
{
    HarmonicImproperContainer *new = NULL ;
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Integer i ;
        new = HarmonicImproperContainer_Allocate ( self->nTerms + other->nTerms, self->nParameters + other->nParameters ) ;
        for ( i = 0 ; i < self->nTerms ; i++ )
        {
            new->terms[i].isActive = self->terms[i].isActive ;
            new->terms[i].atom1   = self->terms[i].atom1   ;
            new->terms[i].atom2   = self->terms[i].atom2   ;
            new->terms[i].atom3   = self->terms[i].atom3   ;
            new->terms[i].atom4   = self->terms[i].atom4   ;
            new->terms[i].type    = self->terms[i].type    ;
        }
        for ( i = 0 ; i < other->nTerms ; i++ )
        {
            new->terms[i+self->nTerms].isActive = other->terms[i].isActive ;
            new->terms[i+self->nTerms].atom1   = other->terms[i].atom1 + atomincrement     ;
            new->terms[i+self->nTerms].atom2   = other->terms[i].atom2 + atomincrement     ;
            new->terms[i+self->nTerms].atom3   = other->terms[i].atom3 + atomincrement     ;
            new->terms[i+self->nTerms].atom4   = other->terms[i].atom4 + atomincrement     ;
            new->terms[i+self->nTerms].type    = other->terms[i].type  + self->nParameters ;
        }
        for ( i = 0 ; i < self->nParameters ; i++ )
        {
            new->parameters[i].eq    = self->parameters[i].eq    ;
            new->parameters[i].fc    = self->parameters[i].fc    ;
            new->parameters[i].coseq = self->parameters[i].coseq ;
            new->parameters[i].sineq = self->parameters[i].sineq ;
        }
        for ( i = 0 ; i < other->nParameters ; i++ )
        {
            new->parameters[i+self->nParameters].eq    = other->parameters[i].eq    ;
            new->parameters[i+self->nParameters].fc    = other->parameters[i].fc    ;
            new->parameters[i+self->nParameters].coseq = other->parameters[i].coseq ;
            new->parameters[i+self->nParameters].sineq = other->parameters[i].sineq ;
        }
        new->isSorted = ( self->isSorted && other->isSorted ) ;
    }
    return new ;
}

/*------------------------------------------------------------------------------
! . Return the number of inactive terms.
!-----------------------------------------------------------------------------*/
int HarmonicImproperContainer_NumberOfInactiveTerms ( const HarmonicImproperContainer *self )
{
    Integer n = 0 ;
    if ( self != NULL )
    {
        auto Integer i ;
	for ( i = 0 ; i < self->nTerms ; i++ ) if ( ! self->terms[i].isActive ) n++ ;
    }
    return n ;
}

/*------------------------------------------------------------------------------
! . Pruning.
! . Only terms are pruned.
!-----------------------------------------------------------------------------*/
HarmonicImproperContainer *HarmonicImproperContainer_Prune ( HarmonicImproperContainer *self, Selection *selection )
{
    HarmonicImproperContainer *new = NULL ;
    if ( ( self != NULL ) && ( selection != NULL ) )
    {
        auto Boolean       *toKeep    ;
        auto BooleanBlock  *flags     ;
        auto IntegerBlock  *positions ;
        auto Integer        i, n      ;
        n         = HarmonicImproperContainer_UpperBound ( self ) ;
        flags     = Selection_MakeFlags     ( selection, n, NULL ) ;
        positions = Selection_MakePositions ( selection, n, NULL ) ;
        toKeep    = Boolean_Allocate        ( self->nTerms, NULL ) ;
	for ( i = n = 0 ; i < self->nTerms ; i++ )
	{
            toKeep[i] = ( Block_Item ( flags, self->terms[i].atom1 ) &&
                          Block_Item ( flags, self->terms[i].atom2 ) &&
                          Block_Item ( flags, self->terms[i].atom3 ) &&
                          Block_Item ( flags, self->terms[i].atom4 ) ) ;
            if ( toKeep[i] ) n++ ;
	}
	if ( n > 0 )
	{
            new = HarmonicImproperContainer_Allocate ( n, self->nParameters ) ;
            for ( i = 0 ; i < self->nParameters ; i++ )
            {
                new->parameters[i].eq    = self->parameters[i].eq    ;
                new->parameters[i].fc    = self->parameters[i].fc    ;
                new->parameters[i].coseq = self->parameters[i].coseq ;
                new->parameters[i].sineq = self->parameters[i].sineq ;
            }
            for ( i = n = 0 ; i < self->nTerms ; i++ )
            {
                if ( toKeep[i] )
                {
        	    new->terms[n].isActive = self->terms[i].isActive ;
        	    new->terms[n].atom1   = Block_Item ( positions, self->terms[i].atom1 ) ;
        	    new->terms[n].atom2   = Block_Item ( positions, self->terms[i].atom2 ) ;
        	    new->terms[n].atom3   = Block_Item ( positions, self->terms[i].atom3 ) ;
        	    new->terms[n].atom4   = Block_Item ( positions, self->terms[i].atom4 ) ;
        	    new->terms[n].type    = self->terms[i].type ;
        	    n++ ;
                }
            }
            new->isSorted = self->isSorted ;
	}
	Boolean_Deallocate ( &toKeep ) ;
    }
    return new ;
}

/*------------------------------------------------------------------------------
! . Sorting.
! . Within a harmonic improper, atom2 > atom3.
! . Within the array, ordering is done with increased values of atom2 and then
! . atom3 and then atom1 and then atom4.
! . Duplicates are not removed.
!-----------------------------------------------------------------------------*/
void HarmonicImproperContainer_Sort ( HarmonicImproperContainer *self )
{
    if ( ( self != NULL ) && ( ! self->isSorted ) )
    {
        auto Integer atom1, atom2, atom3, atom4, i ;
        /* . Order atom2 and atom3 within each term. */
        for ( i = 0 ; i < self->nTerms ; i++ )
        {
            atom1 = self->terms[i].atom1 ;
            atom2 = self->terms[i].atom2 ;
            atom3 = self->terms[i].atom3 ;
            atom4 = self->terms[i].atom4 ;
            if ( atom3 > atom2 )
            {
                self->terms[i].atom1 = atom4 ;
                self->terms[i].atom2 = atom3 ;
                self->terms[i].atom3 = atom2 ;
                self->terms[i].atom4 = atom1 ;
            }
        }
        /* . Order the terms within the container. */
        qsort ( ( void * ) self->terms, ( size_t ) self->nTerms, sizeof ( HarmonicImproper ), ( void * ) HarmonicImproperTerm_Compare ) ;
        self->isSorted    = True ;
    }
}

/*------------------------------------------------------------------------------
! . Return the upper bound for the container.
! . This is the value of the largest index plus one.
!-----------------------------------------------------------------------------*/
int HarmonicImproperContainer_UpperBound ( HarmonicImproperContainer *self )
{
    Integer upperBound = 0 ;
    if ( ( self != NULL ) && ( self->nTerms > 0 ) )
    {
        HarmonicImproperContainer_Sort ( self ) ;
        upperBound  = Maximum ( self->terms[self->nTerms-1].atom1, self->terms[self->nTerms-1].atom2 ) ;
        upperBound  = Maximum ( upperBound, self->terms[self->nTerms-1].atom4 ) ;
        upperBound += 1 ;
    }
    return upperBound ;
}

/*==============================================================================
! . Private procedures.
!============================================================================*/
static Integer HarmonicImproperTerm_Compare ( const void *vTerm1, const void *vTerm2 )
{
    HarmonicImproper *term1, *term2 ;
    Integer i ;
    term1 = ( HarmonicImproper * ) vTerm1 ;
    term2 = ( HarmonicImproper * ) vTerm2 ;
         if ( term1->atom2 < term2->atom2 ) i = -1 ;
    else if ( term1->atom2 > term2->atom2 ) i =  1 ;
    else if ( term1->atom3 < term2->atom3 ) i = -1 ;
    else if ( term1->atom3 > term2->atom3 ) i =  1 ;
    else if ( term1->atom1 < term2->atom1 ) i = -1 ;
    else if ( term1->atom1 > term2->atom1 ) i =  1 ;
    else if ( term1->atom4 < term2->atom4 ) i = -1 ;
    else if ( term1->atom4 > term2->atom4 ) i =  1 ;
    else if ( term1->type  < term2->type  ) i = -1 ;
    else if ( term1->type  > term2->type  ) i =  1 ;
    else if ( ! term1->isActive && term2->isActive ) i = -1 ;
    else if ( term1->isActive && ! term2->isActive ) i =  1 ;
    else i = 0 ;
    return i ;
}
