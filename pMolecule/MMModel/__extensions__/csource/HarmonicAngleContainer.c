/*==============================================================================
!=============================================================================*/

# include <math.h>
# include <stdlib.h>

# include "BooleanUtilities.h"
# include "HarmonicAngleContainer.h"
# include "Memory.h"
# include "NumericalMacros.h"

/* . Use cos(t-t0)/sin(t-t0) in calculation and normalize from c^2+s^2? */

/*------------------------------------------------------------------------------
! . Local procedures.
!-----------------------------------------------------------------------------*/
static Integer HarmonicAngleTerm_Compare ( const void *vTerm1, const void *vTerm2 ) ;

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
/* . The tolerance for calculating angle linearity. */
# define DOT_LIMIT 0.999999

/*==============================================================================
! . Procedures.
!=============================================================================*/
/*------------------------------------------------------------------------------
! . Activate terms.
!-----------------------------------------------------------------------------*/
void HarmonicAngleContainer_ActivateTerms ( HarmonicAngleContainer *self )
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
HarmonicAngleContainer *HarmonicAngleContainer_Allocate ( const Integer nTerms, const Integer nParameters )
{
    HarmonicAngleContainer *self = NULL ;
    if ( ( nTerms != 0 ) && ( nParameters != 0 ) )
    {
        Integer i ;
	self = Memory_AllocateType ( HarmonicAngleContainer ) ;
        self->isSorted     = False       ;
	self->nTerms      = nTerms      ;
	self->nParameters = nParameters ;
	self->terms	  = Memory_AllocateArrayOfTypes ( nTerms     , HarmonicAngle          ) ;
	self->parameters  = Memory_AllocateArrayOfTypes ( nParameters, HarmonicAngleParameter ) ;
	/* . Make all terms inactive. */
	for ( i = 0 ; i < nTerms ; i++ ) self->terms[i].isActive = False ;
    }
    return self ;
}

/*------------------------------------------------------------------------------
! . Cloning.
!-----------------------------------------------------------------------------*/
HarmonicAngleContainer *HarmonicAngleContainer_Clone ( const HarmonicAngleContainer *self )
{
    HarmonicAngleContainer *new = NULL ;
    if ( self != NULL )
    {
        auto Integer i ;
        new = HarmonicAngleContainer_Allocate ( self->nTerms, self->nParameters ) ;
        for ( i = 0 ; i < self->nTerms ; i++ )
        {
            new->terms[i].isActive = self->terms[i].isActive ;
            new->terms[i].atom1   = self->terms[i].atom1   ;
            new->terms[i].atom2   = self->terms[i].atom2   ;
            new->terms[i].atom3   = self->terms[i].atom3   ;
            new->terms[i].type    = self->terms[i].type    ;
        }
        for ( i = 0 ; i < self->nParameters ; i++ )
        {
            new->parameters[i].eq = self->parameters[i].eq ;
            new->parameters[i].fc = self->parameters[i].fc ;
        }
        new->isSorted = self->isSorted ;
    }
    return new ;
}

/*------------------------------------------------------------------------------
! . Deactivate terms that do not involve any atoms in the selection.
! . Already deactivated terms are not affected.
!-----------------------------------------------------------------------------*/
void HarmonicAngleContainer_DeactivateTerms ( HarmonicAngleContainer *self, Selection *selection )
{
    if ( ( self != NULL ) && ( selection != NULL ) )
    {
        auto BooleanBlock *flags ;
        auto Integer       i, n ;
        n     = HarmonicAngleContainer_UpperBound ( self ) ;
        flags = Selection_MakeFlags ( selection, n, NULL ) ;
	for ( i = 0 ; i < self->nTerms ; i++ )
	{
            if ( self->terms[i].isActive )
            {
                self->terms[i].isActive = ( Block_Item ( flags, self->terms[i].atom1 ) ||
                                            Block_Item ( flags, self->terms[i].atom2 ) ||
                                            Block_Item ( flags, self->terms[i].atom3 ) ) ;
            }
	}
    }
}

/*------------------------------------------------------------------------------
! . Deallocation.
!-----------------------------------------------------------------------------*/
void HarmonicAngleContainer_Deallocate ( HarmonicAngleContainer **self )
{
    if ( (*self) != NULL )
    {
        Memory_Deallocate ( (*self)->terms      ) ;
        Memory_Deallocate ( (*self)->parameters ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*------------------------------------------------------------------------------
! . Energy and gradients.
!-----------------------------------------------------------------------------*/
double HarmonicAngleContainer_Energy ( const HarmonicAngleContainer *self, const Coordinates3 *coordinates3, Coordinates3 *gradients3 )
{
    Real energy = 0.0e+00 ;
    if ( ( self != NULL ) && ( coordinates3 != NULL ) )
    {
	auto Boolean   QGRADIENTS ;
        auto Real df, disp, dot, dtdx, dtxi, dtxk, dtyi, dtyk, dtzi, dtzk, rij, rkj, theta, xij, yij, zij, xkj, ykj, zkj ;
        auto Integer    i, j, k, n, t ;
	QGRADIENTS = ( gradients3 != NULL ) ;
	for ( n = 0 ; n < self->nTerms ; n++ )
	{
	    if ( self->terms[n].isActive )
	    {
		i = self->terms[n].atom1 ;
		j = self->terms[n].atom2 ;
	        k = self->terms[n].atom3 ;
		t = self->terms[n].type  ;
	        Coordinates3_DifferenceRow ( coordinates3, i, j, xij, yij, zij ) ;
	        Coordinates3_DifferenceRow ( coordinates3, k, j, xkj, ykj, zkj ) ;
	        rij = sqrt ( xij * xij + yij * yij + zij * zij ) ;
	        rkj = sqrt ( xkj * xkj + ykj * ykj + zkj * zkj ) ;
	        xij /= rij ; yij /= rij ; zij /= rij ;
	        xkj /= rkj ; ykj /= rkj ; zkj /= rkj ;
	        dot   = xij * xkj + yij * ykj + zij * zkj ;
                dot   = Maximum ( - DOT_LIMIT, dot ) ;
	        dot   = Minimum (   DOT_LIMIT, dot ) ;
                theta = acos ( dot ) ;
	        disp  = theta - self->parameters[t].eq ;
	        df    = self->parameters[t].fc * disp  ;
	        energy += ( df * disp ) ;
                if ( QGRADIENTS )
                {
	            dtdx = - 1.0 / sqrt ( 1.0 - dot * dot ) ;
                    df  *= ( 2.0e+00 * dtdx ) ;
                    dtxi = df * ( xkj - dot * xij ) / rij ;
                    dtyi = df * ( ykj - dot * yij ) / rij ;
                    dtzi = df * ( zkj - dot * zij ) / rij ;
                    dtxk = df * ( xij - dot * xkj ) / rkj ;
                    dtyk = df * ( yij - dot * ykj ) / rkj ;
                    dtzk = df * ( zij - dot * zkj ) / rkj ;
	            Coordinates3_IncrementRow ( gradients3, i,   dtxi,            dtyi,            dtzi          ) ;
	            Coordinates3_IncrementRow ( gradients3, k,          dtxk,            dtyk,            dtzk   ) ;
	            Coordinates3_DecrementRow ( gradients3, j, ( dtxi + dtxk ), ( dtyi + dtyk ), ( dtzi + dtzk ) ) ;
        	}
	    }
	}
    }
    return energy ;
}

/*------------------------------------------------------------------------------
! . Merging.
!-----------------------------------------------------------------------------*/
HarmonicAngleContainer *HarmonicAngleContainer_Merge ( const HarmonicAngleContainer *self, const HarmonicAngleContainer *other, const Integer atomincrement )
{
    HarmonicAngleContainer *new = NULL ;
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Integer i ;
        new = HarmonicAngleContainer_Allocate ( self->nTerms + other->nTerms, self->nParameters + other->nParameters ) ;
        for ( i = 0 ; i < self->nTerms ; i++ )
        {
            new->terms[i].isActive = self->terms[i].isActive ;
            new->terms[i].atom1   = self->terms[i].atom1   ;
            new->terms[i].atom2   = self->terms[i].atom2   ;
            new->terms[i].atom3   = self->terms[i].atom3   ;
            new->terms[i].type    = self->terms[i].type    ;
        }
        for ( i = 0 ; i < other->nTerms ; i++ )
        {
            new->terms[i+self->nTerms].isActive = other->terms[i].isActive ;
            new->terms[i+self->nTerms].atom1   = other->terms[i].atom1 + atomincrement     ;
            new->terms[i+self->nTerms].atom2   = other->terms[i].atom2 + atomincrement     ;
            new->terms[i+self->nTerms].atom3   = other->terms[i].atom3 + atomincrement     ;
            new->terms[i+self->nTerms].type    = other->terms[i].type  + self->nParameters ;
        }
        for ( i = 0 ; i < self->nParameters ; i++ )
        {
            new->parameters[i].eq = self->parameters[i].eq ;
            new->parameters[i].fc = self->parameters[i].fc ;
        }
        for ( i = 0 ; i < other->nParameters ; i++ )
        {
            new->parameters[i+self->nParameters].eq = other->parameters[i].eq ;
            new->parameters[i+self->nParameters].fc = other->parameters[i].fc ;
        }
        new->isSorted = ( self->isSorted && other->isSorted ) ;
    }
    return new ;
}

/*------------------------------------------------------------------------------
! . Return the number of inactive terms.
!-----------------------------------------------------------------------------*/
int HarmonicAngleContainer_NumberOfInactiveTerms ( const HarmonicAngleContainer *self )
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
HarmonicAngleContainer *HarmonicAngleContainer_Prune ( HarmonicAngleContainer *self, Selection *selection )
{
    HarmonicAngleContainer *new = NULL ;
    if ( ( self != NULL ) && ( selection != NULL ) )
    {
        auto Boolean       *toKeep    ;
        auto BooleanBlock  *flags     ;
        auto IntegerBlock  *positions ;
        auto Integer        i, n      ;
        n         = HarmonicAngleContainer_UpperBound ( self ) ;
        flags     = Selection_MakeFlags     ( selection, n, NULL ) ;
        positions = Selection_MakePositions ( selection, n, NULL ) ;
        toKeep    = Boolean_Allocate        ( self->nTerms, NULL ) ;
	for ( i = n = 0 ; i < self->nTerms ; i++ )
	{
            toKeep[i] = ( Block_Item ( flags, self->terms[i].atom1 ) &&
                          Block_Item ( flags, self->terms[i].atom2 ) &&
                          Block_Item ( flags, self->terms[i].atom3 ) ) ;
            if ( toKeep[i] ) n++ ;
	}
	if ( n > 0 )
	{
            new = HarmonicAngleContainer_Allocate ( n, self->nParameters ) ;
            for ( i = 0 ; i < self->nParameters ; i++ )
            {
                new->parameters[i].eq = self->parameters[i].eq ;
                new->parameters[i].fc = self->parameters[i].fc ;
            }
            for ( i = n = 0 ; i < self->nTerms ; i++ )
            {
                if ( toKeep[i] )
                {
        	    new->terms[n].isActive = self->terms[i].isActive ;
        	    new->terms[n].atom1   = Block_Item ( positions, self->terms[i].atom1 ) ;
        	    new->terms[n].atom2   = Block_Item ( positions, self->terms[i].atom2 ) ;
        	    new->terms[n].atom3   = Block_Item ( positions, self->terms[i].atom3 ) ;
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
! . Within a harmonicangle, atom1 > atom3.
! . Within the array, ordering is done with increased values of atom2 and then
! . atom1 and then atom3.
! . Duplicate terms are not removed.
!-----------------------------------------------------------------------------*/
void HarmonicAngleContainer_Sort ( HarmonicAngleContainer *self )
{
    if ( ( self != NULL ) && ( ! self->isSorted ) )
    {
        auto Integer atom1, atom3, i ;
        /* . Order atom1 and atom2 within each term. */
        for ( i = 0 ; i < self->nTerms ; i++ )
        {
            atom1 = self->terms[i].atom1 ;
            atom3 = self->terms[i].atom3 ;
            if ( atom3 > atom1 )
            {
                self->terms[i].atom1 = atom3 ;
                self->terms[i].atom3 = atom1 ;
            }
        }
        /* . Order the terms within the container. */
        qsort ( ( void * ) self->terms, ( size_t ) self->nTerms, sizeof ( HarmonicAngle ), ( void * ) HarmonicAngleTerm_Compare ) ;
        self->isSorted = True ;
    }
}

/*------------------------------------------------------------------------------
! . Return the upper bound for the container.
! . This is the value of the largest index plus one.
!-----------------------------------------------------------------------------*/
int HarmonicAngleContainer_UpperBound ( HarmonicAngleContainer *self )
{
    Integer upperBound = 0 ;
    if ( ( self != NULL ) && ( self->nTerms > 0 ) )
    {
        HarmonicAngleContainer_Sort ( self ) ;
        upperBound = Maximum ( self->terms[self->nTerms-1].atom1, self->terms[self->nTerms-1].atom2 ) + 1 ;
    }
    return upperBound ;
}

/*==============================================================================
! . Private procedures.
!============================================================================*/
static Integer HarmonicAngleTerm_Compare ( const void *vTerm1, const void *vTerm2 )
{
    HarmonicAngle *term1, *term2 ;
    Integer i ;
    term1 = ( HarmonicAngle * ) vTerm1 ;
    term2 = ( HarmonicAngle * ) vTerm2 ;
         if ( term1->atom2 < term2->atom2 ) i = -1 ;
    else if ( term1->atom2 > term2->atom2 ) i =  1 ;
    else if ( term1->atom1 < term2->atom1 ) i = -1 ;
    else if ( term1->atom1 > term2->atom1 ) i =  1 ;
    else if ( term1->atom3 < term2->atom3 ) i = -1 ;
    else if ( term1->atom3 > term2->atom3 ) i =  1 ;
    else if ( term1->type  < term2->type  ) i = -1 ;
    else if ( term1->type  > term2->type  ) i =  1 ;
    else if ( ! term1->isActive && term2->isActive ) i = -1 ;
    else if ( term1->isActive && ! term2->isActive ) i =  1 ;
    else i = 0 ;
    return i ;
}
