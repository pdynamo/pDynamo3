/*==============================================================================
!=============================================================================*/

# include <math.h>
# include <stdlib.h>

# include "BooleanUtilities.h"
# include "HarmonicBondContainer.h"
# include "Memory.h"
# include "NumericalMacros.h"

/*------------------------------------------------------------------------------
! . Local procedures.
!-----------------------------------------------------------------------------*/
static Integer HarmonicBondTerm_Compare ( const void *vTerm1, const void *vTerm2 ) ;

/*==============================================================================
! . Standard procedures.
!=============================================================================*/
/*------------------------------------------------------------------------------
! . Activate terms.
!-----------------------------------------------------------------------------*/
void HarmonicBondContainer_ActivateTerms ( HarmonicBondContainer *self )
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
HarmonicBondContainer *HarmonicBondContainer_Allocate ( const Integer nTerms, const Integer nParameters )
{
    HarmonicBondContainer *self = NULL ;
    if ( ( nTerms != 0 ) && ( nParameters != 0 ) )
    {
        Integer i ;
	self = Memory_AllocateType ( HarmonicBondContainer ) ;
        self->isSorted     = False       ;
	self->nTerms      = nTerms      ;
	self->nParameters = nParameters ;
	self->terms	  = Memory_AllocateArrayOfTypes ( nTerms     , HarmonicBond          ) ;
	self->parameters  = Memory_AllocateArrayOfTypes ( nParameters, HarmonicBondParameter ) ;
	/* . Make all terms inactive. */
	for ( i = 0 ; i < nTerms ; i++ ) self->terms[i].isActive = False ;
    }
    return self ;
}

/*------------------------------------------------------------------------------
! . Cloning.
!-----------------------------------------------------------------------------*/
HarmonicBondContainer *HarmonicBondContainer_Clone ( const HarmonicBondContainer *self )
{
    HarmonicBondContainer *new = NULL ;
    if ( self != NULL )
    {
        auto Integer i ;
        new = HarmonicBondContainer_Allocate ( self->nTerms, self->nParameters ) ;
        for ( i = 0 ; i < self->nTerms ; i++ )
        {
            new->terms[i].isActive = self->terms[i].isActive ;
            new->terms[i].atom1   = self->terms[i].atom1   ;
            new->terms[i].atom2   = self->terms[i].atom2   ;
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
void HarmonicBondContainer_DeactivateTerms ( HarmonicBondContainer *self, Selection *selection )
{
    if ( ( self != NULL ) && ( selection != NULL ) )
    {
        auto BooleanBlock *flags ;
        auto Integer       i, n ;
        n     = HarmonicBondContainer_UpperBound ( self ) ;
        flags = Selection_MakeFlags ( selection, n, NULL ) ;
	for ( i = 0 ; i < self->nTerms ; i++ )
	{
            if ( self->terms[i].isActive )
            {
                self->terms[i].isActive = ( Block_Item ( flags, self->terms[i].atom1 ) || Block_Item ( flags, self->terms[i].atom2 ) ) ;
            }
	}
    }
}

/*------------------------------------------------------------------------------
! . Deallocation.
!-----------------------------------------------------------------------------*/
void HarmonicBondContainer_Deallocate ( HarmonicBondContainer **self )
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
double HarmonicBondContainer_Energy ( const HarmonicBondContainer *self, const Coordinates3 *coordinates3, Coordinates3 *gradients3 )
{
    Real energy = 0.0e+00 ;
    if ( ( self != NULL ) && ( coordinates3 != NULL ) )
    {
	auto Boolean   QGRADIENTS ;
	auto Real df, disp, rij, xij, yij, zij ;
	auto Integer    i, j, n, t ;
	QGRADIENTS = ( gradients3 != NULL ) ;
	for ( n = 0 ; n < self->nTerms ; n++ )
	{
	    if ( self->terms[n].isActive )
	    {
		i = self->terms[n].atom1 ;
		j = self->terms[n].atom2 ;
		t = self->terms[n].type  ;
		Coordinates3_DifferenceRow ( coordinates3, i, j, xij, yij, zij ) ;
		rij  = sqrt ( xij * xij + yij * yij + zij * zij ) ;
		disp = rij - self->parameters[t].eq ;
		df   = self->parameters[t].fc * disp ;
		energy += ( df * disp ) ;
        	if ( QGRADIENTS )
        	{
        	    df  *= ( 2.0e+00 / rij ) ;
		    xij *= df ;
		    yij *= df ;
		    zij *= df ;
		    Coordinates3_IncrementRow ( gradients3, i, xij, yij, zij ) ;
		    Coordinates3_DecrementRow ( gradients3, j, xij, yij, zij ) ;
        	}
	    }
	}
    }
    return energy ;
}

/*------------------------------------------------------------------------------
! . Identify boundary atoms.
! . qcAtoms is the selection of pure QC atoms only.
! . No checking/sorting of the returned data need be done.
!-----------------------------------------------------------------------------*/
int HarmonicBondContainer_IdentifyBoundaryAtoms ( HarmonicBondContainer *self, Selection *qcAtoms, Integer **mmBoundary, Integer **qcPartners )
{
    Integer n = 0 ;
    if ( ( self != NULL ) && ( qcAtoms != NULL ) )
    {
        auto BooleanBlock *flags ;
        auto Integer       i ;
        n     = HarmonicBondContainer_UpperBound ( self ) ;
        flags = Selection_MakeFlags ( qcAtoms, n, NULL ) ;
        /* . Loop to find the number of boundary atoms. */
	for ( i = 0 ; i < self->nTerms ; i++ )
	{
            if ( (   Block_Item ( flags, self->terms[i].atom1 ) && ! Block_Item ( flags, self->terms[i].atom2 ) ) ||
                 ( ! Block_Item ( flags, self->terms[i].atom1 ) &&   Block_Item ( flags, self->terms[i].atom2 ) ) ) n++ ;
	}
        /* . Fill the boundary atom arrays. */
        if ( ( n > 0 ) && ( mmBoundary != NULL ) && ( qcPartners != NULL ) )
        {
            (*mmBoundary) = Memory_AllocateArrayOfTypes ( n, Integer ) ;
            (*qcPartners) = Memory_AllocateArrayOfTypes ( n, Integer ) ;
            for ( i = n = 0 ; i < self->nTerms ; i++ )
	    {
                if ( Block_Item ( flags, self->terms[i].atom1 ) && ! Block_Item ( flags, self->terms[i].atom2 ) )
                {
                    (*qcPartners)[n] = self->terms[i].atom1 ;
                    (*mmBoundary)[n] = self->terms[i].atom2 ;
                    n++ ;
                }
                else if ( ! Block_Item ( flags, self->terms[i].atom1 ) && Block_Item ( flags, self->terms[i].atom2 ) )
                {
                    (*mmBoundary)[n] = self->terms[i].atom1 ;
                    (*qcPartners)[n] = self->terms[i].atom2 ;
                    n++ ;
                }
	    }
        }
    }
    return n ;
}

/*------------------------------------------------------------------------------
! . Merging.
!-----------------------------------------------------------------------------*/
HarmonicBondContainer *HarmonicBondContainer_Merge ( const HarmonicBondContainer *self, const HarmonicBondContainer *other, const Integer atomincrement )
{
    HarmonicBondContainer *new = NULL ;
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Integer i ;
        new = HarmonicBondContainer_Allocate ( self->nTerms + other->nTerms, self->nParameters + other->nParameters ) ;
        for ( i = 0 ; i < self->nTerms ; i++ )
        {
            new->terms[i].isActive = self->terms[i].isActive ;
            new->terms[i].atom1   = self->terms[i].atom1   ;
            new->terms[i].atom2   = self->terms[i].atom2   ;
            new->terms[i].type    = self->terms[i].type    ;
        }
        for ( i = 0 ; i < other->nTerms ; i++ )
        {
            new->terms[i+self->nTerms].isActive = other->terms[i].isActive ;
            new->terms[i+self->nTerms].atom1   = other->terms[i].atom1 + atomincrement     ;
            new->terms[i+self->nTerms].atom2   = other->terms[i].atom2 + atomincrement     ;
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
int HarmonicBondContainer_NumberOfInactiveTerms ( const HarmonicBondContainer *self )
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
HarmonicBondContainer *HarmonicBondContainer_Prune ( HarmonicBondContainer *self, Selection *selection )
{
    HarmonicBondContainer *new = NULL ;
    if ( ( self != NULL ) && ( selection != NULL ) )
    {
        auto Boolean       *toKeep    ;
        auto BooleanBlock  *flags     ;
        auto IntegerBlock  *positions ;
        auto Integer        i, n      ;
        n         = HarmonicBondContainer_UpperBound ( self ) ;
        flags     = Selection_MakeFlags     ( selection, n, NULL ) ;
        positions = Selection_MakePositions ( selection, n, NULL ) ;
        toKeep    = Boolean_Allocate        ( self->nTerms, NULL ) ;
	for ( i = n = 0 ; i < self->nTerms ; i++ )
	{
            toKeep[i] = ( Block_Item ( flags, self->terms[i].atom1 ) &&
                          Block_Item ( flags, self->terms[i].atom2 ) ) ;
            if ( toKeep[i] ) n++ ;
	}
	if ( n > 0 )
	{
            new = HarmonicBondContainer_Allocate ( n, self->nParameters ) ;
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
! . Within a harmonicbond, atom1 > atom3.
! . Within the array, ordering is done with increased values of atom2 and then
! . atom1 and then atom3.
! . Duplicate terms are not removed.
!-----------------------------------------------------------------------------*/
void HarmonicBondContainer_Sort ( HarmonicBondContainer *self )
{
    if ( ( self != NULL ) && ( ! self->isSorted ) )
    {
        auto Integer atom1, atom2, i ;
        /* . Order atom1 and atom2 within each term. */
        for ( i = 0 ; i < self->nTerms ; i++ )
        {
            atom1 = self->terms[i].atom1 ;
            atom2 = self->terms[i].atom2 ;
            if ( atom2 > atom1 )
            {
                self->terms[i].atom1 = atom2 ;
                self->terms[i].atom2 = atom1 ;
            }
        }
        /* . Order the terms within the container. */
        qsort ( ( void * ) self->terms, ( size_t ) self->nTerms, sizeof ( HarmonicBond ), ( void * ) HarmonicBondTerm_Compare ) ;
        self->isSorted = True ;
    }
}

/*------------------------------------------------------------------------------
! . Return the upper bound for the container.
! . This is the value of the largest index plus one.
!-----------------------------------------------------------------------------*/
int HarmonicBondContainer_UpperBound ( HarmonicBondContainer *self )
{
    Integer upperBound = 0 ;
    if ( ( self != NULL ) && ( self->nTerms > 0 ) )
    {
        HarmonicBondContainer_Sort ( self ) ;
        upperBound = self->terms[self->nTerms-1].atom1 + 1 ;
    }
    return upperBound ;
}

/*==============================================================================
! . Private procedures.
!============================================================================*/
static Integer HarmonicBondTerm_Compare ( const void *vTerm1, const void *vTerm2 )
{
    HarmonicBond *term1, *term2 ;
    Integer i ;
    term1 = ( HarmonicBond * ) vTerm1 ;
    term2 = ( HarmonicBond * ) vTerm2 ;
         if ( term1->atom1 < term2->atom1 ) i = -1 ;
    else if ( term1->atom1 > term2->atom1 ) i =  1 ;
    else if ( term1->atom2 < term2->atom2 ) i = -1 ;
    else if ( term1->atom2 > term2->atom2 ) i =  1 ;
    else if ( term1->type  < term2->type  ) i = -1 ;
    else if ( term1->type  > term2->type  ) i =  1 ;
    else if ( ! term1->isActive && term2->isActive ) i = -1 ;
    else if ( term1->isActive && ! term2->isActive ) i =  1 ;
    else i = 0 ;
    return i ;
}

