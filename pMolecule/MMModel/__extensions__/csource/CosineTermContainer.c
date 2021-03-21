/*==================================================================================================================================
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "Boolean.h"
# include "BooleanUtilities.h"
# include "CosineTermContainer.h"
# include "Integer.h"
# include "Memory.h"
# include "NumericalMacros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Activate terms.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CosineTermContainer_ActivateTerms ( CosineTermContainer *self )
{
    if ( self != NULL )
    {
        auto Integer i ;
	for ( i = 0 ; i < self->nTerms ; i++ ) self->terms[i].isActive = True ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
CosineTermContainer *CosineTermContainer_Allocate ( const Integer nIndices, const Integer nTerms, const Integer nParameters )
{
    CosineTermContainer *self = NULL ;
    if ( ( nIndices > 0 ) && ( nTerms != 0 ) && ( nParameters != 0 ) )
    {
        Integer i ;
	self = Memory_AllocateType ( CosineTermContainer ) ;
        self->nIndices    = nIndices    ;
	self->nParameters = nParameters ;
 	self->nTerms      = nTerms      ;
	self->parameters  = Memory_AllocateArrayOfTypes ( nParameters, CosineParameter ) ;
	self->terms	  = Memory_AllocateArrayOfTypes ( nTerms     , CosineTerm      ) ;
        for ( i = 0 ; i < nParameters ; i++ ) CosineParameter_Initialize ( &(self->parameters[i])           ) ;
	for ( i = 0 ; i < nTerms      ; i++ ) CosineTerm_Allocate        ( &(self->terms     [i]), nIndices ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
CosineTermContainer *CosineTermContainer_Clone ( const CosineTermContainer *self )
{
    CosineTermContainer *new = NULL ;
    if ( self != NULL )
    {
        auto Integer i ;
        new = CosineTermContainer_Allocate ( self->nIndices, self->nTerms, self->nParameters ) ;
        for ( i = 0 ; i < self->nParameters ; i++ ) CosineParameter_Clone ( &(new->parameters[i]), &(self->parameters[i]) ) ;
        for ( i = 0 ; i < self->nTerms      ; i++ ) CosineTerm_Clone      ( &(new->terms     [i]), &(self->terms     [i]) , self->nIndices ) ;
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deactivate terms that do not involve any atoms in the selection.
! . Already deactivated terms are not affected.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CosineTermContainer_DeactivateTerms ( CosineTermContainer *self, Selection *selection )
{
    if ( ( self != NULL ) && ( selection != NULL ) )
    {
        auto Boolean       isActive ;
        auto BooleanBlock *flags    ;
        auto Integer       i, n, t  ;
        n     = CosineTermContainer_UpperBound ( self ) ;
        flags = Selection_MakeFlags ( selection, n, NULL ) ;
	for ( t = 0 ; t < self->nTerms ; t++ )
	{
            if ( self->terms[t].isActive )
            {
                isActive = False ;
                for ( i = 0 ; i < self->nIndices ; i++ )
                {
                    isActive = ( isActive || Block_Item ( flags, self->terms[t].indices[i] ) ) ; /* . Only one atom need be active. */
                }
                self->terms[t].isActive = isActive ;
            }
	}
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CosineTermContainer_Deallocate ( CosineTermContainer **self )
{
    if ( (*self) != NULL )
    {
        auto Integer i ;
        for ( i = 0 ; i < (*self)->nParameters ; i++ ) CosineParameter_Deallocate ( &((*self)->parameters[i]) ) ;
        for ( i = 0 ; i < (*self)->nTerms      ; i++ ) CosineTerm_Deallocate      ( &((*self)->terms     [i]) ) ;
        Memory_Deallocate ( (*self)->parameters ) ;
        Memory_Deallocate ( (*self)->terms      ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find the maximum period.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer CosineTermContainer_FindMaximumPeriod ( CosineTermContainer *self )
{
    Integer p = -1 ;
    if ( self != NULL )
    {
        auto Integer i, n ;
        for ( i = 0 ; i < self->nParameters ; i++ )
        {
            for ( n = 0 ; n < self->parameters[i].nTerms ; n++ )
            {
                p = Maximum ( p, self->parameters[i].periods[n] ) ;
            }
        }
    }
    return p ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the powers representation of the parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CosineTermContainer_MakePowers ( CosineTermContainer *self )
{
    if ( self != NULL )
    {
        auto Integer i ;
        for ( i = 0 ; i < self->nParameters ; i++ ) CosineParameter_MakePowers ( &(self->parameters[i]) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Return the number of inactive terms.
!---------------------------------------------------------------------------------------------------------------------------------*/
int CosineTermContainer_NumberOfInactiveTerms ( const CosineTermContainer *self )
{
    Integer n = 0 ;
    if ( self != NULL )
    {
        auto Integer i ;
	for ( i = 0 ; i < self->nTerms ; i++ ) if ( ! self->terms[i].isActive ) n++ ;
    }
    return n ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Pruning.
! . Only terms are pruned.
!---------------------------------------------------------------------------------------------------------------------------------*/
CosineTermContainer *CosineTermContainer_Prune ( CosineTermContainer *self, Selection *selection )
{
    CosineTermContainer *new = NULL ;
    if ( ( self != NULL ) && ( selection != NULL ) )
    {
        auto Boolean        flagged, *toKeep ;
        auto BooleanBlock  *flags     ;
        auto IntegerBlock  *positions ;
        auto Integer        i, n, t   ;
        n         = CosineTermContainer_UpperBound ( self ) ;
        flags     = Selection_MakeFlags     ( selection, n, NULL ) ;
        positions = Selection_MakePositions ( selection, n, NULL ) ;
        toKeep    = Boolean_Allocate        ( self->nTerms, NULL ) ;
	for ( t = n = 0 ; t < self->nTerms ; t++ )
	{
            flagged = True ;
            for ( i = 0 ; i < self->nIndices ; i++ )
            {
                flagged = ( flagged && Block_Item ( flags, self->terms[t].indices[i] ) ) ; /* . All atoms must be to keep. */
            }
            toKeep[t] = flagged ;
            if ( toKeep[t] ) n++ ;
	}
	if ( n > 0 )
	{
            new = CosineTermContainer_Allocate ( self->nIndices, n, self->nParameters ) ;
            for ( i = 0 ; i < self->nParameters ; i++ ) CosineParameter_Clone ( &(new->parameters[i]), &(self->parameters[i]) ) ;
            for ( t = n = 0 ; t < self->nTerms ; t++ )
            {
                if ( toKeep[t] )
                {
                    for ( i = 0 ; i < self->nIndices ; i++ )
                    {
        	        new->terms[n].indices[i] = Block_Item ( positions, self->terms[t].indices[i] ) ;
                    }
        	    new->terms[n].isActive = self->terms[t].isActive ;
        	    new->terms[n].type     = self->terms[t].type ;
        	    n++ ;
                }
            }
	}
	Boolean_Deallocate ( &toKeep ) ;
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Return the upper bound for the container.
! . This is the value of the largest index plus one.
!---------------------------------------------------------------------------------------------------------------------------------*/
int CosineTermContainer_UpperBound ( CosineTermContainer *self )
{
    Integer upperBound = -1 ;
    if ( ( self != NULL ) && ( self->nTerms > 0 ) )
    {
        auto Integer i, t ;
	for ( t = 0 ; t < self->nTerms ; t++ )
	{
            for ( i = 0 ; i < self->nIndices ; i++ ) { upperBound = Maximum ( upperBound, self->terms[t].indices[i] ) ; }
        }
        upperBound += 1 ;
    }
    return upperBound ;
}
