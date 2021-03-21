/*==================================================================================================================================
! . This module implements the random number generator interface.
!=================================================================================================================================*/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "Memory.h"
# include "RandomNumberGenerator.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _DefaultSeed 0

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
RandomNumberGenerator *RandomNumberGenerator_Allocate ( const RandomNumberGeneratorType *type )
{
    RandomNumberGenerator *self = Memory_AllocateType ( RandomNumberGenerator ) ;
    if ( self != NULL )
    {
        self->vState = ( type->Allocate ) ( ) ;
        if ( self->vState == NULL ) RandomNumberGenerator_Deallocate ( &self ) ;
        else
        {
            self->gaussian    = 0.0   ;
            self->hasGaussian = False ;
            self->type        = type  ;
            RandomNumberGenerator_SetSeed ( self, _DefaultSeed ) ;
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
RandomNumberGenerator *RandomNumberGenerator_Clone ( const RandomNumberGenerator *self )
{
    RandomNumberGenerator *new = NULL ;
    if ( self != NULL )
    {
        new = RandomNumberGenerator_Allocate ( self->type ) ;
        if ( new != NULL )
        {
            new->gaussian    = self->gaussian    ;
            new->hasGaussian = self->hasGaussian ;
            Memory_CopyTo ( self->vState, new->vState, self->type->size ) ;
        }
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RandomNumberGenerator_Deallocate ( RandomNumberGenerator **self )
{
    if ( (*self) != NULL )
    {   
        ( (*self)->type->Deallocate ) ( &((*self)->vState) ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Next cardinal.
!---------------------------------------------------------------------------------------------------------------------------------*/
Cardinal RandomNumberGenerator_NextCardinal ( const RandomNumberGenerator *self )
{
    return ( self->type->NextCardinal ) ( self->vState ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Next real on the half-open interval [0,1).
!---------------------------------------------------------------------------------------------------------------------------------*/
Real RandomNumberGenerator_NextReal ( const RandomNumberGenerator *self )
{
    return ( self->type->NextReal ) ( self->vState ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Next real on the open interval (0,1).
!---------------------------------------------------------------------------------------------------------------------------------*/
Real RandomNumberGenerator_NextRealOpen ( const RandomNumberGenerator *self )
{
    Real x ;
    do { x = self->type->NextReal ( self->vState ) ; } while ( x == 0.0 ) ;
    return x ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set the seed.
!---------------------------------------------------------------------------------------------------------------------------------*/
void RandomNumberGenerator_SetSeed ( const RandomNumberGenerator *self, Cardinal seed )
{
    ( self->type->SetSeed ) ( self->vState, seed ) ;
}
