/*==================================================================================================================================
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "CosineTerm.h"
# include "Memory.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Term internal allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CosineTerm_Allocate ( CosineTerm *self, const Integer nIndices )
{
    CosineTerm_Initialize ( self ) ;
    if ( ( self != NULL ) && ( nIndices > 0 ) )
    {
        self->indices = Memory_AllocateArrayOfTypes ( nIndices, Integer ) ;
        Memory_Set ( self->indices, nIndices, 0 ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Term internal clone.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CosineTerm_Clone ( CosineTerm *self, const CosineTerm *other, const Integer nIndices )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Integer i ;
        CosineTerm_Allocate ( self, nIndices ) ;
        self->isActive = other->isActive ;
        self->type     = other->type     ;
        for ( i = 0 ; i < nIndices ; i++ ) { self->indices[i] = other->indices[i] ; }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Term internal deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CosineTerm_Deallocate ( CosineTerm *self )
{
    if ( self != NULL )
    {
        self->isActive = False ;
        self->type     = -1    ;
        Memory_Deallocate ( self->indices ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Term internal initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CosineTerm_Initialize ( CosineTerm *self )
{
    if ( self != NULL )
    {
        self->indices  = NULL  ;
        self->isActive = False ;
        self->type     = -1    ;
    }
}
