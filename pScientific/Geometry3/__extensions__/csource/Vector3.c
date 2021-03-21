/*==================================================================================================================================
! . Procedures for vectors of size 3.
!=================================================================================================================================*/

# include <math.h>
# include <stdlib.h>

# include "Vector3.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
Vector3 *Vector3_Allocate ( void ) { return RealArray1D_AllocateWithExtent ( 3, NULL ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cross-product (in-place).
!---------------------------------------------------------------------------------------------------------------------------------*/
void Vector3_CrossProduct ( Vector3 *self, const Vector3 *other )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Real a0, a1, a2, b0, b1, b2 ;
        a0 = Vector3_Item ( self,  0 ) ;
        a1 = Vector3_Item ( self,  1 ) ;
        a2 = Vector3_Item ( self,  2 ) ;
        b0 = Vector3_Item ( other, 0 ) ;
        b1 = Vector3_Item ( other, 1 ) ;
        b2 = Vector3_Item ( other, 2 ) ;
        Vector3_Item ( self, 0 ) = a1 * b2 - a2 * b1 ;
        Vector3_Item ( self, 1 ) = a2 * b0 - a0 * b2 ;
        Vector3_Item ( self, 2 ) = a0 * b1 - a1 * b0 ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Test for equality of contents.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define TOLERANCE 1.0e-6

Boolean Vector3_IsEqual ( const Vector3 *self, const Vector3 *other )
{
    Boolean QEQUAL = False ;
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Integer i ;
        QEQUAL = True ;
        for ( i = 0 ; i < 3 ; i++ )
        {
            if ( fabs ( self->data[i] - other->data[i] ) > TOLERANCE )
            {
                QEQUAL = False ;
                break ;
            }
        }
    }
    return QEQUAL ;
}

# undef TOLERANCE

/*----------------------------------------------------------------------------------------------------------------------------------
! . Test for nullness.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define TOLERANCE 1.0e-6

Boolean Vector3_IsNull ( const Vector3 *self )
{
    Boolean QNULL = False ;
    if ( self != NULL )
    {
        auto Integer i ;
        QNULL = True ;
        for ( i = 0 ; i < 3 ; i++ )
        {
            if ( fabs ( self->data[i] ) > TOLERANCE )
            {
                QNULL = False ;
                break ;
            }
        }
    }
    return QNULL ;
}

# undef TOLERANCE
