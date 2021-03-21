/*==================================================================================================================================
! . A module to handle a container for transformations.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>

# include "IntegerUtilities.h"
# include "Transformation3Container.h"
# include "Memory.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
Transformation3Container *Transformation3Container_Allocate ( const Integer  capacity, Status *status )
{
    Transformation3Container *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( Transformation3Container ) ;
        if ( self != NULL )
        {
	    self->capacity = capacity         ;
            self->identity = -1 ;
            self->inverses = NULL             ;
            self->isOwner  = False            ;
	    self->items    = Memory_AllocateArrayOfReferences ( capacity, Transformation3 ) ;
            if ( self->items == NULL ) Transformation3Container_Deallocate ( &self ) ;
            else
            {
                auto Integer  i ;
                for ( i = 0 ; i < capacity ; i++ ) self->items[i] = NULL ;
            }
        }
        if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Transformation3Container_Deallocate ( Transformation3Container **self )
{
    if ( (*self) != NULL )
    {
        if ( (*self)->isOwner )
        {
            auto Integer  i ;
	    for ( i = 0 ; i < (*self)->capacity ; i++ ) Transformation3_Deallocate ( &((*self)->items[i]) ) ;
        }
        Integer_Deallocate ( &((*self)->inverses) ) ;
        Memory_Deallocate   (   (*self)->items     ) ;
        Memory_Deallocate   (   (*self)            ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find the identity.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Transformation3Container_FindIdentity ( Transformation3Container *self )
{
    if ( ( self != NULL ) && ( self->capacity > 0 ) )
    {
        auto Integer  i ;
        for ( i = 0 ; i < self->capacity ; i++ )
        {
            if ( Transformation3_IsIdentity ( self->items[i] ) )
            {
                self->identity = i ;
                break ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find the inverse translation for a transformation (in terms of a, b, c).
!---------------------------------------------------------------------------------------------------------------------------------*/
# define BIGNUM    -999999
# define TOLERANCE  1.0e-4
void Transformation3Container_FindInverseIntegerTranslation ( const Transformation3Container *self        ,
                                                              const Integer                   t           ,
                                                              const Integer                   a           ,
                                                              const Integer                   b           ,
                                                              const Integer                   c           ,
                                                                    Vector3                  *translation ,
                                                                    Integer                  *aInverse    ,
                                                                    Integer                  *bInverse    ,
                                                                    Integer                  *cInverse    )
{
    (*aInverse) = BIGNUM ;
    (*bInverse) = BIGNUM ;
    (*cInverse) = BIGNUM ;
    if ( ( self != NULL ) && ( translation != NULL ) && ( t >= 0 ) && ( self->inverses[t] >= 0 ) )
    {
        auto Integer  i, tInverse ;
        auto Integer  n[3], v     ;
        tInverse = self->inverses[t] ;
        /* . Find the inverse translation minus the inverse's translation. */
        Vector3_CopyTo ( self->items[t]->translation, translation, NULL ) ;
        translation->data[0] += ( Real ) a ;
        translation->data[1] += ( Real ) b ;
        translation->data[2] += ( Real ) c ;
        Matrix33_ApplyToVector3 ( self->items[tInverse]->rotation, translation ) ;
        Vector3_Scale           ( translation, -1.0e+00 ) ;
        Vector3_Add  ( translation, -1.0e+00, self->items[tInverse]->translation, NULL ) ;
        /* . Convert the result to integers. */
        for ( i = 0 ; i < 3 ; i++ )
        {
            v = Real_RoundToInteger ( translation->data[i] ) ;
            if ( fabs ( translation->data[i] - ( Real ) v ) < TOLERANCE ) n[i] = v      ;
            else                                                          n[i] = BIGNUM ;
        }
        (*aInverse) = n[0] ;
        (*bInverse) = n[1] ;
        (*cInverse) = n[2] ;
    }
}
# undef BIGNUM
# undef TOLERANCE

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find the inverses - only the rotational part!
! . This should probably be changed (have both rotational and full inverses).
!---------------------------------------------------------------------------------------------------------------------------------*/
void Transformation3Container_FindInverses ( Transformation3Container *self, Status *status )
{
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto Matrix33 *m    ;
        Integer_Deallocate ( &(self->inverses) ) ;
        m              = Matrix33_Allocate ( ) ; Matrix33_Set ( m, 0.0e+00 ) ;
        self->inverses = Integer_Allocate ( self->capacity, status ) ;
        if ( ( m != NULL ) && ( self->inverses != NULL ) )
        {
            auto Integer   i, j ;
            Integer_Set ( self->inverses, self->capacity, -1 ) ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                if ( self->inverses[i] == -1 )
                {
                    Matrix33_Invert ( m, self->items[i]->rotation ) ;
                    for ( j = 0 ; j <= i ; j++ )
                    {
                        if ( self->inverses[j] == -1 )
                        {
                            if ( Matrix33_IsEqual ( m, self->items[j]->rotation ) )
                            {
                                self->inverses[i] = j ;
                                self->inverses[j] = i ;
                                break ;
                            }
                        }
                    }
                }
            }
        }
        else Status_Set ( status, Status_OutOfMemory ) ;
        Matrix33_Deallocate ( &m ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Various counters.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer  Transformation3Container_NumberOfNonSelfInverses ( const Transformation3Container *self )
{
    Integer  n = 0 ;
    if ( ( self != NULL ) && ( self->capacity > 0 ) && ( self->inverses != NULL ) )
    {
        auto Integer  i, t ;
        for ( i = 0 ; i < self->capacity ; i++ )
        {
            t = self->inverses[i] ;
            if ( ( t < -1 ) && ( t != i ) ) n += 1 ;
        }
    }
    return n ;
}

Integer  Transformation3Container_NumberOfSelfInverses ( const Transformation3Container *self )
{
    Integer  n = 0 ;
    if ( ( self != NULL ) && ( self->capacity > 0 ) && ( self->inverses != NULL ) )
    {
        auto Integer  i ;
        for ( i = 0 ; i < self->capacity ; i++ ) { if ( self->inverses[i] == i ) n += 1 ; }
    }
    return n ;
}
