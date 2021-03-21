/*==================================================================================================================================
! . A module to handle transformations.
!=================================================================================================================================*/

# include "Memory.h"
# include "Transformation3.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
Transformation3 *Transformation3_Allocate ( Status *status )
{
    Transformation3 *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( Transformation3 ) ;
        if ( self != NULL )
        {
            self->rotation    = NULL ;
            self->translation = NULL ;
        }
        else Status_Set ( status, Status_OutOfMemory ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation with uninitialized rotation and translation.
!---------------------------------------------------------------------------------------------------------------------------------*/
Transformation3 *Transformation3_AllocateFull ( Status *status )
{
    Transformation3 *self = Transformation3_Allocate ( status ) ;
    if ( self != NULL )
    {
        self->rotation    = Matrix33_Allocate ( ) ;
        self->translation = Vector3_Allocate  ( ) ;
        if ( ( self->rotation == NULL ) || ( self->translation == NULL ) )
        {
            Transformation3_Deallocate ( &self ) ;
            Status_Set ( status, Status_OutOfMemory ) ;
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Apply to a vector.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Transformation3_ApplyToVector3 ( const Transformation3 *self, Vector3 *a )
{
    if ( ( self != NULL ) && ( a != NULL ) )
    {
        Matrix33_ApplyToVector3 ( self->rotation, a ) ;
        Vector3_Add  ( a, 1.0e+00, self->translation, NULL ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
Transformation3 *Transformation3_Clone ( const Transformation3 *self, Status *status )
{
    Transformation3 *new = NULL ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        new = Transformation3_AllocateFull ( status ) ;
        if ( new != NULL )
        {
            Matrix33_CopyTo ( self->rotation,    new->rotation   , NULL ) ;
            Vector3_CopyTo  ( self->translation, new->translation, NULL ) ;
        }
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copying.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Transformation3_CopyTo ( const Transformation3 *self, Transformation3 *other )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        Matrix33_CopyTo ( self->rotation   , other->rotation,    NULL ) ;
        Vector3_CopyTo  ( self->translation, other->translation, NULL ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Transformation3_Deallocate ( Transformation3 **self )
{
    if ( (*self) != NULL )
    {
        Matrix33_Deallocate ( &((*self)->rotation)    ) ;
        Vector3_Deallocate  ( &((*self)->translation) ) ;
        Memory_Deallocate   ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Constructor from a rotation and translation.
!---------------------------------------------------------------------------------------------------------------------------------*/
Transformation3 *Transformation3_FromRotationTranslation ( Matrix33 *rotation    ,
                                                           Vector3  *translation ,
                                                           Status   *status      )
{
    Transformation3 *self = NULL ;
    if ( ( rotation != NULL ) && ( translation != NULL) && Status_IsOK ( status ) )
    {
        self = Transformation3_Allocate ( status ) ;
        if ( self != NULL )
        {
            self->rotation    = Matrix33_CloneShallow ( rotation   , status ) ;
            self->translation = Vector3_CloneShallow  ( translation, status ) ;
            if ( ! Status_IsOK ( status ) ) Transformation3_Deallocate ( &self ) ;
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Test for equality of contents.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean Transformation3_IsEqual ( const Transformation3 *self, const Transformation3 *other )
{
    Boolean isEqual = False ;
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        isEqual = ( Matrix33_IsEqual ( self->rotation, other->rotation ) && Vector3_IsEqual ( self->translation, other->translation ) ) ;
    }
    return isEqual ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Test for identity.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean Transformation3_IsIdentity ( const Transformation3 *self )
{
    Boolean isIdentity = False ;
    if ( self != NULL ) isIdentity = ( Matrix33_IsIdentity ( self->rotation ) && Vector3_IsNull ( self->translation ) ) ;
    return isIdentity ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Orthogonalization (transform to M * R * M^-1 + M t).
!---------------------------------------------------------------------------------------------------------------------------------*/
void Transformation3_Orthogonalize ( Transformation3 *self, const Matrix33 *A, const Matrix33 *B )
{
    if ( ( self != NULL ) && ( A != NULL ) && ( B != NULL ) )
    {
        Matrix33_PreMultiplyBy  ( self->rotation, A ) ;
        Matrix33_PostMultiplyBy ( self->rotation, B ) ;
        Matrix33_ApplyToVector3 ( A, self->translation ) ;
    }
}

