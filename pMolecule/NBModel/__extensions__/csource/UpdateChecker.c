/*==================================================================================================================================
! . Update checking.
!=================================================================================================================================*/

# include <math.h>

# include "Integer.h"
# include "Memory.h"
# include "UpdateChecker.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _DefaultBuffer 1.5e+00

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
UpdateChecker *UpdateChecker_Allocate ( Status *status )
{
    UpdateChecker *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( UpdateChecker ) ;
        if ( self != NULL ) self->buffer = _DefaultBuffer ;
        else Status_Set ( status, Status_OutOfMemory ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for an image update (those due to changes in the symmetry parameters).
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _NumberOfMinimumImageTranslations 13
Boolean UpdateChecker_CheckForImageUpdate ( const SymmetryParameters     *set1                ,
                                                  SymmetryParameters     *set2                ,
                                                  ImagePairListContainer *images              ,
                                            const Real                    buffer              ,
                                            const Real                    maximumDisplacement )
{
    Boolean doUpdate = False ;
    if ( ( set1 != NULL ) && ( set2 != NULL ) )
    {
        auto Integer        i  ;
        auto Real           dI ;
        auto Matrix33      *dH           = NULL ;
        auto Vector3       *displacement = NULL ;
        auto Integer _MinimumImageTranslations[_NumberOfMinimumImageTranslations][3] = { { -1, -1, -1 } ,
                                                                                         { -1, -1,  0 } ,
                                                                                         { -1, -1,  1 } ,
                                                                                         { -1,  0, -1 } ,
                                                                                         { -1,  0,  0 } ,
                                                                                         { -1,  0,  1 } ,
                                                                                         { -1,  1, -1 } ,
                                                                                         { -1,  1,  0 } ,
                                                                                         { -1,  1,  1 } ,
                                                                                         {  0, -1, -1 } ,
                                                                                         {  0, -1,  0 } ,
                                                                                         {  0, -1,  1 } ,
                                                                                         {  0,  0, -1 } } ;
        dH           = Matrix33_Allocate ( ) ;
        displacement = Vector3_Allocate  ( ) ;
        if ( ( dH != NULL ) && ( displacement != NULL ) )
        {
            Matrix33_CopyTo ( set1->H, dH, NULL ) ;
            Matrix33_Add    ( dH, - 1.0e+00, set2->H, NULL ) ;
            /* . Assume minimum image. */
            if ( images == NULL )
            {
                for ( i = 0 ; i < _NumberOfMinimumImageTranslations ; i++ )
                {
                    Vector3_Item ( displacement, 0 ) = ( Real ) _MinimumImageTranslations[i][0] ;
                    Vector3_Item ( displacement, 1 ) = ( Real ) _MinimumImageTranslations[i][1] ;
                    Vector3_Item ( displacement, 2 ) = ( Real ) _MinimumImageTranslations[i][2] ;
                    Matrix33_ApplyToVector3 ( dH, displacement ) ;
                    dI = Vector3_Norm2 ( displacement ) ;
                    if ( dI > ( buffer - maximumDisplacement ) ) { doUpdate = True ; break ; }
                }
            }
            /* . General case. */
            else
            {
                auto ImagePairList *image ;
                for ( i = 0 ; i < images->count ; i++ )
                {
                    image = images->records[i] ;
                    Vector3_CopyTo ( image->transformation3->translation, displacement, NULL ) ;
                    Vector3_Item ( displacement, 0 ) += ( Real ) image->a ;
                    Vector3_Item ( displacement, 1 ) += ( Real ) image->b ;
                    Vector3_Item ( displacement, 2 ) += ( Real ) image->c ;
                    Matrix33_ApplyToVector3 ( dH, displacement ) ;
                    dI = Vector3_Norm2 ( displacement ) ;
                    if ( dI > ( buffer - maximumDisplacement ) ) { doUpdate = True ; break ; }
                }
            }
        }
        Matrix33_Deallocate ( &dH           ) ;
        Vector3_Deallocate  ( &displacement ) ;
    }
    return doUpdate ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check for a regular update.
! . Updates occur whenever a particle has moved by more than half the buffer distance.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean UpdateChecker_CheckForUpdate  ( const Coordinates3 *set1                ,
                                              Coordinates3 *set2                ,
                                              Selection    *freeAtoms           ,
                                        const Real          buffer              ,
                                              Real         *maximumDisplacement )
{
    Boolean doUpdate = False ;
    if ( ( set1 != NULL ) && ( set2 != NULL ) )
    {
        auto Integer  i, s ;
        auto Real     buffer2, dX, dY, dZ, maxR2 = 0.0e+00, r2, x1, x2, y1, y2, z1, z2 ;
        buffer2 = 0.25e+00 * buffer * buffer ;
        if ( freeAtoms == NULL )
        {
            for ( i = 0 ; i < Coordinates3_Rows ( set1 ) ; i++ )
            {
                Coordinates3_GetRow ( set1, i, x1, y1, z1 ) ;
                Coordinates3_GetRow ( set2, i, x2, y2, z2 ) ;
                dX    = x1 - x2 ;
                dY    = y1 - y2 ;
                dZ    = z1 - z2 ;
                r2    = dX * dX + dY * dY + dZ * dZ ;
                maxR2 = Maximum ( maxR2, r2 ) ;
                if ( r2 > buffer2 ) { doUpdate = True ; break ; }
            }
        }
        else
        {
            for ( s = 0 ; s < Selection_Capacity ( freeAtoms ) ; s++ )
            {
                i = Selection_Item ( freeAtoms, s ) ;
                Coordinates3_GetRow ( set1, i, x1, y1, z1 ) ;
                Coordinates3_GetRow ( set2, i, x2, y2, z2 ) ;
                dX    = x1 - x2 ;
                dY    = y1 - y2 ;
                dZ    = z1 - z2 ;
                r2    = dX * dX + dY * dY + dZ * dZ ;
                maxR2 = Maximum ( maxR2, r2 ) ;
                if ( r2 > buffer2 ) { doUpdate = True ; break ; }
            }
        }
        (*maximumDisplacement) = sqrt ( maxR2 ) ;
    }
    return doUpdate ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void UpdateChecker_Deallocate ( UpdateChecker **self )
{
    if ( (*self) != NULL ) Memory_Deallocate ( (*self) ) ;
}
