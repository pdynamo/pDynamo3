/*==================================================================================================================================l
! . A module to handle image pair-lists.
!=================================================================================================================================*/

# include <stdio.h>

# include "ImagePairListContainer.h"
# include "Memory.h"
# include "NumericalMacros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _GrowthFactor    2.0e+00
# define _MinimumCapacity 32

/*==================================================================================================================================
! . ImagePairList functions.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
ImagePairList *ImagePairList_Allocate ( Status *status )
{
    ImagePairList *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( ImagePairList ) ;
        if( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
        else
        {
            self->a               = 0 ;
            self->b               = 0 ;
            self->c               = 0 ;
            self->scale           = 0.0e+00 ;
            self->pairList        = NULL ;
            self->transformation3 = NULL ;
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ImagePairList_Deallocate ( ImagePairList **self )
{
/*printf ( "Image Pair List IN> %p\n", *self ) ; fflush ( stdout ) ; */
    if ( (*self) != NULL )
    {
        PairList_Deallocate ( &((*self)->pairList) ) ;
        (*self)->transformation3 = NULL ;
        Memory_Deallocate ( (*self) ) ;
    }
/*printf ( "Image Pair List OUT> %p\n", *self ) ; fflush ( stdout ) ; */
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Constructor from items.
!---------------------------------------------------------------------------------------------------------------------------------*/
ImagePairList *ImagePairList_FromItems ( const Integer                  a               ,
                                         const Integer                  b               ,
                                         const Integer                  c               ,
                                         const Real                     scale           ,
                                               PairList                *pairList        ,
                                               Transformation3         *transformation3 ,
                                               Status                  *status          )
{
    ImagePairList *self = NULL ;
    if ( ( PairList_NumberOfPairs ( pairList ) > 0 ) && Status_IsOK ( status ) )
    {
        self = ImagePairList_Allocate ( status ) ;
        if ( self != NULL )
        {
            self->a               = a ;
            self->b               = b ;
            self->c               = c ;
            self->scale           = scale           ;
            self->pairList        = pairList        ;
            self->transformation3 = transformation3 ;
        }
    }
    return self ;
}

/*==================================================================================================================================
! . ImagePairListContainer functions.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
ImagePairListContainer *ImagePairListContainer_Allocate ( Status *status )
{
    ImagePairListContainer *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( ImagePairListContainer ) ;
        if ( self != NULL )
        {
            self->records = Memory_AllocateArrayOfReferences ( _MinimumCapacity, ImagePairList ) ;
            if ( self->records == NULL ) ImagePairListContainer_Deallocate ( &self ) ;
            else
            {
                auto Integer  i ;
                for ( i = 0 ; i < _MinimumCapacity ; i++ ) self->records[i] = NULL ;
                self->capacity      = _MinimumCapacity ;
                self->count         = 0 ;
                self->numberOfPairs = 0 ;
            }
        }
        if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Append a record.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ImagePairListContainer_Append ( ImagePairListContainer *self   ,
                                     ImagePairList          *record ,
                                     Status                 *status )
{
    if ( ( self != NULL ) && ( record != NULL ) && Status_IsOK ( status ) )
    {
        auto Boolean isOK = True ;
        if ( self->count >= self->capacity )
        {
            isOK = ImagePairListContainer_Reallocate ( self, ( Integer  ) ( self->capacity * _GrowthFactor ), status ) ;
        }
        if ( isOK )
        {
            self->records[self->count] = record ;
            self->count               += 1 ;
            self->numberOfPairs       += PairList_NumberOfPairs ( record->pairList ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Constructor.
!---------------------------------------------------------------------------------------------------------------------------------*/
ImagePairListContainer *ImagePairListContainer_Constructor ( const PairListGenerator        *generator          ,
                                                                   Selection                *atomsA             ,
                                                                   Selection                *atomsB             ,
                                                                   Selection                *freeAtoms          ,
                                                             const Coordinates3             *coordinates3A      ,
                                                             const Coordinates3             *coordinates3B      ,
                                                             const SymmetryParameters       *symmetryParameters ,
                                                             const Transformation3Container *transformations    ,
                                                             const ImageScanContainer       *scanData           ,
                                                                   RegularGrid              *gridA              ,
                                                                   RegularGridOccupancy     *occupancyA         ,
                                                             const Boolean                   checkForInverses   ,
                                                                   Status                   *status             )
{
    ImagePairListContainer *self = NULL ;
    if ( ( generator          != NULL ) &&
         ( generator->cutOff  >  0.0  ) &&
         ( coordinates3A      != NULL ) &&
         ( coordinates3B      != NULL ) &&
         ( scanData           != NULL ) &&
         ( symmetryParameters != NULL ) &&
         ( transformations    != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Boolean          isOK ;
        auto Coordinates3    *iCoordinates3    = NULL ;
        auto Transformation3 *iTransformation3 = NULL ;
        auto Vector3         *displacement     = NULL ;
        iCoordinates3    = Coordinates3_Allocate           ( Coordinates3_Rows ( coordinates3B ), status ) ;
        self             = ImagePairListContainer_Allocate ( status ) ;
        iTransformation3 = Transformation3_AllocateFull    ( status ) ;
        displacement     = Vector3_Allocate                ( ) ;
        isOK             = ( displacement     != NULL ) &&
                           ( iCoordinates3    != NULL ) &&
                           ( self             != NULL ) &&
                           ( iTransformation3 != NULL ) ;
        if ( isOK )
        {
            auto Integer        r, t, tLast = -1 ;
            auto ImagePairList *image  ;
            auto ImageScan     *record ;
            auto PairList      *pairList = NULL ;
            for ( r = 0 ; r < scanData->count ; r++ )
            {
                record = &(scanData->records[r]) ;
                t      = record->t ;
                if ( checkForInverses && record->doSkip ) continue ;
                if ( t != tLast )
                {
                    Transformation3_CopyTo        ( transformations->items[t], iTransformation3 ) ;
                    Transformation3_Orthogonalize ( iTransformation3, symmetryParameters->H, symmetryParameters->inverseH ) ;
                    Coordinates3_CopyTo           ( coordinates3B , iCoordinates3    , NULL   ) ;
                    Coordinates3_Transform        ( iCoordinates3 , iTransformation3 , atomsB ) ; /* . Worth doing or just NULL? */
                    tLast = t ;
                }
                SymmetryParameters_Displacement ( symmetryParameters, record->a, record->b, record->c, displacement ) ;
                Coordinates3_Translate ( iCoordinates3, displacement, NULL ) ;
                pairList = PairListGenerator_CrossPairListFromDoubleCoordinates3 ( generator     ,
                                                                                   coordinates3A ,
                                                                                   iCoordinates3 ,
                                                                                   NULL          ,
                                                                                   NULL          ,
                                                                                   atomsA        ,
                                                                                   atomsB        ,
                                                                                   freeAtoms     ,
                                                                                   freeAtoms     ,
                                                                                   NULL          ,
                                                                                   gridA         ,
                                                                                   occupancyA    ,
                                                                                   status        ) ;
                image = ImagePairList_FromItems ( record->a, record->b, record->c, record->scale, pairList, transformations->items[t], status ) ;
                ImagePairListContainer_Append ( self, image, status ) ;
                if ( ! Status_IsOK ( status ) ) goto FinishUp ;
                Vector3_Scale          ( displacement, -1.0e+00 ) ;
                Coordinates3_Translate ( iCoordinates3, displacement, NULL ) ;
	    }
        }
        else Status_Set ( status, Status_OutOfMemory ) ;
    FinishUp:
        Coordinates3_Deallocate    ( &iCoordinates3    ) ;
        Transformation3_Deallocate ( &iTransformation3 ) ;
        Vector3_Deallocate         ( &displacement     ) ;
        if ( ! Status_IsOK ( status ) ) ImagePairListContainer_Deallocate ( &self ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ImagePairListContainer_Deallocate ( ImagePairListContainer **self )
{
    if ( (*self) != NULL )
    {
        auto Integer  i ;
        for ( i = 0 ; i < (*self)->count ; i++ )
        {
            ImagePairList_Deallocate ( &((*self)->records[i]) ) ;
        }
        Memory_Deallocate ( (*self)->records ) ;
        Memory_Deallocate ( (*self)          ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Various counters.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer ImagePairListContainer_NumberOfImages ( const ImagePairListContainer *self ) { return ( ( self == NULL ) ? 0 : self->count ) ; }

Integer ImagePairListContainer_NumberOfPairs ( const ImagePairListContainer *self ) { return ( ( self == NULL ) ? 0 : self->numberOfPairs ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Reallocate records.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean ImagePairListContainer_Reallocate ( ImagePairListContainer *self, const Integer  capacity, Status *status )
{
    Boolean isOK = True ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer  n ;
             if ( capacity > self->capacity ) n = capacity    ;
        else if ( capacity < self->count    ) n = self->count ;
        n = Maximum ( capacity, _MinimumCapacity ) ;
        if ( n != self->capacity )
        {
            auto ImagePairList **records  ;
            records = Memory_ReallocateArrayOfReferences ( self->records, n, ImagePairList ) ;
            if ( records != NULL )
            {
                if ( n > self->count ) { auto Integer  i ; for ( i = self->count ; i < n ; i++ ) records[i] = NULL ; }
                self->capacity = n ;
                self->records  = records ;
            }
            else isOK = False ;
        }
    }
    return isOK ;
}

/*==================================================================================================================================
! . ImagePairListIterator functions.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Finalization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ImagePairListIterator_Finalize ( ImagePairListIterator *self )
{
    if ( self != NULL )
    {
        self->current = self->target->count + 1 ;
        Coordinates3_Deallocate    ( &(self->iCoordinates3   ) ) ;
        Coordinates3_Deallocate    ( &(self->iGradients3     ) ) ;
        Transformation3_Deallocate ( &(self->iTransformation3) ) ;
        Transformation3_Deallocate ( &(self->xTransformation3) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ImagePairListIterator_Gradients ( ImagePairListIterator *self )
{
    if ( ( self != NULL ) && ( self->doGradients ) )
    {
        SymmetryParameterGradients_ImageDerivatives ( self->symmetryParameterGradients ,
                                                      self->symmetryParameters         ,
                                                      self->xTransformation3           , 
                                                      self->coordinates3               ,
                                                      self->iGradients3                ) ;
        Matrix33_Transpose          ( self->iTransformation3->rotation, NULL ) ;
        Coordinates3_Rotate         ( self->iGradients3, self->iTransformation3->rotation, NULL ) ;
        Coordinates3_Add ( self->gradients3, 1.0e+00, self->iGradients3, NULL ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ImagePairListIterator_Initialize ( ImagePairListIterator      *self                       ,
                                        ImagePairListContainer     *target                     ,
                                        Coordinates3               *coordinates3               ,
                                        SymmetryParameters         *symmetryParameters         ,
                                        Coordinates3               *gradients3                 ,
                                        SymmetryParameterGradients *symmetryParameterGradients ,
                                        Status                     *status                     )
{
    if ( ( self               != NULL ) &&
         ( target             != NULL ) &&
         ( coordinates3       != NULL ) &&
         ( symmetryParameters != NULL ) &&
         Status_IsOK ( status ) )
    {
        /* . Initialization. */
        self->current     = 0 ;
        self->doGradients = ( gradients3 != NULL ) && ( symmetryParameterGradients != NULL ) ;
        self->scale       = 1.0e+00 ;
        /* . Aliases. */
        self->coordinates3                = coordinates3               ;
        self->gradients3                  = gradients3                 ;
        self->target                      = target                     ;
        self->pairList                    = NULL                       ;
        self->symmetryParameterGradients  = symmetryParameterGradients ;
        self->symmetryParameters          = symmetryParameters         ;
        /* . Allocate space. */
        self->iCoordinates3    = Coordinates3_Allocate ( Coordinates3_Rows ( coordinates3 ), status ) ;
        self->iTransformation3 = Transformation3_AllocateFull ( status ) ;
        self->xTransformation3 = Transformation3_AllocateFull ( status ) ;
        if ( self->doGradients ) self->iGradients3 = Coordinates3_Allocate ( Coordinates3_Rows ( coordinates3 ), status ) ;
        else                     self->iGradients3 = NULL ;
        if ( ( self->iCoordinates3    == NULL ) ||
             ( self->iTransformation3 == NULL ) ||
             ( self->xTransformation3 == NULL ) ||
             ( self->doGradients && ( self->iGradients3 == NULL ) ) )
        {
            ImagePairListIterator_Finalize ( self ) ;
            Status_Set ( status, Status_OutOfMemory ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Next iteration.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean ImagePairListIterator_Next ( ImagePairListIterator *self )
{
    Boolean toContinue = False ;
    if ( ( self != NULL ) && ( self->current < self->target->count ) )
    {
        auto ImagePairList *record = self->target->records[self->current] ;
        /* . Basic data. */
        self->pairList = record->pairList ;
        self->scale    = record->scale    ;
        /* . Generate the image transformation in real space. */
        Transformation3_CopyTo        ( record->transformation3, self->xTransformation3 ) ;
        self->xTransformation3->translation->data[0] += ( Real ) record->a ;
        self->xTransformation3->translation->data[1] += ( Real ) record->b ;
        self->xTransformation3->translation->data[2] += ( Real ) record->c ;
        Transformation3_CopyTo        ( self->xTransformation3, self->iTransformation3 ) ;
        Transformation3_Orthogonalize ( self->iTransformation3, self->symmetryParameters->H, self->symmetryParameters->inverseH ) ;
        /* . Generate image coordinates. */
        Coordinates3_CopyTo    ( self->coordinates3 , self->iCoordinates3   , NULL ) ;
        Coordinates3_Transform ( self->iCoordinates3, self->iTransformation3, NULL ) ;
        /* . Gradients. */
        if ( self->doGradients ) Coordinates3_Set ( self->iGradients3, 0.0e+00 ) ;
        /* . Finish up. */
        self->current += 1    ;
        toContinue     = True ;
    }
    return toContinue ;
}
