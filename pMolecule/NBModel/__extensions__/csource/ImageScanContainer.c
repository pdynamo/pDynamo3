/*==================================================================================================================================l
! . A module to handle images to scan when generating image pair-lists.
!=================================================================================================================================*/

# include "ImageScanContainer.h"
# include "Memory.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _GrowthFactor    2.0e+00
# define _MinimumCapacity 32

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
ImageScanContainer *ImageScanContainer_Allocate ( Status *status )
{
    ImageScanContainer *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( ImageScanContainer ) ;
        if ( self != NULL )
        { 
            self->records = Memory_AllocateArrayOfTypes ( _MinimumCapacity, ImageScan ) ;
            if ( self->records == NULL ) ImageScanContainer_Deallocate ( &self ) ;
            else
            {
                self->capacity = _MinimumCapacity ;
                self->count    = 0 ;
                self->cutOff   = 0.0e+00 ;
            }
        }
        if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Append a record.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ImageScanContainer_Append ( ImageScanContainer *self   ,
                                 Boolean             doSkip ,
                                 Integer             t      ,
                                 Integer             a      ,
                                 Integer             b      , 
                                 Integer             c      ,
                                 Real                scale  ,
                                 Status             *status )
{
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto Boolean isOK = True ;
        if ( self->count >= self->capacity )
        {
            isOK = ImageScanContainer_Reallocate ( self, ( Integer  ) ( self->capacity * _GrowthFactor ), status ) ;
        }
        if ( isOK )
        {
            auto ImageScan *record = &(self->records[self->count]) ;
            record->doSkip  = doSkip ;
            record->t       = t      ;
            record->a       = a      ;
            record->b       = b      ;
            record->c       = c      ;
            record->scale   = scale  ;
            self->count    += 1      ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Constructor from coordinates and symmetry information.
!---------------------------------------------------------------------------------------------------------------------------------*/
ImageScanContainer *ImageScanContainer_Constructor ( const Real                      cutOff             ,
                                                     const Coordinates3             *coordinates3       ,
                                                     const SymmetryParameters       *symmetryParameters ,
                                                     const Transformation3Container *transformations    ,
                                                     const Boolean                   checkForInverses   ,
                                                     const Integer                   expandFactor       ,
                                                           Status                   *status             )
{
    ImageScanContainer *self = NULL ;
    if ( ( coordinates3       != NULL ) &&
         ( symmetryParameters != NULL ) &&
         ( transformations    != NULL ) &&
         ( cutOff             >  0.0  ) &&
         Status_IsOK ( status ) )
    {
        auto Boolean          isOK ;
        auto Transformation3 *iTransformation3 = NULL ;
        auto Vector3         *displacement     = NULL ,
                             *iLower           = NULL ,
                             *iUpper           = NULL ,
                             *lower            = NULL ,
                             *upper            = NULL ;
        self             = ImageScanContainer_Allocate  ( status ) ;
        iTransformation3 = Transformation3_AllocateFull ( status ) ;
        displacement     = Vector3_Allocate ( ) ;
        iLower           = Vector3_Allocate ( ) ;
        iUpper           = Vector3_Allocate ( ) ;
        lower            = Vector3_Allocate ( ) ;
        upper            = Vector3_Allocate ( ) ;
        isOK             = ( displacement     != NULL ) &&
                           ( iLower           != NULL ) &&
                           ( iUpper           != NULL ) &&
                           ( iTransformation3 != NULL ) &&
                           ( lower            != NULL ) &&
                           ( self             != NULL ) &&
                           ( upper            != NULL ) ;
        if ( isOK )
        {
            auto Boolean  doSkip, doTSkip ;
            auto Integer  i, t, tInverse ;
            auto Integer  a, aHigh, aInverse, aLow, b, bHigh, bInverse, bLow, c, cHigh, cInverse, cLow ;
	    auto Real     defaultScale, scale, u, v ;
            /* . Initialization. */
            self->cutOff = cutOff ;
            /* . Find the orthorhombic box within which interactions are to be sought. */
            Coordinates3_EnclosingOrthorhombicBox ( coordinates3, NULL, NULL, lower, upper ) ;
            Vector3_Add       ( upper, 1.0e+00, lower, NULL ) ; /* . To get absolute upper coordinates rather than extents. */
            Vector3_Increment ( lower, -cutOff ) ;
            Vector3_Increment ( upper,  cutOff ) ;
            /* . Loop over the transformations. */
            for ( t = 0 ; t < transformations->capacity ; t++ )
            {
                defaultScale = 0.5e+00 ;
                doTSkip      = False  ;
                tInverse     = transformations->inverses[t] ;
                if ( checkForInverses )
                {
                    defaultScale = 1.0e+00 ;          
                    doTSkip      = ( t < tInverse ) ;
                }
                /* . Generate the image transformation in real space. */
                Transformation3_CopyTo        ( transformations->items[t], iTransformation3 ) ;
                Transformation3_Orthogonalize ( iTransformation3, symmetryParameters->H, symmetryParameters->inverseH ) ;
                /* . Generate the coordinates. */
                Vector3_CopyTo ( lower, iLower, NULL ) ; Transformation3_ApplyToVector3 ( iTransformation3, iLower ) ;
                Vector3_CopyTo ( upper, iUpper, NULL ) ; Transformation3_ApplyToVector3 ( iTransformation3, iUpper ) ;
                for ( i = 0 ; i < 3 ; i++ )
                {
                    u = Vector3_Item ( iLower, i ) ; v = Vector3_Item ( iUpper, i ) ;
                    if ( u > v ) { Vector3_Item ( iLower, i ) = v ; Vector3_Item ( iUpper, i ) = u ; }
                } 
                /* . Find the limits for the search. */
                SymmetryParameters_FindBoxSearchLimits ( symmetryParameters, lower, upper, iLower, iUpper, &aLow, &aHigh, &bLow, &bHigh, &cLow, &cHigh ) ;
                /* . Expand the search range (for testing only). */
                if ( expandFactor > 0 )
                {
                    aLow -= expandFactor ; aHigh += expandFactor ;
                    bLow -= expandFactor ; bHigh += expandFactor ;
                    cLow -= expandFactor ; cHigh += expandFactor ;
                }
                /* . Loop over the search directions. */
	        for ( a = aLow ; a <= aHigh ; a++ )
	        {
	            for ( b = bLow ; b <= bHigh ; b++ )
	            {
        	        for ( c = cLow ; c <= cHigh ; c++ )
		        {
		            /* . Exclude self-interaction. */
		            if ( ( a == 0 ) && ( b == 0 ) && ( c == 0 ) && ( transformations->identity == t ) ) continue ;
                            /* . Set default values. */
                            doSkip = doTSkip      ;
                            scale  = defaultScale ;
                            /* . Treat self-inverses. */
                            if ( checkForInverses && ( tInverse == t ) )
                            {
                                /* . Find the inverse translation. */
                                Transformation3Container_FindInverseIntegerTranslation ( transformations, t, a, b, c, displacement, &aInverse, &bInverse, &cInverse ) ;
                                /* . The inverse is in the search range. */
                                /* . Note that this is only valid for the self-inverse (otherwise the other transformation's search range must be checked). */
                                if ( ( aInverse >= aLow ) && ( aInverse <= aHigh ) && ( bInverse >= bLow ) && ( bInverse <= bHigh ) && ( cInverse >= cLow ) && ( cInverse <= cHigh ) )
                                {
                                    /* . Pure self-inverses only occur once. */
                                    if ( ( a == aInverse ) && ( b == bInverse ) && ( c == cInverse ) ) scale = 0.5e+00 ;
                                    /* . Skip this image if its inverse will occur later. */
                                    else
                                    {
                                        doSkip = ( ( a < aInverse ) || ( ( a == aInverse ) && ( b < bInverse ) ) || ( ( a == aInverse ) && ( b == bInverse ) && ( c < cInverse ) ) ) ;
                                        scale  = 1.0e+00 ;
                                    }
                                }
                                /* . The inverse is out of the search range but treat it implicitly anyway (this should never occur). */
                                else scale = 1.0e+00 ;
                            }
                            /* . Displace the coordinates. */
                            SymmetryParameters_Displacement ( symmetryParameters, a, b, c, displacement ) ;
                            Vector3_Add ( iLower, 1.0e+00, displacement, NULL ) ;
                            Vector3_Add ( iUpper, 1.0e+00, displacement, NULL ) ;
                            /* . Do a final box check to skip boxes with no overlap. */
                            if ( ( Vector3_Item( iLower, 0 ) <= Vector3_Item( upper, 0 ) ) &&
                                 ( Vector3_Item( iLower, 1 ) <= Vector3_Item( upper, 1 ) ) &&
                                 ( Vector3_Item( iLower, 2 ) <= Vector3_Item( upper, 2 ) ) &&
                                 ( Vector3_Item( iUpper, 0 ) >= Vector3_Item( lower, 0 ) ) &&
                                 ( Vector3_Item( iUpper, 1 ) >= Vector3_Item( lower, 1 ) ) &&
                                 ( Vector3_Item( iUpper, 2 ) >= Vector3_Item( lower, 2 ) ) )
			    {
                                ImageScanContainer_Append ( self, doSkip, t, a, b, c, scale, status ) ;
                                if ( ! Status_IsOK ( status ) ) goto FinishUp ;
			    }
                            /* . Move back the coordinates. */
                            Vector3_Scale ( displacement, -1.0e+00 ) ;
                            Vector3_Add   ( iLower, 1.0e+00, displacement, NULL ) ;
                            Vector3_Add   ( iUpper, 1.0e+00, displacement, NULL ) ;
        	        }
                    }
	        }
	    }
        }
        else Status_Set ( status, Status_OutOfMemory ) ;
    FinishUp:
        Transformation3_Deallocate ( &iTransformation3 ) ;
        Vector3_Deallocate         ( &displacement     ) ;
        Vector3_Deallocate         ( &iLower           ) ;
        Vector3_Deallocate         ( &iUpper           ) ;
        Vector3_Deallocate         ( &lower            ) ;
        Vector3_Deallocate         ( &upper            ) ;
        if ( ! Status_IsOK ( status ) ) ImageScanContainer_Deallocate ( &self ) ;
    }
    return self ;
}


/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void ImageScanContainer_Deallocate ( ImageScanContainer **self )
{
    if ( (*self) != NULL )
    {
        Memory_Deallocate ( (*self)->records ) ;
        Memory_Deallocate ( (*self)          ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Reallocate records.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean ImageScanContainer_Reallocate ( ImageScanContainer *self, const Integer  capacity, Status *status )
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
            auto ImageScan *records  ;
            records = Memory_ReallocateArrayOfTypes ( self->records, n, ImageScan ) ;
            if ( records != NULL )
            {
                self->capacity = n ;
                self->records  = records ;
            }
            else isOK = False ;
        }
    }
    return isOK ;
}
