/*==================================================================================================================================
! . Symmetry parameter functions.
!=================================================================================================================================*/

# include <math.h>

# include "Memory.h"
# include "NumericalMacros.h"
# include "SymmetryParameters.h"
# include "Units.h"

/* . Notes:

   r = real space coordinates
   f = fractional coordinates

   r = H * f
   f = inverseH * r

   primary image - all f are in the range [0,1)

   minimum image vectors - all f are in the range [-1/2,1/2).

   isOrthogonal - alpha = beta = gamma = 90.0
                  this also implies that H and inverseH are diagonal (and thus symmetric)

*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
SymmetryParameters *SymmetryParameters_Allocate ( Status *status )
{
    SymmetryParameters *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( SymmetryParameters ) ;
        if ( self != NULL )
        {
            self->a            = 0.0e+00 ;
            self->b            = 0.0e+00 ;
            self->c            = 0.0e+00 ;
            self->alpha        = 0.0e+00 ;
            self->beta         = 0.0e+00 ;
            self->gamma        = 0.0e+00 ;
            self->H            = NULL    ;
            self->inverseH     = NULL    ;
            self->isOrthogonal = False   ;
        }
        else Status_Set ( status, Status_OutOfMemory ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation with matrices.
!---------------------------------------------------------------------------------------------------------------------------------*/
SymmetryParameters *SymmetryParameters_AllocateFull ( Status *status )
{
    SymmetryParameters *self = SymmetryParameters_Allocate ( status ) ;
    if ( self != NULL )
    {
        self->H        = Matrix33_Allocate ( ) ; Matrix33_Set ( self->H       , 0.0e+00 ) ;
        self->inverseH = Matrix33_Allocate ( ) ; Matrix33_Set ( self->inverseH, 0.0e+00 ) ;
        if ( ( self->H == NULL ) || ( self->inverseH == NULL ) )
        {
            SymmetryParameters_Deallocate ( &self ) ;
            Status_Set ( status, Status_OutOfMemory ) ;
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation given H and inverseH.
!---------------------------------------------------------------------------------------------------------------------------------*/
SymmetryParameters *SymmetryParameters_AllocateWithMatrices ( const Matrix33 *H        ,
                                                              const Matrix33 *inverseH ,
                                                                    Status   *status   )
{
    SymmetryParameters *self = NULL ;
    if ( ( H != NULL ) && ( inverseH != NULL ) && Status_IsOK ( status ) )
    {
        self = SymmetryParameters_Allocate ( status ) ;
        if ( self != NULL )
        {
            self->H        = Matrix33_CloneShallow (        H, status ) ; Matrix33_Set ( self->H       , 0.0e+00 ) ;
            self->inverseH = Matrix33_CloneShallow ( inverseH, status ) ; Matrix33_Set ( self->inverseH, 0.0e+00 ) ;
        }
        if ( ! Status_IsOK ( status ) ) SymmetryParameters_Deallocate ( &self ) ;
        if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Center coordinates within the primary image by free isolate.
! . By default all isolates are free.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameters_CenterCoordinates3ByFreeIsolate ( const SymmetryParameters  *self         ,
                                                          const SelectionContainer  *isolates     ,
                                                          const BooleanBlock        *freeIsolates ,
                                                                Coordinates3        *coordinates3 ,
                                                                Status              *status       )
{
    if ( ( self != NULL ) && ( isolates != NULL ) && ( coordinates3 != NULL ) && Status_IsOK ( status ) )
    {
        auto Vector3 *center, *translation ;
        center      = Vector3_Allocate ( ) ;
        translation = Vector3_Allocate ( ) ;
        if ( ( center != NULL ) && ( translation != NULL ) )
        {
            auto Boolean    checkWhetherFree = ( freeIsolates != NULL ) ;
            auto Integer    i ;
            auto Selection *isolate ;
            for ( i = 0 ; i < isolates->capacity ; i++ )
            {
                if ( ( ! checkWhetherFree ) || ( checkWhetherFree && Block_Item ( freeIsolates, i ) ) )
                {
                    isolate = isolates->items[i] ;
                    Coordinates3_Center ( coordinates3, isolate, NULL, &center ) ;
                    SymmetryParameters_FindCenteringTranslation ( self, center, translation ) ;
                    Coordinates3_Translate ( coordinates3, translation, isolate ) ;
                }
            }
        }
        else Status_Set ( status, Status_OutOfMemory ) ;
        Vector3_Deallocate ( &center      ) ;
        Vector3_Deallocate ( &translation ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Center coordinates within the primary image by index.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameters_CenterCoordinates3ByIndex ( const SymmetryParameters *self         ,
                                                    const Selection          *selection    ,
                                                          Coordinates3       *coordinates3 ,
                                                          Status             *status       )
{
    if ( ( self != NULL ) && ( coordinates3 != NULL ) && Status_IsOK ( status ) )
    {
        auto Coordinates3 *fractional ;
        auto Integer       n          ;
        /* . Allocate space. */
        if ( selection == NULL ) n = Coordinates3_Rows ( coordinates3 ) ;
        else                     n = selection->capacity ;
        fractional = Coordinates3_Allocate ( n, status ) ;
        if ( fractional != NULL )
        {
            auto Real x, y, z ;
            auto Integer    i ;
            /* . Get the fractional coordinates. */
            Coordinates3_Gather ( fractional, coordinates3  , selection ) ;
            Coordinates3_Rotate ( fractional, self->inverseH, NULL      ) ;
            /* . Shift the coordinates. */
            for ( i = 0 ; i < Coordinates3_Rows ( fractional ) ; i++ )
            {
                Coordinates3_GetRow ( fractional, i, x, y, z ) ;
                x -= floor ( x ) ;
                y -= floor ( y ) ;
                z -= floor ( z ) ;
                Coordinates3_SetRow ( fractional, i, x, y, z ) ;
            }
            /* . Back transform and copy back. */
            Coordinates3_Rotate  ( fractional, self->H     , NULL      ) ;
            Coordinates3_Scatter ( fractional, coordinates3, selection ) ;
            /* . Finish up. */
            Coordinates3_Deallocate ( &fractional ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Center coordinates within the primary image by isolate.
! . Only isolates with selected members are centered.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameters_CenterCoordinates3ByIsolate ( const SymmetryParameters *self         ,
                                                      const SelectionContainer *isolates     ,
                                                            Selection          *selection    ,
                                                            Coordinates3       *coordinates3 ,
                                                            Status             *status       )
{
    if ( ( self != NULL ) && ( isolates != NULL ) && ( coordinates3 != NULL ) && Status_IsOK ( status ) )
    {
        auto Boolean     *flags = NULL, hasSelection, isOK = True ;
        auto Integer      i, s ;
        auto RealArray1D *center, *translation ;
        auto Selection   *isolate ;
        /* . Allocate space. */
        center      = Vector3_Allocate ( ) ;
        translation = Vector3_Allocate ( ) ;
        /* . Check for a selection. */
        hasSelection = ( selection != NULL ) ;
        if ( hasSelection )
        {
            auto BooleanBlock *localFlags = Selection_MakeFlags ( selection, SelectionContainer_UpperBound ( isolates ), status ) ;
            if ( localFlags != NULL ) flags = Block_Items ( localFlags ) ;
        }
        /* . Check memory. */
        if ( ( center != NULL ) && ( translation != NULL ) && ( ( ! hasSelection ) || ( hasSelection && ( flags != NULL ) ) ) )
        {
            /* . Loop over isolates. */
            for ( i = 0 ; i < isolates->capacity ; i++, isOK = True )
            {
                isolate = isolates->items[i] ;
                /* . Check whether to include the isolate. */
                if ( hasSelection  )
                {
                    for ( s = 0 ; s < isolate->capacity ; s++ )
                    {
                        if ( ! flags[isolate->indices[s]] ) { isOK = False ; break ; }
                    }
                }
                /* . Do the centering. */
                if ( isOK )
                {
                    Coordinates3_Center ( coordinates3, isolate, NULL, &center ) ;
                    SymmetryParameters_FindCenteringTranslation ( self, center, translation ) ;
                    Coordinates3_Translate ( coordinates3, translation, isolate ) ;
                }
            }
        }
        else Status_Set ( status, Status_OutOfMemory ) ;
        /* . Finish up. */
        Vector3_Deallocate ( &center      ) ;
        Vector3_Deallocate ( &translation ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Clear the M/inverseH representation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameters_ClearH ( SymmetryParameters *self )
{
    if ( self != NULL )
    {
        Matrix33_Deallocate ( &(self->H       ) ) ;
        Matrix33_Deallocate ( &(self->inverseH) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Transfer data from one structure to another.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameters_CopyTo ( const SymmetryParameters *self, SymmetryParameters *other )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        other->a            = self->a            ;
        other->b            = self->b            ;
        other->c            = self->c            ;
        other->alpha        = self->alpha        ;
        other->beta         = self->beta         ;
        other->gamma        = self->gamma        ;
        other->isOrthogonal = self->isOrthogonal ;
        SymmetryParameters_MakeH ( other ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameters_Deallocate ( SymmetryParameters **self )
{
    if ( (*self) != NULL )
    {
        SymmetryParameters_ClearH ( (*self) ) ;
        Memory_Deallocate         ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate a displacement (given in terms of a, b, c which are the columns
! . of the matrix M).
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameters_Displacement ( const SymmetryParameters *self         ,
                                       const Integer             a            ,
                                       const Integer             b            ,
                                       const Integer             c            ,
                                             Vector3            *displacement )
{
    if ( ( self != NULL ) && ( displacement != NULL ) )
    {
        auto Integer i ;
        auto Real    nA, nB, nC ;
        nA = ( Real ) a ; nB = ( Real ) b ; nC = ( Real ) c ;
        for ( i = 0 ; i < 3 ; i++ ) { Vector3_Item ( displacement, i ) = ( nA * Matrix33_Item ( self->H, i, 0 ) +
                                                                           nB * Matrix33_Item ( self->H, i, 1 ) +
                                                                           nC * Matrix33_Item ( self->H, i, 2 ) ) ; }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find the range of a, b and c values for which an image box overlaps with a central box.
! . This procedure makes use of the special properties of H and the boxes are assumed to be orthorhombic.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . This needs to be generalized to other orientations for H. */
/* . Limits method. */
static Boolean GetLimits ( const Real bl, const Real bu, Real il, Real iu, const Real t, Integer *low, Integer *high )
{
    Boolean QOK = False ;
    Integer  n   = 0     ;
    /* . Move left until iu < bl. */
    while ( iu >= bl ) { il -= t ; iu -= t ; n-- ; }
    /* . Move right until iu >= bl. */
    while ( iu <  bl ) { il += t ; iu += t ; n++ ; }
    /* . Check that there is some overlap. */
    if ( il <= bu ) /* . && iu >= bl. */
    {
        (*low)  = n ;
        /* . Move right until il > bu. */
        while ( il <= bu ) { il += t ; iu += t ; n++ ; }
        (*high) = n - 1 ;
        /* . Everything is OK. */
        QOK = True ;
    }
    return QOK ;
}

/* . Box method. */
void SymmetryParameters_FindBoxSearchLimits ( const SymmetryParameters *self   ,
                                              const Vector3            *lower  ,
                                              const Vector3            *upper  ,
                                              const Vector3            *iLower ,
                                              const Vector3            *iUpper ,
                                                    Integer            *aLow   ,
                                                    Integer            *aHigh  ,
                                                    Integer            *bLow   ,
                                                    Integer            *bHigh  ,
                                                    Integer            *cLow   ,
                                                    Integer            *cHigh  )
{
    /* . Initialization. */
    (*aLow)  = (*bLow)  = (*cLow)  =  0 ;
    (*aHigh) = (*bHigh) = (*cHigh) = -1 ;
    if ( ( self != NULL ) && ( lower != NULL ) && ( upper != NULL ) && ( iLower != NULL ) && ( iUpper != NULL ) )
    {
        auto Boolean   QOK ;
        auto Real bl, bu, d1, d2 ;
        /* . Do c first as this is the only lattice vector that contributes to z. */
        QOK = GetLimits ( lower->data[2], upper->data[2], iLower->data[2], iUpper->data[2], Matrix33_Item ( self->H, 2, 2 ), cLow, cHigh ) ;
        /* . Now do b which contributes to y (along with c). */
        if ( QOK )
        {
            d1 = (*cLow)  * Matrix33_Item ( self->H, 1, 2 ) ;
            d2 = (*cHigh) * Matrix33_Item ( self->H, 1, 2 ) ;
            bl = lower->data[1] - Maximum ( d1, d2 ) ;
            bu = upper->data[1] - Minimum ( d1, d2 ) ;
            QOK = GetLimits ( bl, bu, iLower->data[1], iUpper->data[1], Matrix33_Item ( self->H, 1, 1 ), bLow, bHigh ) ;
            /* . Now do a which contributes to x (along with b and c). */
            if ( QOK )
            {
                d1  = (*cLow ) * Matrix33_Item ( self->H, 0, 2 ) ;
                d2  = (*cHigh) * Matrix33_Item ( self->H, 0, 2 ) ;
                bl  = lower->data[0] - Maximum ( d1, d2 ) ;
                bu  = upper->data[0] - Minimum ( d1, d2 ) ;
                d1  = (*bLow ) * Matrix33_Item ( self->H, 0, 1 ) ;
                d2  = (*bHigh) * Matrix33_Item ( self->H, 0, 1 ) ;
                bl -= Maximum ( d1, d2 ) ;
                bu -= Minimum ( d1, d2 ) ;
                QOK = GetLimits ( bl, bu, iLower->data[0], iUpper->data[0], Matrix33_Item ( self->H, 0, 0 ), aLow, aHigh ) ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find the translation that puts the real space point inside the primary image.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameters_FindCenteringTranslation ( const SymmetryParameters *self        ,
                                                   const Vector3            *point       ,
                                                         Vector3            *translation )
{
    if ( ( self != NULL ) && ( point != NULL ) && ( translation != NULL ) )
    {
        auto Integer a, b, c ;
        Vector3_CopyTo ( point, translation, NULL ) ;
        Matrix33_ApplyToVector3 ( self->inverseH, translation ) ;
        a = ( Integer ) floor ( Vector3_Item ( translation, 0 ) ) ;
        b = ( Integer ) floor ( Vector3_Item ( translation, 1 ) ) ;
        c = ( Integer ) floor ( Vector3_Item ( translation, 2 ) ) ;
        SymmetryParameters_Displacement ( self, -a, -b, -c, translation ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Check to see if the minimum image convention is satisfied given a length.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean SymmetryParameters_IsMinimumImageConventionSatisfied ( const SymmetryParameters *self, const Real length )
{
    Boolean isOK = False ;
    if ( self != NULL )
    {
        auto Real d = 2.0e+00 * length, widths[3] ;
        SymmetryParameters_PerpendicularWidths ( self, widths ) ;
        isOK = ( d <= widths[0] ) && ( d <= widths[1] ) && ( d <= widths[2] ) ;
    }
    return isOK ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Orthogonality of unit cell.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean SymmetryParameters_IsOrthogonal ( const SymmetryParameters *self ) { return ( ( self == NULL ) ? False : self->isOrthogonal ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Scale the symmetry parameters isotropically.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameters_IsotropicScale ( SymmetryParameters *self, const Real scale )
{
    if ( self != NULL )
    {
        self->a *= scale ;
        self->b *= scale ;
        self->c *= scale ;
        Matrix33_Scale ( self->H, scale ) ;
        if ( scale == 0.0e+00 ) Matrix33_Set   ( self->inverseH, 0.0e+00         ) ;
        else                    Matrix33_Scale ( self->inverseH, 1.0e+00 / scale ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make fractional from real coordinates ensuring that the result is in the primary image.
!---------------------------------------------------------------------------------------------------------------------------------*/
Coordinates3 *SymmetryParameters_MakeFractionalCoordinates ( const SymmetryParameters *self         ,
                                                             const Coordinates3       *coordinates3 ,
                                                                   Status             *status       )
{
    Coordinates3 *fractional = NULL ;
    if ( ( self != NULL ) && ( coordinates3 != NULL ) && Status_IsOK ( status ) )
    {
        fractional = Coordinates3_Allocate ( Coordinates3_Rows ( coordinates3 ), status ) ;
        if ( fractional != NULL )
        {
            auto Integer c, r ;
            auto Real    v ;
            RealArray2D_MatrixMultiply ( False, True, 1.0e+00, coordinates3, self->inverseH, 0.0e+00, fractional, status ) ;
            for ( r = 0 ; r < Coordinates3_Rows ( fractional ) ; r++ )
            {
                for ( c = 0 ; c < 3 ; c++ )
                {
                    v = floor ( Coordinates3_Item ( fractional, r, c ) ) ;
                    Coordinates3_Item ( fractional, r, c ) -= v ;
                }
            }
        }
    }
    return fractional ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the metric matrix G = H^T * H.
! . This is independent of the orientation used for constructing H.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameters_MakeG ( const SymmetryParameters *self, Matrix33 *G )
{
    if ( ( self != NULL ) && ( G != NULL ) )
    {
        auto Real a, abG, acB, b, bcA, c ;
        a   = self->a ;
        b   = self->b ;
        c   = self->c ;
        abG = a * b * cos ( self->gamma * Units_Angle_Degrees_To_Radians ) ;
        acB = a * c * cos ( self->beta  * Units_Angle_Degrees_To_Radians ) ;
        bcA = b * c * cos ( self->alpha * Units_Angle_Degrees_To_Radians ) ;
        Matrix33_Set  ( G, 0.0e+00 ) ;
        Matrix33_Item ( G, 0, 0 ) = a * a ; Matrix33_Item ( G, 0, 1 ) = abG   ; Matrix33_Item ( G, 0, 2 ) = acB   ;
        Matrix33_Item ( G, 1, 0 ) = abG   ; Matrix33_Item ( G, 1, 1 ) = b * b ; Matrix33_Item ( G, 1, 2 ) = bcA   ;
        Matrix33_Item ( G, 2, 0 ) = acB   ; Matrix33_Item ( G, 2, 1 ) = bcA   ; Matrix33_Item ( G, 2, 2 ) = c * c ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the H and inverseH matrices.
! . The columns of the matrix H correspond to the lattice vectors a, b, c.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameters_MakeH ( SymmetryParameters *self )
{
    if ( self != NULL )
    {
        auto Real alpha, beta, cosAlpha, cosBeta, cosGamma, gamma, sinGamma ;
        /* . Some factors. */
        alpha    = self->alpha * Units_Angle_Degrees_To_Radians ;
        beta     = self->beta  * Units_Angle_Degrees_To_Radians ;
        gamma    = self->gamma * Units_Angle_Degrees_To_Radians ;
        cosAlpha = cos ( alpha ) ;
        cosBeta  = cos ( beta  ) ;
        cosGamma = cos ( gamma ) ;
        sinGamma = sin ( gamma ) ;
        /* . Create the H matrix - standard orientation. */
        Matrix33_Set  ( self->H, 0.0e+00 ) ;
        Matrix33_Item ( self->H, 0, 0 ) = self->a ;
        Matrix33_Item ( self->H, 0, 1 ) = self->b * cosGamma ;
        Matrix33_Item ( self->H, 1, 1 ) = self->b * sinGamma ;
        Matrix33_Item ( self->H, 0, 2 ) = self->c * cosBeta  ;
        Matrix33_Item ( self->H, 1, 2 ) = self->c * ( cosAlpha - cosBeta * cosGamma ) / sinGamma ;
        Matrix33_Item ( self->H, 2, 2 ) = self->c * sqrt ( 1.0e+00 - cosAlpha * cosAlpha -
                                                                     cosBeta  * cosBeta  -
                                                                     cosGamma * cosGamma +
                                                           2.0e+00 * cosAlpha * cosBeta * cosGamma ) / sinGamma ;
        /* . Invert it. */
        Matrix33_Invert ( self->inverseH, self->H ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Apply the minimum image convention to an interaction vector.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameters_MakeMinimumImageVector ( const SymmetryParameters *self, Real *r, Real *dR )
{
    if ( ( self != NULL ) && ( r != NULL ) )
    {
        auto Boolean doDR = ( dR != NULL ) ;
        auto Integer i ;
        if ( self->isOrthogonal )
        {
            auto Real d ;
            for ( i = 0 ; i < 3 ; i++ )
            {
                d     = - Matrix33_Item ( self->H, i, i ) * Real_RoundToInteger ( Matrix33_Item ( self->inverseH, i, i ) * r[i] ) ;
                r[i] += d ; 
                if ( doDR ) dR[i] = d ;
            }
        }
        else
        {
            auto Real d, t[3] ;
            for ( i = 0 ; i < 3 ; i++ )
            {
                d    = Matrix33_Item ( self->inverseH, i, 0 ) * r[0] +
                       Matrix33_Item ( self->inverseH, i, 1 ) * r[1] +
                       Matrix33_Item ( self->inverseH, i, 2 ) * r[2] ;
                t[i] = - Real_RoundToInteger ( d ) ;
            }
            for ( i = 0 ; i < 3 ; i++ )
            {
                d     = Matrix33_Item ( self->H, i, 0 ) * t[0] +
                        Matrix33_Item ( self->H, i, 1 ) * t[1] +
                        Matrix33_Item ( self->H, i, 2 ) * t[2] ;
                r[i] += d ;
                if ( doDR ) dR[i] = d ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Perpendicular widths.
!---------------------------------------------------------------------------------------------------------------------------------*/
void SymmetryParameters_PerpendicularWidths ( const SymmetryParameters *self, Real *widths )
{
    if ( ( self != NULL ) && ( widths != NULL ) )
    {
        if ( self->isOrthogonal )
        {
           widths[0] = Matrix33_Item ( self->H, 0, 0 ) ;
           widths[1] = Matrix33_Item ( self->H, 1, 1 ) ;
           widths[2] = Matrix33_Item ( self->H, 2, 2 ) ;
        }
        else
        {
            auto Integer  i, j, x, y ;
            auto Vector3 *u, *v ;
            u = Vector3_Allocate ( ) ;
            v = Vector3_Allocate ( ) ;
            for ( i = 0 ; i < 3 ; i++ )
            {
                x = ( i + 1 ) % 3 ;
                y = ( i + 2 ) % 3 ;
                for ( j = 0 ; j < 3 ; j++ )
                {
                    Vector3_Item ( u, j ) = Matrix33_Item ( self->H, j, x ) ;
                    Vector3_Item ( v, j ) = Matrix33_Item ( self->H, j, y ) ;
                }
                Vector3_CrossProduct ( u, v ) ;
                for ( j = 0 ; j < 3 ; j++ ) { Vector3_Item ( v, j ) = Matrix33_Item ( self->H, j, i ) ; }
                widths[i] = fabs ( Vector3_Dot ( u, v, NULL ) ) / Vector3_Norm2 ( u ) ;
            }
            Vector3_Deallocate ( &u ) ;
            Vector3_Deallocate ( &v ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set the parameters appropriate for a crystal.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _Ninety                 90.0e+00
# define _OrthogonalityTolerance 1.0e-04
void SymmetryParameters_SetCrystalParameters (       SymmetryParameters *self  ,
                                               const Real                a     ,
                                               const Real                b     ,
                                               const Real                c     ,
                                               const Real                alpha ,
                                               const Real                beta  ,
                                               const Real                gamma )
{
    if ( self != NULL )
    {
        self->a            = a     ;
        self->b            = b     ;
        self->c            = c     ;
        self->alpha        = alpha ;
        self->beta         = beta  ;
        self->gamma        = gamma ;
        self->isOrthogonal = ( fabs ( alpha - _Ninety ) <= _OrthogonalityTolerance ) &&
                             ( fabs ( beta  - _Ninety ) <= _OrthogonalityTolerance ) &&
                             ( fabs ( gamma - _Ninety ) <= _OrthogonalityTolerance ) ;
        SymmetryParameters_MakeH ( self ) ;
    }
}
# undef _Ninety
# undef _OrthogonalityTolerance

/*----------------------------------------------------------------------------------------------------------------------------------
! . Volume.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real SymmetryParameters_Volume ( const SymmetryParameters *self ) { return ( ( self == NULL ) ? 0.0e+00 : Matrix33_Determinant ( self->H ) ) ; }
