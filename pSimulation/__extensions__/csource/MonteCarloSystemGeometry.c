/*==================================================================================================================================
! . This module implements Monte Carlo system geometry procedures.
!=================================================================================================================================*/

# include <math.h>

# include "Memory.h"
# include "MonteCarloSystemGeometry.h"
# include "NumericalMacros.h"
# include "Units.h"

/*==================================================================================================================================
! . Private procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Check to see whether a move has been rejected.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _ExponentialUnderFlow 75.0e+00
static Boolean IsMoveRejected ( Real deltaEb, Real random )
{
    auto Boolean isAccepted = False ;
    if ( deltaEb < _ExponentialUnderFlow )
    {
        isAccepted = ( deltaEb <= 0.0e+00 ) || ( exp ( - deltaEb ) > random ) ;
    }
    return ! isAccepted ;
}
# undef _ExponentialUnderFlow

/*==================================================================================================================================
! . Public procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
MonteCarloSystemGeometry *MonteCarloSystemGeometry_Allocate ( const Integer  numberOfParticles ,
                                                              const Integer  numberOfRandom    )
{
    MonteCarloSystemGeometry *self = NULL ;
    if ( ( numberOfParticles > 0 ) && ( numberOfRandom > 0 ) )
    {
        self = Memory_AllocateType ( MonteCarloSystemGeometry ) ;
        if ( self != NULL )
        {
            auto Boolean isOK ;
            /* . Counters. */
            self->blocks             = 0 ;
            self->moves              = 0 ;
            self->nReject            = 0 ;
            self->nRejectM           = 0 ;
            self->nRejectT           = 0 ;
            self->nRejectV           = 0 ;
            self->nTryM              = 0 ;
            self->nTryV              = 0 ;
            /* . Current values and other factors. */
            self->beta               = 0.0e+00 ;
            self->dielectric         = 1.0e+00 ;
            self->eCurrent           = 0.0e+00 ;
            self->pressure           = 0.0e+00 ;
            self->tFactor            = 0.0e+00 ;
            self->volume             = 0.0e+00 ;
            /* . Move sizes. */
            self->acceptanceRatio    = 0.0e+00 ;
            self->rMax               = 0.0e+00 ;
            self->tMax               = 0.0e+00 ;
            self->vMax               = 0.0e+00 ;
            /* . Statistics. */
            self->eAv                = 0.0e+00 ;
            self->eAv2               = 0.0e+00 ;
            self->eTot               = 0.0e+00 ;
            self->eTot2              = 0.0e+00 ;
            self->eTotB              = 0.0e+00 ;
            self->eTotB2             = 0.0e+00 ;
            self->hAv                = 0.0e+00 ;
            self->hAv2               = 0.0e+00 ;
            self->hTot               = 0.0e+00 ;
            self->hTot2              = 0.0e+00 ;
            self->hTotB              = 0.0e+00 ;
            self->hTotB2             = 0.0e+00 ;
            self->vAv                = 0.0e+00 ;
            self->vAv2               = 0.0e+00 ;
            self->vTot               = 0.0e+00 ;
            self->vTot2              = 0.0e+00 ;
            self->vTotB              = 0.0e+00 ;
            self->vTotB2             = 0.0e+00 ;
            /* . Aliases. */
            self->charges               = NULL ;
            self->coordinates3          = NULL ;
            self->isolates              = NULL ;
            self->ljParameters          = NULL ;
            self->ljTypes               = NULL ;
            self->pairwiseInteraction   = NULL ;
            self->symmetryParameters    = NULL ;
            /* . Arrays to allocate. */
            self->random                = Memory_AllocateArrayOfTypes ( numberOfRandom, Real ) ;
            self->oldCoordinates3       = Coordinates3_Allocate ( numberOfParticles, NULL ) ;
            self->rotation              = Matrix33_Allocate               (      ) ;
            self->oldSymmetryParameters = SymmetryParameters_AllocateFull ( NULL ) ;
            self->translation           = Vector3_Allocate                (      ) ;
            /* . Deallocate if there is not enough memory. */
            isOK = ( self->random                != NULL ) &&
                   ( self->oldCoordinates3       != NULL ) &&
                   ( self->rotation              != NULL ) &&
                   ( self->oldSymmetryParameters != NULL ) &&
                   ( self->translation           != NULL ) ;
            if ( ! isOK ) MonteCarloSystemGeometry_Deallocate ( &self ) ;
            else Memory_Set ( self->random, numberOfRandom, 0.0e+00 ) ;
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MonteCarloSystemGeometry_Deallocate ( MonteCarloSystemGeometry **self )
{
    if ( (*self) != NULL )
    {
        Coordinates3_Deallocate       ( &((*self)->oldCoordinates3      ) ) ;
        Matrix33_Deallocate           ( &((*self)->rotation             ) ) ;
        SymmetryParameters_Deallocate ( &((*self)->oldSymmetryParameters) ) ;
        Vector3_Deallocate            ( &((*self)->translation          ) ) ;
        Memory_Deallocate             (   (*self)->random                 ) ;
        Memory_Deallocate             (   (*self)                         ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Adjust the move sizes.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _Down 0.95
# define _Up   1.05
void MonteCarloSystemGeometry_AdjustMoveSizes ( MonteCarloSystemGeometry *self )
{
    if ( self != NULL )
    {
        /* . Adjust the rotation and translation move sizes. */
        if ( self->nTryM > 0 )
        {
            if ( ( ( Real ) ( self->nTryM - self->nRejectM ) / ( Real ) ( self->nTryM ) ) > self->acceptanceRatio )
            {
                self->rMax *=   _Up ;
                self->tMax *=   _Up ;
            }
            else
            {
                self->rMax *= _Down ;
                self->tMax *= _Down ;
            }
            self->nRejectM = 0 ; self->nTryM = 0 ;
        }
        /* . Adjust the volume move size. */
        if ( self->nTryV > 0 )
        {
            if ( ( ( Real ) ( self->nTryV - self->nRejectV ) / ( Real ) ( self->nTryV ) ) > self->acceptanceRatio )
            {
                self->vMax *=   _Up ;
            }
            else
            {
                self->vMax *= _Down ;
            }
            self->nRejectV = 0 ; self->nTryV = 0 ;
        }
    }
}
# undef _Down
# undef _Up

/*----------------------------------------------------------------------------------------------------------------------------------
! . Do an isolate move.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status MonteCarloSystemGeometry_MoveIsolate ( MonteCarloSystemGeometry *self )
{
    Status status = Status_OK ;
    if ( self != NULL )
    {
        auto Integer       axis, chosen, i, rIndex = 0 ;
        auto Real          angle, eAfter, eBefore, oldEnergy ;
        auto Coordinates3 *coordinates3 ;
        auto Matrix33     *rotation     ;
        auto Selection    *cSelection   ;
        auto Vector3      *translation  ;
        /* . Set some aliases. */
        coordinates3 = self->coordinates3 ;
        rotation     = self->rotation     ;
        translation  = self->translation  ;
        /* . Increment the number of tries. */
        self->nTryM += 1 ;
        /* . Choose a isolate to move. */
        chosen     = ( Integer ) floor ( ( Real ) self->isolates->capacity * self->random[rIndex] ) ; rIndex += 1 ;
        cSelection = self->isolates->items[chosen] ;
        /* . Calculate the energy of the isolate at the old configuration. */
        eBefore = PairwiseInteractionMonteCarlo_MMMMIsolateEnergy ( self->pairwiseInteraction  ,
                                                                    chosen                     ,
                                                                    self->charges              ,
                                                                    self->ljTypes              ,
                                                                    self->ljParameters         ,
                                                                    1.0e+00 / self->dielectric ,
                                                                    1.0e+00                    ,
                                                                    self->isolates             ,
                                                                    NULL                       ,
                                                                    coordinates3               ,
                                                                    self->symmetryParameters   ,
                                                                    NULL                       ,
                                                                    NULL                       ,
                                                                    &status                    ) ;
        /* . SAve the old configuration. */
        oldEnergy = self->eCurrent ;
        Coordinates3_Gather ( self->oldCoordinates3, coordinates3, cSelection ) ;
        /* . Calculate the center of the isolate and translate to the origin. */
        Coordinates3_Center    ( coordinates3, cSelection, NULL, &translation ) ;
        Vector3_Scale          ( translation, -1.0e+00 ) ;
        Coordinates3_Translate ( coordinates3, translation, cSelection ) ;
        /* . Do a rotation but only if the isolate has more than one particle. */
        if ( cSelection->capacity > 1 )
        {
            angle = 2.0e+00 * self->rMax * ( self->random[rIndex] - 0.5e+00 ) * Units_Angle_Degrees_To_Radians ; rIndex += 1 ;
            axis  = ( Integer ) floor ( 3.0e+00 * self->random[rIndex] ) ; rIndex += 1 ;
            switch ( axis )
            {
                case 0: Matrix33_RotationAboutAxis ( &rotation, angle, 1.0e+00, 0.0e+00, 0.0e+00 ) ; break ;
                case 1: Matrix33_RotationAboutAxis ( &rotation, angle, 0.0e+00, 1.0e+00, 0.0e+00 ) ; break ;
                case 2: Matrix33_RotationAboutAxis ( &rotation, angle, 0.0e+00, 0.0e+00, 1.0e+00 ) ; break ;
            }
            Coordinates3_Rotate ( coordinates3, rotation, cSelection ) ;
        }
        /* . Calculate the translation for the isolate within the minimum image convention. */
        Vector3_Scale ( translation, -1.0e+00 ) ;
        for ( i = 0 ; i < 3 ; i++, rIndex++ ) translation->data[i] += 2.0e+00 * self->tMax * ( self->random[rIndex] - 0.5e+00 ) ;
        SymmetryParameters_MakeMinimumImageVector ( self->symmetryParameters, Array1D_Data ( translation ), NULL ) ;
        Coordinates3_Translate ( coordinates3, translation, cSelection ) ;
        /* . Calculate the energy of the isolate at the new configuration. */
        eAfter = PairwiseInteractionMonteCarlo_MMMMIsolateEnergy ( self->pairwiseInteraction  ,
                                                                   chosen                     ,
                                                                   self->charges              ,
                                                                   self->ljTypes              ,
                                                                   self->ljParameters         ,
                                                                   1.0e+00 / self->dielectric ,
                                                                   1.0e+00                    ,
                                                                   self->isolates             ,
                                                                   NULL                       ,
                                                                   coordinates3               ,
                                                                   self->symmetryParameters   ,
                                                                   NULL                       ,
                                                                   NULL                       ,
                                                                   &status                    ) ;
        /* . Calculate the total energy of the new configuration. */
        self->eCurrent = oldEnergy + eAfter - eBefore ;
        /* . Check to see if the move is rejected. */
        if ( IsMoveRejected ( self->beta * ( self->eCurrent - oldEnergy ), self->random[rIndex] ) )
        {
            /* . Increment nReject and nRejectM. */
            self->nReject  += 1 ;
            self->nRejectM += 1 ;
            /* . Reactivate the old configuration. */
            self->eCurrent = oldEnergy ;
            Coordinates3_Scatter ( self->oldCoordinates3, coordinates3, cSelection ) ;
        }
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Do a volume move.
! . The volume is changed isotropically.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status MonteCarloSystemGeometry_MoveVolume ( MonteCarloSystemGeometry *self )
{
    Status status = Status_OK ;
    if ( self != NULL )
    {
        auto Integer             i, rIndex = 0 ;
        auto Real                delta, gamma, oldEnergy, newVolume, oldVolume ;
        auto Coordinates3       *coordinates3       ;
        auto Selection          *iSelection         ;
        auto SymmetryParameters *symmetryParameters ;
        auto Vector3            *translation        ;
        /* . Set some aliases. */
        coordinates3       = self->coordinates3       ;
        translation        = self->translation        ;
        symmetryParameters = self->symmetryParameters ;
        /* . Increment the number of tries. */
        self->nTryV += 1 ;
        /* . Save the old configuration. */
        oldEnergy = self->eCurrent ;
        oldVolume = self->volume   ;
        Coordinates3_CopyTo       ( coordinates3      , self->oldCoordinates3, NULL ) ;
        SymmetryParameters_CopyTo ( symmetryParameters, self->oldSymmetryParameters ) ;
        /* . Calculate the new volume and scale the symmetry parameters accordingly. */
        newVolume = oldVolume + 2.0e+00 * self->vMax * ( self->random[rIndex] - 0.5e+00 ) ; rIndex += 1 ;
        gamma     = pow ( ( newVolume / oldVolume ), 1.0e+00 / 3.0e+00 )   ;
        SymmetryParameters_IsotropicScale ( symmetryParameters, gamma )    ;
        self->volume = SymmetryParameters_Volume ( symmetryParameters ) ;
        /* . Check the minimum image convention. */
        if ( SymmetryParameters_IsMinimumImageConventionSatisfied ( symmetryParameters, self->pairwiseInteraction->cutOff ) )
        {
            /* . Translate the coordinates of the particles in each isolate. */
            gamma -= 1.0e+00 ;
            for ( i = 0 ; i < self->isolates->capacity ; i++ )
            {
                iSelection = self->isolates->items[i] ;
                Coordinates3_Center    ( coordinates3, iSelection, NULL, &translation ) ;
                Vector3_Scale          ( translation, gamma ) ;
                Coordinates3_Translate ( coordinates3, translation, iSelection ) ;
            }
            /* . Calculate the total energy for the configuration. */
            self->eCurrent = PairwiseInteractionMonteCarlo_MMMMEnergy ( self->pairwiseInteraction  ,
                                                                        self->charges              ,
                                                                        self->ljTypes              ,
                                                                        self->ljParameters         ,
                                                                        1.0e+00 / self->dielectric ,
                                                                        1.0e+00                    ,
                                                                        self->isolates             ,
                                                                        NULL                       ,
                                                                        coordinates3               ,
                                                                        symmetryParameters         ,
                                                                        NULL                       ,
                                                                        NULL                       ,
                                                                        &status                    ) ;
            /* . Check to see if the move is rejected. */
            delta = self->beta * ( self->eCurrent - oldEnergy + self->pressure * ( self->volume - oldVolume ) - self->tFactor * log ( self->volume / oldVolume ) ) ;
            if ( IsMoveRejected ( delta, self->random[rIndex] ) )
            {
                /* . Increment nReject and nRejectV. */
                self->nReject  += 1 ;
                self->nRejectV += 1 ;
                /* . Reactivate the old configuration. */
                self->eCurrent = oldEnergy ;
                self->volume   = oldVolume ;
                Coordinates3_CopyTo       ( self->oldCoordinates3      , coordinates3, NULL ) ;
                SymmetryParameters_CopyTo ( self->oldSymmetryParameters, symmetryParameters ) ;
            }
        }
        else status = Status_AlgorithmError ;
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Accumulate block statistics.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MonteCarloSystemGeometry_StatisticsBlockAccumulate ( MonteCarloSystemGeometry *self )
{
    if ( self != NULL )
    {
        auto Real e, h, v ;
        e = self->eCurrent ;
        v = self->volume   ;
        h = e + self->pressure * v ;
        self->eAv += e ; self->eAv2 += e * e ;
        self->hAv += h ; self->hAv2 += h * h ;
        self->vAv += v ; self->vAv2 += v * v ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Start block statistics.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MonteCarloSystemGeometry_StatisticsBlockStart ( MonteCarloSystemGeometry *self )
{
    if ( self != NULL )
    {
        self->nReject = 0 ;
        self->eAv = 0.0e+00 ; self->eAv2 = 0.0e+00 ;
        self->hAv = 0.0e+00 ; self->hAv2 = 0.0e+00 ;
        self->vAv = 0.0e+00 ; self->vAv2 = 0.0e+00 ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Stop block statistics.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MonteCarloSystemGeometry_StatisticsBlockStop ( MonteCarloSystemGeometry *self )
{
    if ( self != NULL )
    {
        auto Real e2, h2, n, v2 ;
        /* . Accumulate run statistics. */
        self->nRejectT += self->nReject ;
        self->eTot += self->eAv ; self->eTot2 += self->eAv2 ;
        self->hTot += self->hAv ; self->hTot2 += self->hAv2 ;
        self->vTot += self->vAv ; self->vTot2 += self->vAv2 ;
        /* . Calculate block statistics. */
        n = ( Real ) self->moves ;
        self->eAv /= n ; e2 = self->eAv2 / n - self->eAv * self->eAv ; self->eAv2 = Maximum ( e2, 0.0e+00 ) ;
        self->hAv /= n ; h2 = self->hAv2 / n - self->hAv * self->hAv ; self->hAv2 = Maximum ( h2, 0.0e+00 ) ;
        self->vAv /= n ; v2 = self->vAv2 / n - self->vAv * self->vAv ; self->vAv2 = Maximum ( v2, 0.0e+00 ) ;
        /* . Accumulate run block statistics. */
        self->eTotB += self->eAv ; self->eTotB2 += self->eAv * self->eAv ;
        self->hTotB += self->hAv ; self->hTotB2 += self->hAv * self->hAv ;
        self->vTotB += self->vAv ; self->vTotB2 += self->vAv * self->vAv ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Start statistics.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MonteCarloSystemGeometry_StatisticsStart ( MonteCarloSystemGeometry *self )
{
    if ( self != NULL )
    {
        self->nRejectT = 0 ;
        self->nRejectM = 0 ; self->nTryM = 0 ;
        self->nRejectV = 0 ; self->nTryV = 0 ;
        /* . Run. */
        self->eTot  = 0.0e+00 ; self->eTot2  = 0.0e+00 ;
        self->hTot  = 0.0e+00 ; self->hTot2  = 0.0e+00 ;
        self->vTot  = 0.0e+00 ; self->vTot2  = 0.0e+00 ;
        /* . Block. */
        self->eTotB = 0.0e+00 ; self->eTotB2 = 0.0e+00 ;
        self->hTotB = 0.0e+00 ; self->hTotB2 = 0.0e+00 ;
        self->vTotB = 0.0e+00 ; self->vTotB2 = 0.0e+00 ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Start statistics.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MonteCarloSystemGeometry_StatisticsStop ( MonteCarloSystemGeometry *self )
{
    if ( ( self != NULL ) && ( self->blocks > 1 ) )
    {
        auto Real e2, h2, n, v2 ;
        /* . Run. */
        n = ( Real ) ( self->blocks * self->moves ) ;
        self->eTot /= n ; e2 = self->eTot2 / n - self->eTot * self->eTot ; self->eTot2 = Maximum ( e2, 0.0e+00 ) ;
        self->hTot /= n ; h2 = self->hTot2 / n - self->hTot * self->hTot ; self->hTot2 = Maximum ( h2, 0.0e+00 ) ;
        self->vTot /= n ; v2 = self->vTot2 / n - self->vTot * self->vTot ; self->vTot2 = Maximum ( v2, 0.0e+00 ) ;
        /* . Block - see A&T page 192. */
        n = ( Real ) self->blocks ;
        self->eTotB /= n ; e2 = self->eTotB2 / n + self->eTot * ( self->eTot - 2.0e+00 * self->eTotB ) ; self->eTotB2 = Maximum ( e2, 0.0e+00 ) ;
        self->hTotB /= n ; h2 = self->hTotB2 / n + self->hTot * ( self->hTot - 2.0e+00 * self->hTotB ) ; self->hTotB2 = Maximum ( h2, 0.0e+00 ) ;
        self->vTotB /= n ; v2 = self->vTotB2 / n + self->vTot * ( self->vTot - 2.0e+00 * self->vTotB ) ; self->vTotB2 = Maximum ( v2, 0.0e+00 ) ;
    }
}
