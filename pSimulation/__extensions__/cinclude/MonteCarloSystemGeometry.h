# ifndef _MONTECARLOSYSTEMGEOMETRY
# define _MONTECARLOSYSTEMGEOMETRY

# include "Coordinates3.h"
# include "Integer.h"
# include "IntegerArray1D.h"
# include "Matrix33.h"
# include "PairwiseInteractionMonteCarlo.h"
# include "Real.h"
# include "RealArray1D.h"
# include "SelectionContainer.h"
# include "Status.h"
# include "SymmetryParameters.h"
# include "Vector3.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The Monte Carlo system geometry type. */
typedef struct {
    /* . Counters. */
    Integer                        blocks                ;
    Integer                        moves                 ;
    Integer                        nReject               ;
    Integer                        nRejectM              ;
    Integer                        nRejectT              ;
    Integer                        nRejectV              ;
    Integer                        nTryM                 ;
    Integer                        nTryV                 ;
    /* . Current values and other factors. */
    Real                           beta                  ;
    Real                           dielectric            ;
    Real                           eCurrent              ;
    Real                           pressure              ;
    Real                           tFactor               ;
    Real                           volume                ;
    /* . Move data. */
    Real                           acceptanceRatio       ;
    Real                           rMax                  ;
    Real                           tMax                  ;
    Real                           vMax                  ;
    /* . Statistics. */
    Real                           eAv                   ;
    Real                           eAv2                  ;
    Real                           eTot                  ;
    Real                           eTot2                 ;
    Real                           eTotB                 ;
    Real                           eTotB2                ;
    Real                           hAv                   ;
    Real                           hAv2                  ;
    Real                           hTot                  ;
    Real                           hTot2                 ;
    Real                           hTotB                 ;
    Real                           hTotB2                ;
    Real                           vAv                   ;
    Real                           vAv2                  ;
    Real                           vTot                  ;
    Real                           vTot2                 ;
    Real                           vTotB                 ;
    Real                           vTotB2                ;
    /* . Arrays to allocate. */
    Real                          *random                ;
    Coordinates3                  *oldCoordinates3       ;
    Matrix33                      *rotation              ;
    SymmetryParameters            *oldSymmetryParameters ;
    Vector3                       *translation           ;
    /* . Aliases. */
    IntegerArray1D                *ljTypes               ;
    Coordinates3                  *coordinates3          ;
    LJParameterContainer          *ljParameters          ;
    PairwiseInteractionMonteCarlo *pairwiseInteraction   ;
    RealArray1D                   *charges               ;
    SelectionContainer            *isolates              ;
    SymmetryParameters            *symmetryParameters    ;
} MonteCarloSystemGeometry ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern MonteCarloSystemGeometry *MonteCarloSystemGeometry_Allocate                  ( const Integer  numberOfParticles ,
                                                                                      const Integer  numberOfRandom    ) ;
extern void                      MonteCarloSystemGeometry_Deallocate                ( MonteCarloSystemGeometry **self ) ; 
extern void                      MonteCarloSystemGeometry_AdjustMoveSizes           ( MonteCarloSystemGeometry  *self ) ; 
extern Status                    MonteCarloSystemGeometry_MoveIsolate               ( MonteCarloSystemGeometry  *self ) ; 
extern Status                    MonteCarloSystemGeometry_MoveVolume                ( MonteCarloSystemGeometry  *self ) ; 
extern void                      MonteCarloSystemGeometry_StatisticsBlockAccumulate ( MonteCarloSystemGeometry  *self ) ; 
extern void                      MonteCarloSystemGeometry_StatisticsBlockStart      ( MonteCarloSystemGeometry  *self ) ; 
extern void                      MonteCarloSystemGeometry_StatisticsBlockStop       ( MonteCarloSystemGeometry  *self ) ; 
extern void                      MonteCarloSystemGeometry_StatisticsStop            ( MonteCarloSystemGeometry  *self ) ; 
extern void                      MonteCarloSystemGeometry_StatisticsStart           ( MonteCarloSystemGeometry  *self ) ; 

# endif
