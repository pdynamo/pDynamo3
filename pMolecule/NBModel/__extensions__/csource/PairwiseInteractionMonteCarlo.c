/*==================================================================================================================================
! . A simple pairwise interaction model for Monte Carlo simulations.
! . There are no derivatives, no intraisolate interactions and no QC!
!=================================================================================================================================*/

# include <math.h>
# include "Memory.h"
# include "NumericalMacros.h"
# include "PairwiseInteractionMonteCarlo.h"
# include "Units.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _DefaultBuffer       0.5e+00
# define _DefaultChargeScale  1.0e+00
# define _DefaultCutOff       8.5e+00
# define _DefaultEpsilonScale 1.0e+00
# define _DefaultSigmaScale   1.0e+00
# define _DefaultUnderFlowL   0.2e+00
# define _DefaultUnderFlowQ   0.5e+00

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairwiseInteractionMonteCarlo *PairwiseInteractionMonteCarlo_Allocate ( Status *status )
{
    PairwiseInteractionMonteCarlo *self = Memory_AllocateType ( PairwiseInteractionMonteCarlo ) ;
    if ( self != NULL )
    {
        self->isolateScale = 0 ;
        self->buffer       = _DefaultBuffer       ;
        self->chargeScale  = _DefaultChargeScale  ;
        self->cutOff       = _DefaultCutOff       ;
        self->epsilonScale = _DefaultEpsilonScale ;
        self->sigmaScale   = _DefaultSigmaScale   ;
        self->underFlowL   = _DefaultUnderFlowL   ;
        self->underFlowQ   = _DefaultUnderFlowQ   ;
    }
    else Status_Set ( status, Status_OutOfMemory ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
PairwiseInteractionMonteCarlo *PairwiseInteractionMonteCarlo_Clone ( const PairwiseInteractionMonteCarlo *self   ,
                                                                           Status                        *status )
{
    PairwiseInteractionMonteCarlo *new = NULL ;
    if ( self != NULL )
    {
        new = PairwiseInteractionMonteCarlo_Allocate ( status ) ;
        new->buffer       = self->buffer       ;
        new->chargeScale  = self->chargeScale  ;
        new->cutOff       = self->cutOff       ;
        new->epsilonScale = self->epsilonScale ;
        new->sigmaScale   = self->sigmaScale   ;
        new->underFlowL   = self->underFlowL   ;
        new->underFlowQ   = self->underFlowQ   ;
        new->buffer     = self->buffer     ;
        new->cutOff     = self->cutOff     ;
    }
    else Status_Set ( status, Status_OutOfMemory ) ;
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionMonteCarlo_Deallocate ( PairwiseInteractionMonteCarlo **self )
{
    if ( (*self) != NULL ) Memory_Deallocate ( (*self) ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Interactions.
!---------------------------------------------------------------------------------------------------------------------------------*/
void PairwiseInteractionMonteCarlo_Interactions ( const PairwiseInteractionMonteCarlo  *self          ,
                                                  const RealArray1D                    *r             ,
                                                        RealArray1D                    *electrostatic ,
                                                        RealArray1D                    *lennardJonesA ,
                                                        RealArray1D                    *lennardJonesB )
{
    Boolean doElectrostatic = ( electrostatic != NULL ) ,
            doLennardJonesA = ( lennardJonesA != NULL ) ,
            doLennardJonesB = ( lennardJonesB != NULL ) ;
    if ( ( self != NULL ) && ( r != NULL ) && ( doElectrostatic || doLennardJonesA || doLennardJonesB ) )
    {
        auto Integer     i ;
        auto Real        fA, fB, fQ, s2, s6, underFlowL2, underFlowQ2, x, x2, xL2, xQ2 ;
        underFlowL2 = pow ( self->underFlowL, 2 ) ;
        underFlowQ2 = pow ( self->underFlowQ, 2 ) ;
        RealArray1D_Set ( electrostatic, 0.0e+00 ) ;
        RealArray1D_Set ( lennardJonesA, 0.0e+00 ) ;
        RealArray1D_Set ( lennardJonesB, 0.0e+00 ) ;
        for ( i = 0 ; i < View1D_Extent ( r ) ; i++ )
        {
            x   = Array1D_Item ( r, i ) ;
            x2  = x * x ;
            xL2 = Maximum ( x2, underFlowL2 ) ;
            xQ2 = Maximum ( x2, underFlowQ2 ) ;
	    s2  = 1.0e+00 / xL2 ;
            s6  = s2 * s2 * s2 ;
            fA  = s6 * s6 ;
            fB  = - s6 ;
            fQ  = 1.0e+00 / sqrt ( xQ2 ) ;
            if ( doElectrostatic ) Array1D_Item ( electrostatic, i ) = fQ ;
            if ( doLennardJonesA ) Array1D_Item ( lennardJonesA, i ) = fA ;
            if ( doLennardJonesB ) Array1D_Item ( lennardJonesB, i ) = fB ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . MM/MM energy.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real PairwiseInteractionMonteCarlo_MMMMEnergy ( const PairwiseInteractionMonteCarlo  *self               ,
                                                      RealArray1D                    *charges            ,
                                                      IntegerArray1D                 *ljTypes            ,
                                                      LJParameterContainer           *ljParameters       ,
                                                const Real                            electrostaticScale ,
                                                const Real                            lennardJonesScale  ,
                                                      SelectionContainer             *isolates           ,
                                                      BooleanBlock                   *isFree             ,
                                                      Coordinates3                   *coordinates3       ,
                                                      SymmetryParameters             *symmetryParameters ,
                                                      Real                           *eElectrostatic     ,
                                                      Real                           *eLennardJones      ,
                                                      Status                         *status             )
{
    Real eL = 0.0e+00, eQ = 0.0e+00 ;
/*printf ( "MCNB> %p %p %p %p %p %p %p %d\n", self, charges, coordinates3, isolates, ljParameters, ljTypes, symmetryParameters, (*status) ) ;*/
    if ( ( self               != NULL ) &&
         ( charges            != NULL ) &&
         ( coordinates3       != NULL ) &&
         ( isolates           != NULL ) &&
         ( ljParameters       != NULL ) &&
         ( ljTypes            != NULL ) &&
         ( symmetryParameters != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Boolean       allFree = ( isFree == NULL ), isFreeI ;
        auto Integer       i, j, m, n, nLJTypes, s, t, tI, tIJ ;
        auto Real          cut2, eL0, eQ0, eScale, lower2, qI, qScale,
                           r2, rQ2, rIJ2, rL2, scale, sScale, s2, s6, underFlowL2, underFlowQ2, xIJ, yIJ, zIJ ;
        auto Real          dV[3], v[3] ;
        auto Coordinates3 *centers = NULL ;
        auto Selection    *iSelection, *jSelection ;
        /* . Initialization. */
        cut2        = self->cutOff * self->cutOff ;
        lower2      = pow ( self->cutOff - self->buffer, 2 ) ;
        nLJTypes    = ljParameters->ntypes        ;
        underFlowL2 = pow ( self->underFlowL, 2 ) ;
        underFlowQ2 = pow ( self->underFlowQ, 2 ) ;
        centers     = Coordinates3_Allocate ( isolates->capacity, status ) ;
        if ( centers == NULL ) goto FinishUp ;
        /* . Double loop over isolates - at least one isolate must be free. */
        for ( i = 0 ; i < isolates->capacity ; i++ )
        {
            iSelection = isolates->items[i] ;
            isFreeI    = ( allFree || Block_Item ( isFree, i ) ) ;
            Coordinates3_CenterRaw ( coordinates3, iSelection, NULL, Coordinates3_RowPointer ( centers, i ), NULL ) ;
            for ( j = 0 ; j < i ; j++ )
            {
                if ( isFreeI || Block_Item ( isFree, j ) )
                {
                    jSelection = isolates->items[j] ;
                    /* . Calculate the distance between the isolate centers (applying the minimum image convention). */
                    Coordinates3_DifferenceRow ( centers, i, j, v[0], v[1], v[2] ) ;
                    SymmetryParameters_MakeMinimumImageVector ( symmetryParameters, v, dV ) ;
                    rIJ2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2] ;
                    /* . Check to see whether the molecules are within the cutoff distance. */
                    if ( rIJ2 < cut2 )
                    {
                        /* . Check for scaling. */
                        if ( ( self->isolateScale == i ) || ( self->isolateScale == j ) )
                        {
                            qScale = self->chargeScale  ;
                            eScale = self->epsilonScale ;
                            sScale = self->sigmaScale   ;
                        }
                        else { qScale = eScale = sScale = 1.0e+00 ; }
                        /* . Double loop over atoms. */
                        for ( s = 0, eQ0 = eL0 = 0.0e+00 ; s < iSelection->capacity ; s++ )
                        {
                            m  = iSelection->indices[s] ;
  	                    qI =            Array1D_Item   ( charges, m ) ;
                            tI = nLJTypes * Array1D_Item ( ljTypes, m ) ;
                            for ( t = 0 ; t < jSelection->capacity ; t++ )
	                    {
                                n = jSelection->indices[t] ;
                                /* . Coordinate displacement. */
	                        Coordinates3_DifferenceRow ( coordinates3, m, n, xIJ, yIJ, zIJ ) ;
                                xIJ += dV[0] ;
                                yIJ += dV[1] ;
                                zIJ += dV[2] ;
                                r2   = xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ;
                                rL2  = Maximum ( r2, underFlowL2 ) ;
                                rQ2  = Maximum ( r2, underFlowQ2 ) ;
                                /* . Distance factors. */
	                        s2   = 1.0e+00 / rL2 ;
                                s6   = s2 * s2 * s2 * sScale ;
                                /* . Electrostatic. */
                                eQ0 += qI * Array1D_Item ( charges, n ) / sqrt ( rQ2 ) ;
                                /* . LJ. */
	                        tIJ  = ljParameters->tableindex[tI+Array1D_Item ( ljTypes, n )] ;
	                        eL0 += ( ljParameters->tableA[tIJ] * s6 - ljParameters->tableB[tIJ] ) * s6 ;
                            }
                        }
                        /* . Calculate the scale factor. */
                        if ( rIJ2 > lower2 ) scale = ( cut2 - rIJ2 ) / ( cut2 - lower2 ) ;
                        else                 scale = 1.0e+00                             ;
                        /* . Accumulate the total energies. */
                        eL += eL0 * eScale * scale ;
                        eQ += eQ0 * qScale * scale ;
                    }
                }
            }
        }
        eL *= lennardJonesScale ;
        eQ *= ( electrostaticScale * Units_Energy_E2Angstroms_To_Kilojoules_Per_Mole ) ;
        Coordinates3_Deallocate ( &centers ) ;
    }
FinishUp:
    if ( eElectrostatic != NULL ) (*eElectrostatic) = eQ ;
    if ( eLennardJones  != NULL ) (*eLennardJones ) = eL ;
/*printf ( "%f %f\n", eL, eQ ) ;*/
    return ( eL + eQ ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Single isolate electrostatic and LJ MM/MM energy and gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real PairwiseInteractionMonteCarlo_MMMMIsolateEnergy ( const PairwiseInteractionMonteCarlo  *self               ,
                                                       const Integer                         isolate            ,
                                                             RealArray1D                    *charges            ,
                                                             IntegerArray1D                 *ljTypes            ,
                                                             LJParameterContainer           *ljParameters       ,
                                                       const Real                            electrostaticScale ,
                                                       const Real                            lennardJonesScale  ,
                                                             SelectionContainer             *isolates           ,
                                                             BooleanBlock                   *isFree             ,
                                                             Coordinates3                   *coordinates3       ,
                                                             SymmetryParameters             *symmetryParameters ,
                                                             Real                           *eElectrostatic     ,
                                                             Real                           *eLennardJones      ,
                                                             Status                         *status             )
{
    Real eL = 0.0e+00, eQ = 0.0e+00 ;
    if ( ( self               != NULL ) &&
         ( charges            != NULL ) &&
         ( coordinates3       != NULL ) &&
         ( isolates           != NULL ) &&
         ( ljParameters       != NULL ) &&
         ( ljTypes            != NULL ) &&
         ( symmetryParameters != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Boolean       allFree = ( isFree == NULL ), isFreeI ;
        auto Integer       i, j, m, n, nLJTypes, s, t, tI, tIJ ;
        auto Real          cut2, eL0, eQ0, eScale, lower2, qI, qScale,
                           r2, rQ2, rIJ2, rL2, scale, sScale, s2, s6, underFlowL2, underFlowQ2, xIJ, yIJ, zIJ ;
        auto Real          dV[3], v[3] ;
        auto Coordinates3 *centers = NULL ;
        auto Selection    *iSelection, *jSelection ;
        /* . Initialization. */
        cut2        = self->cutOff * self->cutOff ;
        lower2      = pow ( self->cutOff - self->buffer, 2 ) ;
        nLJTypes    = ljParameters->ntypes        ;
        underFlowL2 = pow ( self->underFlowL, 2 ) ;
        underFlowQ2 = pow ( self->underFlowQ, 2 ) ;
        centers     = Coordinates3_Allocate ( 2, status ) ;
        if ( centers == NULL ) goto FinishUp ;
        /* . Loop over isolates (!= isolate) - at least one isolate must be free. */
        i          = isolate ;
        iSelection = isolates->items[i] ;
        isFreeI    = ( allFree || Block_Item ( isFree, i ) ) ;
        Coordinates3_CenterRaw ( coordinates3, iSelection, NULL, Coordinates3_RowPointer ( centers, 0 ), NULL ) ;
        for ( j = 0 ; j < isolates->capacity ; j++ )
        {
            if ( ( j != i ) && ( isFreeI || Block_Item ( isFree, j ) ) )
            {
                jSelection = isolates->items[j] ;
                Coordinates3_CenterRaw ( coordinates3, jSelection, NULL, Coordinates3_RowPointer ( centers, 1 ), NULL ) ;
                /* . Calculate the distance between the isolate centers (applying the minimum image convention). */
                Coordinates3_DifferenceRow ( centers, 0, 1, v[0], v[1], v[2] ) ;
                SymmetryParameters_MakeMinimumImageVector ( symmetryParameters, v, dV ) ;
                rIJ2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2] ;
                /* . Check to see whether the molecules are within the cutoff distance. */
                if ( rIJ2 < cut2 )
                {
                    /* . Check for scaling. */
                    if ( ( self->isolateScale == i ) || ( self->isolateScale == j ) )
                    {
                        qScale = self->chargeScale  ;
                        eScale = self->epsilonScale ;
                        sScale = self->sigmaScale   ;
                    }
                    else { qScale = eScale = sScale = 1.0e+00 ; }
                    /* . Double loop over atoms. */
                    for ( s = 0, eQ0 = eL0 = 0.0e+00 ; s < iSelection->capacity ; s++ )
                    {
                        m  = iSelection->indices[s] ;
  	                qI =            Array1D_Item   ( charges, m ) ;
                        tI = nLJTypes * Array1D_Item ( ljTypes, m ) ;
                        for ( t = 0 ; t < jSelection->capacity ; t++ )
	                {
                            n = jSelection->indices[t] ;
                            /* . Coordinate displacement. */
	                    Coordinates3_DifferenceRow ( coordinates3, m, n, xIJ, yIJ, zIJ ) ;
                            xIJ += dV[0] ;
                            yIJ += dV[1] ;
                            zIJ += dV[2] ;
                            r2   = xIJ * xIJ + yIJ * yIJ + zIJ * zIJ ;
                            rL2  = Maximum ( r2, underFlowL2 ) ;
                            rQ2  = Maximum ( r2, underFlowQ2 ) ;
                            /* . Distance factors. */
	                    s2   = 1.0e+00 / rL2 ;
                            s6   = s2 * s2 * s2 * sScale ;
                            /* . Electrostatic. */
                            eQ0 += qI * Array1D_Item ( charges, n ) / sqrt ( rQ2 ) ;
                            /* . LJ. */
	                    tIJ  = ljParameters->tableindex[tI+Array1D_Item ( ljTypes, n )] ;
	                    eL0 += ( ljParameters->tableA[tIJ] * s6 - ljParameters->tableB[tIJ] ) * s6 ;
                        }
                    }
                    /* . Calculate the scale factor. */
                    if ( rIJ2 > lower2 ) scale = ( cut2 - rIJ2 ) / ( cut2 - lower2 ) ;
                    else                 scale = 1.0e+00                             ;
                    /* . Accumulate the total energies. */
                    eL += eL0 * eScale * scale ;
                    eQ += eQ0 * qScale * scale ;
                }
            }
        }
        eL *= lennardJonesScale ;
        eQ *= ( electrostaticScale * Units_Energy_E2Angstroms_To_Kilojoules_Per_Mole ) ;
        Coordinates3_Deallocate ( &centers ) ;
    }
FinishUp:
    if ( eElectrostatic != NULL ) (*eElectrostatic) = eQ ;
    if ( eLennardJones  != NULL ) (*eLennardJones ) = eL ;
    return ( eL + eQ ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Scale the isolate interaction parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* It might be preferable to scale the force field parameters directly instead of doing individual scalings.
   This is straightforward for the charges, but more complicated for the LJ parameters as would need to expand
   the number of types.
*/
void PairwiseInteractionMonteCarlo_ScaleIsolateInteractionParameters (       PairwiseInteractionMonteCarlo  *self         ,
                                                                       const Integer                         isolate      ,
                                                                       const Real                            chargeScale  ,
                                                                       const Real                            epsilonScale ,
                                                                       const Real                            sigmaScale   )
{
    if ( self != NULL )
    {
        self->isolateScale = isolate ;
        self->chargeScale  = chargeScale                    ; /* . Scaling for qi * qj.           */
        self->epsilonScale = sqrt ( fabs ( epsilonScale ) ) ; /* . Scaling for sqrt ( ei * ej ).  */
        self->sigmaScale   = pow ( sigmaScale, 3.0 )        ; /* . Scaling for ( si**3 * sj**3 ). */
    }
}
