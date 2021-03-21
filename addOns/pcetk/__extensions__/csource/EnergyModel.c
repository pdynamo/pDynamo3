#include "Array_Macros.h"
#include "EnergyModel.h"
#include "Memory.h"
#include "RealIterator.h"

/*
 * Allocate the energy model.
 * Other attributes of the model (nstates, ninstances, temperature) are set from the Cython level.
 */
EnergyModel *EnergyModel_Allocate ( const Integer nsites, const Integer ninstances, Status *status )
{
    EnergyModel *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( EnergyModel ) ;
        if ( self != NULL )
        {
            Boolean isOK = True ;
            self->vector          = NULL ;
            self->models          = NULL ;
            self->protons         = NULL ;
            self->intrinsic       = NULL ;
            self->interactions    = NULL ;
            self->probabilities   = NULL ;
            self->symmetricmatrix = NULL ;
            if ( nsites > 0 )
            {
                self->vector = StateVector_Allocate ( nsites, status ) ;
                isOK = ( self->vector != NULL ) ;
            }
            if ( ninstances > 0 )
            {
                self->protons         = IntegerArray1D_AllocateWithExtent  ( ninstances,             status ) ;
                self->intrinsic       = RealArray1D_AllocateWithExtent     ( ninstances,             status ) ;
                self->models          = RealArray1D_AllocateWithExtent     ( ninstances,             status ) ;
                self->probabilities   = RealArray1D_AllocateWithExtent     ( ninstances,             status ) ;
                self->interactions    = RealArray2D_AllocateWithExtents    ( ninstances, ninstances, status ) ;
                self->symmetricmatrix = SymmetricMatrix_AllocateWithExtent ( ninstances,             status ) ;
                isOK = isOK && ( ( self->protons         != NULL ) && 
                                 ( self->intrinsic       != NULL ) && 
                                 ( self->models          != NULL ) && 
                                 ( self->probabilities   != NULL ) && 
                                 ( self->interactions    != NULL ) && 
                                 ( self->symmetricmatrix != NULL ) ) ;
            }
            if ( ! isOK ) EnergyModel_Deallocate ( &self ) ;
        }
        if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    }
    return self;
}

/*
 * Deallocate the energy model.
 */
void EnergyModel_Deallocate (EnergyModel **self)
{
    if ( ( self != NULL ) && ( (*self) != NULL) )
    {
        SymmetricMatrix_Deallocate ( &((*self)->symmetricmatrix) ) ;
        RealArray1D_Deallocate     ( &((*self)->probabilities  ) ) ;
        RealArray2D_Deallocate     ( &((*self)->interactions   ) ) ;
        RealArray1D_Deallocate     ( &((*self)->intrinsic      ) ) ;
        IntegerArray1D_Deallocate  ( &((*self)->protons        ) ) ;
        RealArray1D_Deallocate     ( &((*self)->models         ) ) ;
        StateVector_Deallocate     ( &((*self)->vector         ) ) ;
        Memory_Deallocate ((*self));
    }
}

/*
 * Check if the array of interactions is symmetric within the given tolerance (kcal/mol).
 */
Boolean EnergyModel_CheckInteractionsSymmetric (const EnergyModel *self, Real tolerance, Real *maxDeviation) {
    return RealArray2D_IsSymmetric (self->interactions, &tolerance, maxDeviation);
}

/*
 * Symmetrize the array of interactions into a symmetric matrix.
 */
void EnergyModel_SymmetrizeInteractions (const EnergyModel *self, Status *status) {
    SymmetricMatrix_CopyFromRealArray2D (self->symmetricmatrix, self->interactions, status);
}

/*
 * Set all interactions to zero.
 */
void EnergyModel_ResetInteractions (const EnergyModel *self) {
    SymmetricMatrix_Set (self->symmetricmatrix, 0.0f);
}

/*
 * Scale interactions.
 */
void EnergyModel_ScaleInteractions (const EnergyModel *self, Real scale) {
    SymmetricMatrix_Scale (self->symmetricmatrix, scale);
}

/*
 * Generate the lowest energy state vector.
 * If "vector" is NULL, use the EnergyModel's private vector.
 */
void EnergyModel_StateVectorFromProbabilities ( const EnergyModel *self, StateVector *vector, Status *status )
{
    if ( ( self   != NULL ) &&
         Status_IsOK ( status ) )
    {
        Boolean   isOK = False ;
        TitrSite *ts;
        Integer   i, index, maxi ;
        Real      probability, maxp ;
        if ( vector == NULL )
        {
            if ( self->vector != NULL )
            {
                ts = self->vector->sites  ;
                i  = self->vector->nsites ;
                isOK = True ;
            }
        }
        else
        {
            ts   = vector->sites  ;
            i    = vector->nsites ;
            isOK = ( i == self->vector->nsites ) ;
        }
        if ( isOK )
        {
            for (; i >= 0; i--, ts++) {
                index = ts->indexFirst;
                maxi  = index;
                maxp  = -1.0f;
                do {
                    probability = Array1D_Item ( self->probabilities, index ) ;
                    if ( probability > maxp ) {
                        maxi = index;
                        maxp = probability;
                    }
                } while (++index <= ts->indexLast);
                ts->indexActive = maxi;
            }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*
 * Getters.
 */
Real EnergyModel_GetGmodel (const EnergyModel *self, const Integer instIndexGlobal) {
    return Array1D_Item (self->models, instIndexGlobal);
}

Real EnergyModel_GetGintr (const EnergyModel *self, const Integer instIndexGlobal) {
    return Array1D_Item (self->intrinsic, instIndexGlobal);
}

Integer EnergyModel_GetProtons (const EnergyModel *self, const Integer instIndexGlobal) {
    return Array1D_Item (self->protons, instIndexGlobal);
}

Real EnergyModel_GetProbability (const EnergyModel *self, const Integer instIndexGlobal) {
    return Array1D_Item (self->probabilities, instIndexGlobal);
}

Real EnergyModel_GetInteraction (const EnergyModel *self, const Integer instIndexGlobalA, const Integer instIndexGlobalB) {
    return Array2D_Item (self->interactions, instIndexGlobalA, instIndexGlobalB);
}

Real EnergyModel_GetInterSymmetric (const EnergyModel *self, const Integer instIndexGlobalA, const Integer instIndexGlobalB) {
    return EnergyModel_GetW (self, instIndexGlobalA, instIndexGlobalB);
}

Real EnergyModel_GetDeviation (const EnergyModel *self, const Integer i, const Integer j) {
    Real wij, wji, deviation;
    wij = Array2D_Item (self->interactions, i, j);
    wji = Array2D_Item (self->interactions, j, i);
    deviation = (wij + wji) * .5 - wij;
    return deviation;
}

/*
 * Setters.
 */
void EnergyModel_SetGmodel (const EnergyModel *self, const Integer instIndexGlobal, const Real value) {
    Array1D_Item (self->models, instIndexGlobal) = value;
}

void EnergyModel_SetGintr (const EnergyModel *self, const Integer instIndexGlobal, const Real value) {
    Array1D_Item (self->intrinsic, instIndexGlobal) = value;
}

void EnergyModel_SetProtons (const EnergyModel *self, const Integer instIndexGlobal, const Integer value) {
    Array1D_Item (self->protons, instIndexGlobal) = value;
}

void EnergyModel_SetProbability (const EnergyModel *self, const Integer instIndexGlobal, const Real value) {
    Array1D_Item (self->probabilities, instIndexGlobal) = value;
}

void EnergyModel_SetInteraction (const EnergyModel *self, const Integer instIndexGlobalA, const Integer instIndexGlobalB, const Real value) {
    Array2D_Item (self->interactions, instIndexGlobalA, instIndexGlobalB) = value;
}

/*
 * Calculate the energy of a microstate defined by the state vector.
 */
Real EnergyModel_CalculateMicrostateEnergy ( const EnergyModel *self, const StateVector *vector, const Real pH )
{
    Real Gintr = 0.0f ;
    if ( ( self != NULL ) && ( vector != NULL ) )
    {
        Real      W, *interact;
        Integer   nprotons, i, j;
        TitrSite *site, *siteInner;
        W        = 0.0f;
        nprotons = 0;
        site     = vector->sites;
        for (i = 0; i < vector->nsites; i++, site++)
        {
            Gintr     += Array1D_Item (self->intrinsic , site->indexActive);
            nprotons  += Array1D_Item (self->protons   , site->indexActive);
            interact   = EnergyModel_RowPointer (self, site->indexActive);
            siteInner  = vector->sites;
            for (j = 0; j < i; j++, siteInner++) { W += *(interact + (siteInner->indexActive)); }
        }
        Gintr += ( ( ( Real ) nprotons ) * ( Constant_Molar_Gas_Kcalories_Per_Mole * self->temperature * Constant_Ln10 * pH ) + W ) ;
    }
    return Gintr ;
}

/*
 * Calculate the energy of a microstate in an unfolded (=denaturated) protein.
 * In the unfolded state, Gintr become Gmodel and all interactions are set to zero.
 * 
 * Reference: Yang A.-S., Honig B., J. Mol. Biol. 1993, 231, 459-474
 */
Real EnergyModel_CalculateMicrostateEnergyUnfolded ( const EnergyModel *self, const StateVector *vector, const Real pH )
{
    Real Gmodel = 0.0f ;
    if ( ( self != NULL ) && ( vector != NULL ) )
    {
        Integer   nprotons, i;
        TitrSite *site;
        nprotons = 0;
        site     = vector->sites;
        for ( i = 0; i < vector->nsites; i++, site++ )
        {
            Gmodel   += Array1D_Item ( self->models  , site->indexActive ) ;
            nprotons += Array1D_Item ( self->protons , site->indexActive ) ;
        }
        Gmodel -= ( ( Real ) nprotons ) * ( Constant_Molar_Gas_Kcalories_Per_Mole * self->temperature * Constant_Ln10 * pH ) ;
    }
    return Gmodel ;
}

/*
 * Calculate partition function and Boltzmann factors using a custom energy function.
 *
 * Note: bfactors should be allocated beforehand.
 */
Real EnergyModel_CalculateZ ( const EnergyModel *self     , 
                                    Real       (*EnergyFunction)(const EnergyModel*, const StateVector*, const Real), 
                              const Real         pH       ,
                              const Real         Gzero    , 
                                    RealArray1D *bfactors )
{
    Real Z = 0.0f;
    if ( ( self != NULL ) && ( bfactors != NULL ) )
    {
        Real     *bfactor, G, Gmin ;
        Integer   i ;
        StateVector_Reset ( self->vector ) ;
        i       = self->nstates;
        bfactor = Array_DataPointer (bfactors);
        Gmin    = EnergyFunction ( self, self->vector, pH ) - Gzero ;
        for (; i > 0; i--, bfactor++) {
            /* TODO: Optimize.
             * Do not calculate Gmicro after every increment of the state vector, 
             * calculate deltas like in MC moves.
             */
            G = EnergyFunction ( self, self->vector, pH ) - Gzero ;
            if ( G < Gmin ) Gmin = G ;
            *bfactor = G;
            StateVector_Increment ( self->vector ) ;
        }
        RealArray1D_Increment ( bfactors, -Gmin );
        RealArray1D_Scale     ( bfactors, -1.0f / ( Constant_Molar_Gas_Kcalories_Per_Mole * self->temperature ) ) ;
        /* . Exponential. */
        {
            auto Iterator *bIterator ;
            bIterator = View1D_MakeIterator ( ( View1D * ) bfactors, NULL ) ;
            RealIterator_Exponential ( bIterator, Array_DataPointer ( bfactors ), NULL ) ;
            Iterator_Deallocate ( &bIterator ) ;
        }
        Z = RealArray1D_Sum (bfactors);
    }
    return Z;
}

/*
 * Calculate the statistical mechanical partition function of an unfolded protein.
 */
Real EnergyModel_CalculateZunfolded ( const EnergyModel *self, const Real pH, const Real Gzero, Status *status )
{
    Real Z = -1.0f ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        RealArray1D *bfactors = RealArray1D_AllocateWithExtent ( self->nstates, status ) ;
        if ( bfactors != NULL )
        {
            Z = EnergyModel_CalculateZ ( self, EnergyModel_CalculateMicrostateEnergyUnfolded, pH, Gzero, bfactors ) ;
            RealArray1D_Deallocate ( &bfactors ) ;
        }
    }
    return Z ;
}

/*
 * Calculate the statistical mechanical partition function of a folded protein.
 */
Real EnergyModel_CalculateZfolded ( const EnergyModel *self, const Real pH, const Real Gzero, Status *status )
{
    Real Z = -1.0f ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        RealArray1D *bfactors = RealArray1D_AllocateWithExtent ( self->nstates, status ) ;
        if ( bfactors != NULL )
        {
            Z = EnergyModel_CalculateZ ( self, EnergyModel_CalculateMicrostateEnergy, pH, Gzero, bfactors ) ;
            RealArray1D_Deallocate ( &bfactors ) ;
        }
    }
    return Z ;
}

/*
 * Calculate protonation state probabilities from the statistical mechanical partition function.
 */
void EnergyModel_CalculateProbabilitiesFromZ ( const EnergyModel *self, const Real Z, const RealArray1D *bfactors )
{
    if ( ( self != NULL ) && ( bfactors != NULL ) )
    {
        Real      *bfactor;
        Integer    i, j;
        TitrSite  *ts;
        RealArray1D_Set ( self->probabilities, 0.0f );
        StateVector_Reset ( self->vector );
        i = self->nstates;
        bfactor = Array_DataPointer ( bfactors );
        for (; i > 0; i--, bfactor++) {
            j  = self->vector->nsites ;
            ts = self->vector->sites  ;
            for (; j > 0; j--, ts++) {
                Array1D_Item ( self->probabilities, ts->indexActive) += *bfactor ;
            }
            StateVector_Increment (self->vector);
        }
        RealArray1D_Scale ( self->probabilities, 1.0f / Z );
    }
}

/*
 * Analytic evaluation of protonation state probabilities.
 */
void EnergyModel_CalculateProbabilitiesAnalytically ( const EnergyModel *self, const Real pH, Status *status )
{
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        RealArray1D *bfactors = RealArray1D_AllocateWithExtent ( self->nstates, status ) ;
        if ( bfactors != NULL )
        {
            Real Z = EnergyModel_CalculateZ ( self, EnergyModel_CalculateMicrostateEnergy, pH, 0.0f, bfactors ) ;
            EnergyModel_CalculateProbabilitiesFromZ ( self, Z, bfactors ) ;
            RealArray1D_Deallocate ( &bfactors ) ;
        }
    }
}

/*
 * Analytic evaluation of protonation state probabilities (unfolded protein).
 */
void EnergyModel_CalculateProbabilitiesAnalyticallyUnfolded ( const EnergyModel *self, const Real pH, Status *status )
{
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        RealArray1D *bfactors = RealArray1D_AllocateWithExtent ( self->nstates, status ) ;
        if ( bfactors != NULL )
        {
            Real Z = EnergyModel_CalculateZ ( self, EnergyModel_CalculateMicrostateEnergyUnfolded, pH, 0.0f, bfactors ) ;
            EnergyModel_CalculateProbabilitiesFromZ ( self, Z, bfactors ) ;
            RealArray1D_Deallocate ( &bfactors ) ;
        }
    }
}
