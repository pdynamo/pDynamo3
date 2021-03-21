#ifndef _ENERGYMODEL
#define _ENERGYMODEL

/* Data types */
#include "Boolean.h"
#include "Cardinal.h"
#include "Integer.h"
#include "Real.h"

/* Arrays */
#include "IntegerArray1D.h"
#include "RealArray1D.h"
#include "RealArray2D.h"
#include "SymmetricMatrix.h"

/* Other */
#include "Constants.h"
#include "Status.h"
#include "Units.h"

/* Own modules */
#include "StateVector.h"

/* . Units. */
# define Constant_Molar_Gas_Kcalories_Per_Mole ( Constant_Molar_Gas / ( Units_Energy_Calories_To_Joules * 1000.0e+00 ) )

/* . Old value for testing - doesn't seem to make much difference. */
/*# define Constant_Molar_Gas_Kcalories_Per_Mole 0.001987165392*/

/* Macros */
#define EnergyModel_RowPointer(self, i) (&self->symmetricmatrix->data[(i * (i + 1) >> 1)])

#define EnergyModel_GetW(self, i, j) (i >= j ? self->symmetricmatrix->data[(i * (i + 1) >> 1) + j] : self->symmetricmatrix->data[(j * (j + 1) >> 1) + i])

typedef struct {
    /* Number of bound protons of each instance */
    IntegerArray1D   *protons;
    /* Gmodel of each instance (needed for energies of unfolded proteins) */
    RealArray1D      *models;
    /* Gintr of each instance */
    RealArray1D      *intrinsic;
    /* Interactions between instances before symmetrization */
    RealArray2D      *interactions;
    /* Symmetrized interactions */
    SymmetricMatrix  *symmetricmatrix;
    /* Probability of occurrence of each instance */
    RealArray1D      *probabilities;
    /* Private state vector of the energy model */
    StateVector      *vector;
    /* Total number of possible protonation states, no greater than ANALYTIC_STATES */
    Integer           nstates;
    /* Total number of instances */
    Integer           ninstances;
    /* Temperature at which the MEAD part was done */
    Real              temperature;
} EnergyModel;


/* Allocation and deallocation */
extern EnergyModel *EnergyModel_Allocate   (const Integer nsites, const Integer ninstances, Status *status);
extern void         EnergyModel_Deallocate (EnergyModel **self);

/* Miscellaneous functions */
extern void    EnergyModel_SymmetrizeInteractions       (const EnergyModel *self, Status *status);
extern Boolean EnergyModel_CheckInteractionsSymmetric   (const EnergyModel *self, Real tolerance, Real *maxDeviation);
extern void    EnergyModel_ResetInteractions            (const EnergyModel *self);
extern void    EnergyModel_ScaleInteractions            (const EnergyModel *self, Real scale);
extern void    EnergyModel_StateVectorFromProbabilities (const EnergyModel *self, StateVector *vector, Status *status);

/* Calculation of energies */
extern Real EnergyModel_CalculateMicrostateEnergy         (const EnergyModel *self, const StateVector *vector, const Real pH);
extern Real EnergyModel_CalculateMicrostateEnergyUnfolded (const EnergyModel *self, const StateVector *vector, const Real pH);

/* Calculation of partition functions */
extern Real EnergyModel_CalculateZ         (const EnergyModel *self, Real (*EnergyFunction)(const EnergyModel*, const StateVector*, const Real), const Real pH, const Real Gzero, RealArray1D *bfactors);
extern Real EnergyModel_CalculateZunfolded (const EnergyModel *self, const Real pH, const Real Gzero, Status *status);
extern Real EnergyModel_CalculateZfolded   (const EnergyModel *self, const Real pH, const Real Gzero, Status *status);

/* Calculation of probabilities */
extern void EnergyModel_CalculateProbabilitiesFromZ                (const EnergyModel *self, const Real Z, const RealArray1D *bfactors);
extern void EnergyModel_CalculateProbabilitiesAnalytically         (const EnergyModel *self, const Real pH, Status *status);
extern void EnergyModel_CalculateProbabilitiesAnalyticallyUnfolded (const EnergyModel *self, const Real pH, Status *status);

/* Functions for accessing items */
extern Real    EnergyModel_GetGmodel         (const EnergyModel *self, const Integer instIndexGlobal);
extern Real    EnergyModel_GetGintr          (const EnergyModel *self, const Integer instIndexGlobal);
extern Integer EnergyModel_GetProtons        (const EnergyModel *self, const Integer instIndexGlobal);
extern Real    EnergyModel_GetProbability    (const EnergyModel *self, const Integer instIndexGlobal);
extern Real    EnergyModel_GetInteraction    (const EnergyModel *self, const Integer instIndexGlobalA, const Integer instIndexGlobalB);
extern Real    EnergyModel_GetInterSymmetric (const EnergyModel *self, const Integer instIndexGlobalA, const Integer instIndexGlobalB);
extern Real    EnergyModel_GetDeviation      (const EnergyModel *self, const Integer instIndexGlobalA, const Integer instIndexGlobalB);

extern void    EnergyModel_SetGmodel         (const EnergyModel *self, const Integer instIndexGlobal, const Real value);
extern void    EnergyModel_SetGintr          (const EnergyModel *self, const Integer instIndexGlobal, const Real value);
extern void    EnergyModel_SetProtons        (const EnergyModel *self, const Integer instIndexGlobal, const Integer value);
extern void    EnergyModel_SetProbability    (const EnergyModel *self, const Integer instIndexGlobal, const Real value);
extern void    EnergyModel_SetInteraction    (const EnergyModel *self, const Integer instIndexGlobalA, const Integer instIndexGlobalB, const Real value);

#endif
