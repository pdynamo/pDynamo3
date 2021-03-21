#ifndef _MCMODELDEFAULT
#define _MCMODELDEFAULT

/* Needed for random seed */
#include <time.h>
/* Needed for exp */
#include <math.h>

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
#include "RandomNumberGenerator.h"
#include "Status.h"

/* Own modules */
#include "EnergyModel.h"
#include "StateVector.h"

/* Taken from GMCT */
#define TOO_SMALL -500.0

typedef struct {
    /* Energy limit for double moves */
    Real                   limit;
    /* Number of equlibration scans */
    Integer                nequil;
    /* Number of production scans */
    Integer                nprod;
    /* Pointer to the energy model */
    EnergyModel           *energyModel;
    /* Private state vector of the Monte Carlo model */
    StateVector           *vector;
    /* Mersenne Twister generator */
    RandomNumberGenerator *generator;
} MCModelDefault;


/* Allocation and deallocation */
extern MCModelDefault *MCModelDefault_Allocate          (const Real limit, const Integer nequil, const Integer nprod, const Integer randomSeed, Status *status);
extern void            MCModelDefault_Deallocate        (MCModelDefault **self);
extern void            MCModelDefault_LinkToEnergyModel (MCModelDefault *self, EnergyModel *energyModel, Status *status);

/* Monte Carlo-related functions */
extern Boolean MCModelDefault_Metropolis             (const Real GdeltaRT, const RandomNumberGenerator *generator);
extern Boolean MCModelDefault_Move                   (const MCModelDefault *self, const Real pH, const Real G, Real *Gnew);
extern Boolean MCModelDefault_DoubleMove             (const MCModelDefault *self, const Real pH, const Real G, Real *Gnew);
extern Real    MCModelDefault_MCScan                 (const MCModelDefault *self, const Real pH, Integer nmoves, Integer *movesDone, Integer *movesAccepted, Integer *flipsDone, Integer *flipsAccepted);
extern void    MCModelDefault_UpdateProbabilities    (const MCModelDefault *self);
extern Real    MCModelDefault_FindMaxInteraction     (const MCModelDefault *self, const TitrSite *site, const TitrSite *other);
extern Integer MCModelDefault_FindPairs              (const MCModelDefault *self, const Integer npairs, Status *status);
extern void    MCModelDefault_Equilibration          (const MCModelDefault *self, const Real pH);
extern void    MCModelDefault_Production             (const MCModelDefault *self, const Real pH);

#endif
