#ifndef _STATEVECTOR
#define _STATEVECTOR

/* Needed for calloc */
#include <stdlib.h>
/* Needed for memcpy */
#include <string.h>

#include "Boolean.h"
#include "Integer.h"
#include "RandomNumberGenerator.h"
#include "Real.h"
#include "Status.h"

typedef struct {
    /* Site belongs to a substate */
    Boolean isSubstate;
    /* Index of the site itself */
    Integer indexSite;
    /* Global index of the currently active instance of the site */
    Integer indexActive;
    /* Minimum and maximum values of indexActive */
    Integer indexFirst, indexLast;
} TitrSite;

typedef struct {
    /* Pointers to sites that make up a pair */
    TitrSite *a, *b;
    /* Maximum absolute energy of interaction */
    Real Wmax;
} PairSite;

typedef struct {
    TitrSite  *sites, **substateSites;
    Integer    nsites, nssites;
    /* Handled by MCModelDefault */
    PairSite  *pairs;
    Integer    npairs;
} StateVector;


/* Allocation and deallocation */
extern StateVector *StateVector_Allocate          (const Integer nsites, Status *status);
extern void         StateVector_AllocateSubstate  (      StateVector  *self, const Integer nssites, Status *status);
extern void         StateVector_AllocatePairs     (      StateVector  *self, const Integer npairs, Status *status);
extern void         StateVector_ReallocatePairs   (      StateVector  *self, const Integer npairs, Status *status);
extern void         StateVector_Deallocate        (      StateVector **self);

/* Copying and cloning */
extern StateVector *StateVector_Clone             (const StateVector *self, Status *status);
extern void         StateVector_CopyTo            (const StateVector *self, StateVector *other, Status *status);

/* Functions for setting all items at once */
extern void         StateVector_Reset             (const StateVector *self);
extern void         StateVector_ResetSubstate     (const StateVector *self);
extern void         StateVector_ResetToMaximum    (const StateVector *self);
extern void         StateVector_Randomize         (const StateVector *self, const RandomNumberGenerator *generator);

/* Functions for accessing items */
extern void         StateVector_SetSite           (const StateVector *self, const Integer indexSite, const Integer indexFirst, const Integer indexLast, Status *status);
extern void         StateVector_SetPair           (const StateVector *self, const Integer indexPair, const Integer indexFirstSite, const Integer indexSecondSite, const Real Wmax, Status *status);
extern void         StateVector_GetPair           (const StateVector *self, const Integer indexPair, Integer *indexFirstSite, Integer *indexSecondSite, Real *Wmax, Status *status);
extern Boolean      StateVector_IsSubstate        (const StateVector *self, const Integer siteIndex, Status *status);
extern Integer      StateVector_GetItem           (const StateVector *self, const Integer siteIndex, Status *status);
extern void         StateVector_SetItem           (const StateVector *self, const Integer siteIndex, const Integer instanceLocalIndex, Status *status);
extern Integer      StateVector_GetActualItem     (const StateVector *self, const Integer siteIndex, Status *status);
extern void         StateVector_SetActualItem     (const StateVector *self, const Integer siteIndex, const Integer instanceGlobalIndex, Status *status);
extern Integer      StateVector_GetSubstateItem   (const StateVector *self, const Integer index, Status *status);
extern void         StateVector_SetSubstateItem   (const StateVector *self, const Integer selectedSiteIndex, const Integer index, Status *status);

/* Incrementation */
extern Boolean      StateVector_Increment         (const StateVector *self);
extern Boolean      StateVector_IncrementSubstate (const StateVector *self);

#endif
