# ifndef _HARMONICBONDCONTAINER
# define _HARMONICBONDCONTAINER

# include "Boolean.h"
# include "Coordinates3.h"
# include "Integer.h"
# include "Real.h"
# include "Selection.h"

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
typedef struct {
    Boolean isActive ;
    Integer atom1    ;
    Integer atom2    ;
    Integer type     ;
} HarmonicBond ;

typedef struct {
    Real eq ;
    Real fc ;
} HarmonicBondParameter ;

typedef struct {
    Boolean                isSorted    ;
    Integer                nParameters ;
    Integer                nTerms      ;
    HarmonicBond          *terms       ;
    HarmonicBondParameter *parameters  ;
} HarmonicBondContainer ;

/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
extern void                   HarmonicBondContainer_ActivateTerms         (       HarmonicBondContainer  *self ) ;
extern HarmonicBondContainer *HarmonicBondContainer_Allocate              ( const Integer nTerms, const Integer nParameters ) ;
extern HarmonicBondContainer *HarmonicBondContainer_Clone                 ( const HarmonicBondContainer  *self ) ;
extern void                   HarmonicBondContainer_DeactivateTerms       (       HarmonicBondContainer  *self, Selection *selection ) ;
extern void                   HarmonicBondContainer_Deallocate            (       HarmonicBondContainer **self ) ;
extern Real                   HarmonicBondContainer_Energy                ( const HarmonicBondContainer  *self, const Coordinates3 *coordinates3, Coordinates3 *gradients3 ) ;
extern Integer                HarmonicBondContainer_IdentifyBoundaryAtoms (       HarmonicBondContainer  *self, Selection *qcAtoms, Integer **mmboundary, Integer **qcpartners ) ;
extern HarmonicBondContainer *HarmonicBondContainer_Merge                 ( const HarmonicBondContainer  *self, const HarmonicBondContainer *other, const Integer atomincrement ) ;
extern Integer                HarmonicBondContainer_NumberOfInactiveTerms ( const HarmonicBondContainer  *self ) ;
extern HarmonicBondContainer *HarmonicBondContainer_Prune                 (       HarmonicBondContainer  *self, Selection *selection ) ;
extern void                   HarmonicBondContainer_Sort                  (       HarmonicBondContainer  *self ) ;
extern Integer                HarmonicBondContainer_UpperBound            (       HarmonicBondContainer  *self ) ;

# endif
