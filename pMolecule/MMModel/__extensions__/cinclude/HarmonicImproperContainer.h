# ifndef _HARMONICIMPROPERCONTAINER
# define _HARMONICIMPROPERCONTAINER

# include "Boolean.h"
# include "Coordinates3.h"
# include "Integer.h"
# include "Real.h"
# include "Selection.h"

/*------------------------------------------------------------------------------
! . A harmonic improper energy term is defined as follows:
!
!   E = fc * ( phi - phi0 )^2
!
!-----------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
typedef struct {
    Boolean isActive ;
    Integer atom1    ;
    Integer atom2    ;
    Integer atom3    ;
    Integer atom4    ;
    Integer type     ;
} HarmonicImproper ;

typedef struct {
    Real eq    ;
    Real fc    ;
    Real coseq ;
    Real sineq ;
} HarmonicImproperParameter ;

typedef struct {
    Boolean                    isSorted    ;
    Integer                    nParameters ;
    Integer                    nTerms      ;
    HarmonicImproper          *terms       ;
    HarmonicImproperParameter *parameters  ;
} HarmonicImproperContainer ;

/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
extern void                       HarmonicImproperContainer_ActivateTerms         (       HarmonicImproperContainer  *self ) ;
extern HarmonicImproperContainer *HarmonicImproperContainer_Allocate              ( const Integer nTerms, const Integer nParameters ) ;
extern HarmonicImproperContainer *HarmonicImproperContainer_Clone                 ( const HarmonicImproperContainer  *self ) ;
extern void                       HarmonicImproperContainer_DeactivateTerms       (       HarmonicImproperContainer  *self, Selection *selection ) ;
extern void                       HarmonicImproperContainer_Deallocate            (       HarmonicImproperContainer **self ) ;
extern Real                       HarmonicImproperContainer_Energy                ( const HarmonicImproperContainer  *self, const Coordinates3 *coordinates3, Coordinates3 *gradients3 ) ;
extern void                       HarmonicImproperContainer_FillCosSinValues      (       HarmonicImproperContainer  *self ) ;
extern HarmonicImproperContainer *HarmonicImproperContainer_Merge                 ( const HarmonicImproperContainer  *self, const HarmonicImproperContainer *other, const Integer atomincrement ) ;
extern Integer                    HarmonicImproperContainer_NumberOfInactiveTerms ( const HarmonicImproperContainer  *self ) ;
extern HarmonicImproperContainer *HarmonicImproperContainer_Prune                 (       HarmonicImproperContainer  *self, Selection *selection ) ;
extern void                       HarmonicImproperContainer_Sort                  (       HarmonicImproperContainer  *self ) ;
extern Integer                    HarmonicImproperContainer_UpperBound            (       HarmonicImproperContainer  *self ) ;

# endif
