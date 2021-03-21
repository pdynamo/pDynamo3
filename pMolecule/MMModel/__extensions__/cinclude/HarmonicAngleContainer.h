# ifndef _HARMONICANGLECONTAINER
# define _HARMONICANGLECONTAINER

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
    Integer atom3    ;
    Integer type     ;
} HarmonicAngle ;

typedef struct {
    Real eq ;
    Real fc ;
} HarmonicAngleParameter ;

typedef struct {
    Boolean                 isSorted    ;
    Integer                 nParameters ;
    Integer                 nTerms      ;
    HarmonicAngle          *terms       ;
    HarmonicAngleParameter *parameters  ;
} HarmonicAngleContainer ;

/*------------------------------------------------------------------------------
! . Procedure declarations.
!-----------------------------------------------------------------------------*/
extern void                    HarmonicAngleContainer_ActivateTerms         (       HarmonicAngleContainer  *self ) ;
extern HarmonicAngleContainer *HarmonicAngleContainer_Allocate              ( const Integer nTerms, const Integer nParameters ) ;
extern HarmonicAngleContainer *HarmonicAngleContainer_Clone                 ( const HarmonicAngleContainer  *self ) ;
extern void                    HarmonicAngleContainer_DeactivateTerms       (       HarmonicAngleContainer  *self, Selection *selection ) ;
extern void                    HarmonicAngleContainer_Deallocate            (       HarmonicAngleContainer **self ) ;
extern Real                    HarmonicAngleContainer_Energy                ( const HarmonicAngleContainer  *self, const Coordinates3 *coordinates3, Coordinates3 *gradients3 ) ;
extern HarmonicAngleContainer *HarmonicAngleContainer_Merge                 ( const HarmonicAngleContainer  *self, const HarmonicAngleContainer *other, const Integer atomincrement ) ;
extern Integer                 HarmonicAngleContainer_NumberOfInactiveTerms ( const HarmonicAngleContainer  *self ) ;
extern HarmonicAngleContainer *HarmonicAngleContainer_Prune                 (       HarmonicAngleContainer  *self, Selection *selection ) ;
extern void                    HarmonicAngleContainer_Sort                  (       HarmonicAngleContainer  *self ) ;
extern Integer                 HarmonicAngleContainer_UpperBound            (       HarmonicAngleContainer  *self ) ;

# endif
