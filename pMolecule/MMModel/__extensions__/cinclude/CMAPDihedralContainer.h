# ifndef _CMAPDIHEDRALCONTAINER
# define _CMAPDIHEDRALCONTAINER

# include "BicubicSpline.h"
# include "Boolean.h"
# include "Coordinates3.h"
# include "Integer.h"
# include "Real.h"
# include "Selection.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . A CMAP dihedral term is for CHARMM force fields. It is a table-based energy
! . term that is a function of two dihedral angles.
!---------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
typedef struct {
    Boolean isActive ;
    Integer atom1    ;
    Integer atom2    ;
    Integer atom3    ;
    Integer atom4    ;
    Integer atom5    ;
    Integer atom6    ;
    Integer atom7    ;
    Integer atom8    ;
    Integer type     ;
} CMAPDihedral ;

typedef struct {
    Boolean         isSorted    ;
    Integer         nParameters ;
    Integer         nTerms      ;
    CMAPDihedral   *terms       ;
    BicubicSpline **parameters  ;
} CMAPDihedralContainer ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void                   CMAPDihedralContainer_ActivateTerms         (       CMAPDihedralContainer  *self ) ;
extern CMAPDihedralContainer *CMAPDihedralContainer_Allocate              ( const Integer nTerms, const Integer nParameters ) ;
extern CMAPDihedralContainer *CMAPDihedralContainer_Clone                 ( const CMAPDihedralContainer  *self ) ;
extern void                   CMAPDihedralContainer_DeactivateTerms       (       CMAPDihedralContainer  *self, Selection *selection ) ;
extern void                   CMAPDihedralContainer_Deallocate            (       CMAPDihedralContainer **self ) ;
extern Real                   CMAPDihedralContainer_Energy                ( const CMAPDihedralContainer  *self, const Coordinates3 *coordinates3, Coordinates3 *gradients3 ) ;
extern CMAPDihedralContainer *CMAPDihedralContainer_Merge                 ( const CMAPDihedralContainer  *self, const CMAPDihedralContainer *other, const int atomincrement ) ;
extern Integer                CMAPDihedralContainer_NumberOfInactiveTerms ( const CMAPDihedralContainer  *self ) ;
extern CMAPDihedralContainer *CMAPDihedralContainer_Prune                 (       CMAPDihedralContainer  *self, Selection *selection ) ;
extern void                   CMAPDihedralContainer_Sort                  (       CMAPDihedralContainer  *self ) ;
extern Integer                CMAPDihedralContainer_UpperBound            (       CMAPDihedralContainer  *self ) ;

# endif
