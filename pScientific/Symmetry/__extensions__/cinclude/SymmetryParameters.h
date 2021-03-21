# ifndef _SYMMETRYPARAMETERS
# define _SYMMETRYPARAMETERS

# include "Boolean.h"
# include "Coordinates3.h"
# include "Integer.h"
# include "Matrix33.h"
# include "Real.h"
# include "Selection.h"
# include "SelectionContainer.h"
# include "Status.h"
# include "Vector3.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The symmetry parameters type. */
typedef struct {
    Boolean   isOrthogonal ;
    Real      a            ;
    Real      b            ;
    Real      c            ;
    Real      alpha        ;
    Real      beta         ;
    Real      gamma        ;
    Matrix33 *H            ; /* . f -> r. */
    Matrix33 *inverseH     ; /* . r -> f. */
} SymmetryParameters ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern SymmetryParameters *SymmetryParameters_Allocate                          (       Status              *status       ) ;
extern SymmetryParameters *SymmetryParameters_AllocateFull                      (       Status              *status       ) ;
extern SymmetryParameters *SymmetryParameters_AllocateWithMatrices              ( const Matrix33            *H            ,
                                                                                  const Matrix33            *inverseH     ,
                                                                                        Status              *status       ) ;
extern void                SymmetryParameters_CenterCoordinates3ByFreeIsolate   ( const SymmetryParameters  *self         ,
                                                                                  const SelectionContainer  *isolates     ,
                                                                                  const BooleanBlock        *freeIsolates ,
                                                                                        Coordinates3        *coordinates3 ,
                                                                                        Status              *status       ) ;
extern void                SymmetryParameters_CenterCoordinates3ByIndex         ( const SymmetryParameters  *self         ,
                                                                                  const Selection           *selection    ,
                                                                                        Coordinates3        *coordinates3 ,
                                                                                        Status              *status       ) ;
extern void                SymmetryParameters_CenterCoordinates3ByIsolate       ( const SymmetryParameters  *self         ,
                                                                                  const SelectionContainer  *isolates     ,
                                                                                        Selection           *selection    ,
                                                                                        Coordinates3        *coordinates3 ,
                                                                                        Status              *status       ) ;
extern void                SymmetryParameters_ClearH                            (       SymmetryParameters  *self         ) ;
extern void                SymmetryParameters_CopyTo                            ( const SymmetryParameters  *self         ,
                                                                                        SymmetryParameters  *other        ) ;
extern void                SymmetryParameters_Deallocate                        (       SymmetryParameters **self         ) ;
extern void                SymmetryParameters_Displacement                      ( const SymmetryParameters  *self         ,
                                                                                  const Integer              a            ,
                                                                                  const Integer              b            ,
                                                                                  const Integer              c            ,
                                                                                        Vector3             *displacement ) ;
extern void                SymmetryParameters_FindBoxSearchLimits               ( const SymmetryParameters  *self         ,
                                                                                  const Vector3             *lower        ,
                                                                                  const Vector3             *upper        ,
                                                                                  const Vector3             *iLower       ,
                                                                                  const Vector3             *iUpper       ,
                                                                                        Integer             *aLow         ,
                                                                                        Integer             *aHigh        ,
                                                                                        Integer             *bLow         ,
                                                                                        Integer             *bHigh        ,
                                                                                        Integer             *cLow         ,
                                                                                        Integer             *cHigh        ) ;
extern void                SymmetryParameters_FindCenteringTranslation          ( const SymmetryParameters  *self         ,
                                                                                  const Vector3             *point        ,
                                                                                        Vector3             *translation  ) ;
extern Boolean             SymmetryParameters_IsMinimumImageConventionSatisfied ( const SymmetryParameters  *self         ,
                                                                                  const Real                 length       ) ;
extern Boolean             SymmetryParameters_IsOrthogonal                      ( const SymmetryParameters  *self         ) ;
extern void                SymmetryParameters_IsotropicScale                    (       SymmetryParameters  *self         ,
                                                                                  const Real                 scale        ) ;
extern Coordinates3       *SymmetryParameters_MakeFractionalCoordinates         ( const SymmetryParameters  *self         ,
                                                                                  const Coordinates3        *coordinates3 ,
                                                                                        Status              *status       ) ;
extern void                SymmetryParameters_MakeG                             ( const SymmetryParameters  *self         ,
                                                                                        Matrix33            *G            ) ;
extern void                SymmetryParameters_MakeH                             (       SymmetryParameters  *self         ) ;
extern void                SymmetryParameters_MakeMinimumImageVector            ( const SymmetryParameters  *self         ,
                                                                                        Real                *r            ,
                                                                                        Real                *dR           ) ;
extern void                SymmetryParameters_PerpendicularWidths               ( const SymmetryParameters  *self         ,
                                                                                        Real                *widths       ) ;
extern void                SymmetryParameters_SetCrystalParameters              (       SymmetryParameters  *self         ,
                                                                                  const Real                 a            ,
                                                                                  const Real                 b            ,
                                                                                  const Real                 c            ,
                                                                                  const Real                 alpha        ,
                                                                                  const Real                 beta         ,
                                                                                  const Real                 gamma        ) ;
extern Real                SymmetryParameters_Volume                            ( const SymmetryParameters  *self         ) ;

# endif
