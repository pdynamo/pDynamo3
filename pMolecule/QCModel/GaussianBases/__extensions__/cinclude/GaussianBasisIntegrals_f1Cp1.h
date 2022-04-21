# ifndef _GAUSSIANBASISINTEGRALS_F1CP1
# define _GAUSSIANBASISINTEGRALS_F1CP1

# include "Coordinates3.h"
# include "GaussianBasis.h"
# include "Real.h"
# include "RealArray1D.h"
# include "RealArray2D.h"
# include "Selection.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void GaussianBasisIntegrals_f1Cm1R1 ( const GaussianBasis *self       ,
                                             const Real          *rS         ,
                                             const RealArray1D   *charges    ,
                                             const RealArray1D   *widthsE    ,
                                             const RealArray1D   *widthsN    ,
                                             const Coordinates3  *rN         ,
                                             const Selection     *selectionN ,
                                             const RealArray1D   *dOneF      ,
                                             const Integer        s1         ,
                                                   Integer       *iWork      ,
                                                   Real          *rWork      ,
                                                   Real          *gF         ,
                                                   Coordinates3  *gN         ) ;
extern void GaussianBasisIntegrals_f1Cm1V  ( const GaussianBasis *self       ,
                                             const Real          *rS         ,
                                             const RealArray1D   *charges    ,
                                             const RealArray1D   *widthsE    ,
                                             const RealArray1D   *widthsN    ,
                                             const Coordinates3  *rN         ,
                                             const Selection     *selectionN ,
                                             const Integer        s1         ,
                                                   Integer       *iWork      ,
                                                   Real          *rWork      ,
                                                   RealArray1D   *integrals  ) ;
extern void GaussianBasisIntegrals_f1Cp1V  ( const GaussianBasis *self       ,
                                             const Real          *rS         ,
                                             const RealArray1D   *widthsE    ,
                                             const RealArray1D   *widthsN    ,
                                             const Coordinates3  *rN         ,
                                             const Selection     *selectionN ,
                                             const RealArray1D   *dOneF      ,
                                             const Integer        s1         ,
                                                   Integer       *iWork      ,
                                                   Real          *rWork      ,
                                                   RealArray1D   *potentials ) ;
# endif
