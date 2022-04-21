# ifndef _GAUSSIANBASISINTEGRALS_F2CP1
# define _GAUSSIANBASISINTEGRALS_F2CP1

# include "Coordinates3.h"
# include "GaussianBasis.h"
# include "Real.h"
# include "RealArray1D.h"
# include "RealArray2D.h"
# include "Selection.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void GaussianBasisIntegrals_f2Cm1R1 ( const GaussianBasis *self       ,
                                             const Real          *rS         ,
                                             const GaussianBasis *other      ,
                                             const Real          *rO         ,
                                             const RealArray1D   *charges    ,
                                             const RealArray1D   *widthsE    ,
                                             const RealArray1D   *widthsN    ,
                                             const Coordinates3  *rN         ,
                                             const Selection     *selectionN ,
                                             const RealArray2D   *dOneIJ     ,
                                             const Integer        s2         ,
                                                   Integer       *iWork      ,
                                                   Real          *rWork      ,
                                                   Real          *dRi        ,
                                                   Real          *dRj        ,
                                                   Coordinates3  *gN         ) ;
extern void GaussianBasisIntegrals_f2Cm1V  ( const GaussianBasis *self       ,
                                             const Real          *rS         ,
                                             const GaussianBasis *other      ,
                                             const Real          *rO         ,
                                             const RealArray1D   *charges    ,
                                             const RealArray1D   *widthsE    ,
                                             const RealArray1D   *widthsN    ,
                                             const Coordinates3  *rN         ,
                                             const Selection     *selectionN ,
                                             const Integer        s2         ,
                                                   Integer       *iWork      ,
                                                   Real          *rWork      ,
                                                   RealArray2D   *integrals  ) ;
extern void GaussianBasisIntegrals_f2Cp1V  ( const GaussianBasis *self       ,
                                             const Real          *rS         ,
                                             const GaussianBasis *other      ,
                                             const Real          *rO         ,
                                             const RealArray1D   *widthsE    ,
                                             const RealArray1D   *widthsN    ,
                                             const Coordinates3  *rN         ,
                                             const Selection     *selectionN ,
                                             const RealArray2D   *dOneIJ     ,
                                             const Integer        s2         ,
                                                   Integer       *iWork      ,
                                                   Real          *rWork      ,
                                                   RealArray1D   *potentials ) ;
# endif
