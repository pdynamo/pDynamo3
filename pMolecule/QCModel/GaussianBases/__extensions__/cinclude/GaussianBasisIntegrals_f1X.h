# ifndef _GAUSSIANBASISINTEGRALS_F1X
# define _GAUSSIANBASISINTEGRALS_F1X

# include "GaussianBasis.h"
# include "RealArray1D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void GaussianBasisIntegrals_f1Di ( const GaussianBasis *iBasis  ,
                                          const Real          *rI      ,
                                          const Real          *rC      ,
                                          const Integer        s1      ,
                                                Real          *rWork   ,
                                                RealArray1D   *dX      ,
                                                RealArray1D   *dY      ,
                                                RealArray1D   *dZ      ) ;
extern void GaussianBasisIntegrals_f1Oi ( const GaussianBasis *self    ,
                                          const Integer        s1      ,
                                                Real          *rWork   ,
                                                RealArray1D   *overlap ) ;
extern void GaussianBasisIntegrals_f1Qi ( const GaussianBasis *iBasis  ,
                                          const Real          *rI      ,
                                          const Real          *rC      ,
                                          const Integer        s1      ,
                                                Real          *rWork   ,
                                                RealArray1D   *qXX     ,
                                                RealArray1D   *qYY     ,
                                                RealArray1D   *qZZ     ,
                                                RealArray1D   *qXY     ,
                                                RealArray1D   *qXZ     ,
                                                RealArray1D   *qYZ     ) ;
# endif
