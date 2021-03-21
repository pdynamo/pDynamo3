# ifndef _GAUSSIANBASISINTEGRALS_B2E1N0
# define _GAUSSIANBASISINTEGRALS_B2E1N0

# include "GaussianBasis.h"
# include "Real.h"
# include "RealArray2D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void GaussianBasisIntegrals_2Coulomb         ( const GaussianBasis *self   ,
                                                      const Real          *rS     ,
                                                      const GaussianBasis *other  ,
                                                      const Real          *rO     ,
                                                            RealArray2D   *s      ) ;
extern void GaussianBasisIntegrals_2CoulombD        ( const GaussianBasis *self   ,
                                                      const Real          *rS     ,
                                                      const GaussianBasis *other  ,
                                                      const Real          *rO     ,
                                                            RealArray2D   *sX     ,
                                                            RealArray2D   *sY     ,
                                                            RealArray2D   *sZ     ) ;
extern void GaussianBasisIntegrals_2Overlap         ( const GaussianBasis *self   ,
                                                      const Real          *rS     ,
                                                      const GaussianBasis *other  ,
                                                      const Real          *rO     ,
                                                            RealArray2D   *s      ) ;
extern void GaussianBasisIntegrals_2OverlapD        ( const GaussianBasis *self   ,
                                                      const Real          *rS     ,
                                                      const GaussianBasis *other  ,
                                                      const Real          *rO     ,
                                                            RealArray2D   *sX     ,
                                                            RealArray2D   *sY     ,
                                                            RealArray2D   *sZ     ) ;
extern void GaussianBasisIntegrals_Dipole           ( const GaussianBasis *self   ,
                                                      const Real          *rS     ,
                                                      const GaussianBasis *other  ,
                                                      const Real          *rO     ,
                                                      const Real          *center ,
                                                            RealArray2D   *sX     ,
                                                            RealArray2D   *sY     ,
                                                            RealArray2D   *sZ     ) ;
extern void GaussianBasisIntegrals_Kinetic2Overlap  ( const GaussianBasis *self   ,
                                                      const Real          *rS     ,
                                                      const GaussianBasis *other  ,
                                                      const Real          *rO     ,
                                                            RealArray2D   *s      ,
                                                            RealArray2D   *t      ) ;
extern void GaussianBasisIntegrals_Kinetic2OverlapD ( const GaussianBasis *self   ,
                                                      const Real          *rS     ,
                                                      const GaussianBasis *other  ,
                                                      const Real          *rO     ,
                                                            RealArray2D   *sX     ,
                                                            RealArray2D   *sY     ,
                                                            RealArray2D   *sZ     ,
                                                            RealArray2D   *tX     ,
                                                            RealArray2D   *tY     ,
                                                            RealArray2D   *tZ     ) ;
# endif
