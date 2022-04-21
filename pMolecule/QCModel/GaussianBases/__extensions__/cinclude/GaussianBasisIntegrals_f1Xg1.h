# ifndef _GAUSSIANBASISINTEGRALS_F1XG1
# define _GAUSSIANBASISINTEGRALS_F1XG1

# include "GaussianBasis.h"
# include "Real.h"
# include "RealArray2D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void GaussianBasisIntegrals_f1Ag1i   ( const GaussianBasis *self   ,
                                              const Real          *rS     ,
                                              const GaussianBasis *other  ,
                                              const Real          *rO     ,
                                              const Integer        s2     ,
                                                    Integer       *iWork  ,
                                                    Real          *rWork  ,
                                                    RealArray2D   *s      ) ;
extern void GaussianBasisIntegrals_f1Ag1r1  ( const GaussianBasis *self   ,
                                              const Real          *rS     ,
                                              const GaussianBasis *other  ,
                                              const Real          *rO     ,
                                              const Integer        s2     ,
                                                    Integer       *iWork  ,
                                                    Real          *rWork  ,
                                                    RealArray2D   *sX     ,
                                                    RealArray2D   *sY     ,
                                                    RealArray2D   *sZ     ) ;
extern void GaussianBasisIntegrals_f1Cg1i   ( const GaussianBasis *self   ,
                                              const Real          *rS     ,
                                              const GaussianBasis *other  ,
                                              const Real          *rO     ,
                                              const Integer        s2     ,
                                                    Integer       *iWork  ,
                                                    Real          *rWork  ,
                                                    RealArray2D   *s      ) ;
extern void GaussianBasisIntegrals_f1Cg1r1  ( const GaussianBasis *self   ,
                                              const Real          *rS     ,
                                              const GaussianBasis *other  ,
                                              const Real          *rO     ,
                                              const Integer        s2     ,
                                                    Integer       *iWork  ,
                                                    Real          *rWork  ,
                                                    RealArray2D   *sX     ,
                                                    RealArray2D   *sY     ,
                                                    RealArray2D   *sZ     ) ;
extern void GaussianBasisIntegrals_f1Df1i   ( const GaussianBasis *self   ,
                                              const Real          *rS     ,
                                              const GaussianBasis *other  ,
                                              const Real          *rO     ,
                                              const Real          *center ,
                                              const Integer        s2     ,
                                                    Real          *rWork  ,
                                                    RealArray2D   *sX     ,
                                                    RealArray2D   *sY     ,
                                                    RealArray2D   *sZ     ) ;
extern void GaussianBasisIntegrals_f1KOg1i  ( const GaussianBasis *self   ,
                                              const Real          *rS     ,
                                              const GaussianBasis *other  ,
                                              const Real          *rO     ,
                                              const Integer        s2     ,
                                                    Real          *rWork  ,
                                                    RealArray2D   *s      ,
                                                    RealArray2D   *t      ) ;
extern void GaussianBasisIntegrals_f1KOg1r1 ( const GaussianBasis *self   ,
                                              const Real          *rS     ,
                                              const GaussianBasis *other  ,
                                              const Real          *rO     ,
                                              const Integer        s2     ,
                                                    Real          *rWork  ,
                                                    RealArray2D   *sX     ,
                                                    RealArray2D   *sY     ,
                                                    RealArray2D   *sZ     ,
                                                    RealArray2D   *tX     ,
                                                    RealArray2D   *tY     ,
                                                    RealArray2D   *tZ     ) ;
extern void GaussianBasisIntegrals_f1Og1i   ( const GaussianBasis *self   ,
                                              const Real          *rS     ,
                                              const GaussianBasis *other  ,
                                              const Real          *rO     ,
                                              const Integer        s2     ,
                                                    Real          *rWork  ,
                                                    RealArray2D   *s      ) ;
extern void GaussianBasisIntegrals_f1Og1r1  ( const GaussianBasis *self   ,
                                              const Real          *rS     ,
                                              const GaussianBasis *other  ,
                                              const Real          *rO     ,
                                              const Integer        s2     ,
                                                    Real          *rWork  ,
                                                    RealArray2D   *sX     ,
                                                    RealArray2D   *sY     ,
                                                    RealArray2D   *sZ     ) ;
extern void GaussianBasisIntegrals_f1Qf1i   ( const GaussianBasis *self   ,
                                              const Real          *rS     ,
                                              const GaussianBasis *other  ,
                                              const Real          *rO     ,
                                              const Real          *center ,
                                              const Integer        s2     ,
                                                    Real          *rWork  ,
                                                    RealArray2D   *sXX    ,
                                                    RealArray2D   *sYY    ,
                                                    RealArray2D   *sZZ    ,
                                                    RealArray2D   *sXY    ,
                                                    RealArray2D   *sXZ    ,
                                                    RealArray2D   *sYZ    ) ;
# endif
