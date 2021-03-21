# ifndef _GAUSSIANBASISINTEGRALS_B1E0N1
# define _GAUSSIANBASISINTEGRALS_B1E0N1

# include "Coordinates3.h"
# include "GaussianBasis.h"
# include "Real.h"
# include "RealArray2D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void GaussianBasisIntegrals_Grid   ( const GaussianBasis *iBasis ,
                                            const Real          *rI     ,
                                            const Coordinates3  *rG     ,
                                                  RealArray2D   *f      ) ;
extern void GaussianBasisIntegrals_GridD  ( const GaussianBasis *iBasis ,
                                            const Real          *rI     ,
                                            const Coordinates3  *rG     ,
                                                  RealArray2D   *f      ,
                                                  RealArray2D   *fX     ,
                                                  RealArray2D   *fY     ,
                                                  RealArray2D   *fZ     ) ;
extern void GaussianBasisIntegrals_GridD2 ( const GaussianBasis *iBasis ,
                                            const Real          *rI     ,
                                            const Coordinates3  *rG     ,
                                                  RealArray2D   *f      ,
                                                  RealArray2D   *fX     ,
                                                  RealArray2D   *fY     ,
                                                  RealArray2D   *fZ     ,
                                                  RealArray2D   *fXX    ,
                                                  RealArray2D   *fXY    ,
                                                  RealArray2D   *fXZ    ,
                                                  RealArray2D   *fYY    ,
                                                  RealArray2D   *fYZ    ,
                                                  RealArray2D   *fZZ    ) ;
extern void GaussianBasisIntegrals_GridD3 ( const GaussianBasis *iBasis ,
                                            const Real          *rI     ,
                                            const Coordinates3  *rG     ,
                                                  RealArray2D   *f      ,
                                                  RealArray2D   *fX     ,
                                                  RealArray2D   *fY     ,
                                                  RealArray2D   *fZ     ,
                                                  RealArray2D   *fXX    ,
                                                  RealArray2D   *fXY    ,
                                                  RealArray2D   *fXZ    ,
                                                  RealArray2D   *fYY    ,
                                                  RealArray2D   *fYZ    ,
                                                  RealArray2D   *fZZ    ,
                                                  RealArray2D   *fXXX   ,
                                                  RealArray2D   *fXXY   ,
                                                  RealArray2D   *fXXZ   ,
                                                  RealArray2D   *fXYY   ,
                                                  RealArray2D   *fXYZ   ,
                                                  RealArray2D   *fXZZ   ,
                                                  RealArray2D   *fYYY   ,
                                                  RealArray2D   *fYYZ   ,
                                                  RealArray2D   *fYZZ   ,
                                                  RealArray2D   *fZZZ   ) ;
# endif
