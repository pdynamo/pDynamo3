# ifndef _GAUSSIANBASISINTEGRALS_B3E1N0
# define _GAUSSIANBASISINTEGRALS_B3E1N0

# include "BlockStorage.h"
# include "GaussianBasis.h"
# include "Real.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void GaussianBasisIntegrals_ElectronFit  ( const GaussianBasis *iBasis ,
                                                  const Real          *rI     ,
                                                  const GaussianBasis *jBasis ,
                                                  const Real          *rJ     ,
                                                  const Real          *rIJ    ,
                                                  const Real           rIJ2   ,
                                                  const GaussianBasis *fBasis ,
                                                  const Real          *rF     ,
                                                        Block         *block  ) ;
extern void GaussianBasisIntegrals_ElectronFitD ( const GaussianBasis *iBasis ,
                                                  const Real          *rI     ,
                                                  const GaussianBasis *jBasis ,
                                                  const Real          *rJ     ,
                                                  const Real          *rIJ    ,
                                                  const Real           rIJ2   ,
                                                  const GaussianBasis *fBasis ,
                                                  const Real          *rF     ,
                                                        Block         *block  ) ;
# endif
