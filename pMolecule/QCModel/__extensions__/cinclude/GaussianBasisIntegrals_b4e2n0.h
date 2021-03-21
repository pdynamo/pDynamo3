# ifndef _GAUSSIANBASISINTEGRALS_B4E2N0
# define _GAUSSIANBASISINTEGRALS_B4E2N0

# include "BlockStorage.h"
# include "GaussianBasis.h"
# include "Real.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void GaussianBasisIntegrals_TEIs  ( const GaussianBasis *iBasis     ,
                                           const Real          *rI         ,
                                           const GaussianBasis *jBasis     ,
                                           const Real          *rJ         ,
                                           const Real          *rIJ        ,
                                           const Real           rIJ2       ,
                                           const GaussianBasis *kBasis     ,
                                           const Real          *rK         ,
                                           const GaussianBasis *lBasis     ,
                                           const Real          *rL         ,
                                           const Real          *rKL        ,
                                           const Real           rKL2       ,
                                           const Boolean        jLessThanL ,
                                                 Block         *block      ) ;
extern void GaussianBasisIntegrals_TEIsD ( const GaussianBasis *iBasis     ,
                                           const Real          *rI         ,
                                           const GaussianBasis *jBasis     ,
                                           const Real          *rJ         ,
                                           const Real          *rIJ        ,
                                           const Real           rIJ2       ,
                                           const GaussianBasis *kBasis     ,
                                           const Real          *rK         ,
                                           const GaussianBasis *lBasis     ,
                                           const Real          *rL         ,
                                           const Real          *rKL        ,
                                           const Real           rKL2       ,
                                           const Boolean        jLessThanL ,
                                                 Block         *block      ) ;
# endif
