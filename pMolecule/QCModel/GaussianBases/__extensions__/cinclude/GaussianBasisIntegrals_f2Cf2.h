# ifndef _GAUSSIANBASISINTEGRALS_F2CF2
# define _GAUSSIANBASISINTEGRALS_F2CF2

# include "BlockStorage.h"
# include "GaussianBasis.h"
# include "Real.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void GaussianBasisIntegrals_f2Cf2i  ( const GaussianBasis *iBasis     ,
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
                                             const Integer        s4         ,
                                                   Integer       *iWork      ,
                                                   Real          *rWork      ,
                                                   Block         *block      ) ;
extern void GaussianBasisIntegrals_f2Cf2r1 ( const GaussianBasis *iBasis     ,
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
                                             const Integer        s4         ,
                                                   Integer       *iWork      ,
                                                   Real          *rWork      ,
                                                   Block         *block      ) ;
# endif
