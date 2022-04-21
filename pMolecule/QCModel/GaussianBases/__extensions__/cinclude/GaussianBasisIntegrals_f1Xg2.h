# ifndef _GAUSSIANBASISINTEGRALS_F1XG2
# define _GAUSSIANBASISINTEGRALS_F1XG2

# include "BlockStorage.h"
# include "GaussianBasis.h"
# include "Real.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void GaussianBasisIntegrals_f1Ag2i  ( const GaussianBasis *iBasis ,
                                             const Real          *rI     ,
                                             const GaussianBasis *jBasis ,
                                             const Real          *rJ     ,
                                             const Real          *rIJ    ,
                                             const Real           rIJ2   ,
                                             const GaussianBasis *fBasis ,
                                             const Real          *rF     ,
                                             const Integer        s3     ,
                                                   Integer       *iWork  ,
                                                   Real          *rWork  ,
                                                   Block         *block  ) ;
extern void GaussianBasisIntegrals_f1Ag2r1 ( const GaussianBasis *iBasis ,
                                             const Real          *rI     ,
                                             const GaussianBasis *jBasis ,
                                             const Real          *rJ     ,
                                             const Real          *rIJ    ,
                                             const Real           rIJ2   ,
                                             const GaussianBasis *fBasis ,
                                             const Real          *rF     ,
                                             const Integer        s3     ,
                                                   Integer       *iWork  ,
                                                   Real          *rWork  ,
                                                   Block         *block  ) ;
extern void GaussianBasisIntegrals_f1Cg2i  ( const GaussianBasis *iBasis ,
                                             const Real          *rI     ,
                                             const GaussianBasis *jBasis ,
                                             const Real          *rJ     ,
                                             const Real          *rIJ    ,
                                             const Real           rIJ2   ,
                                             const GaussianBasis *fBasis ,
                                             const Real          *rF     ,
                                             const Integer        s3     ,
                                                   Integer       *iWork  ,
                                                   Real          *rWork  ,
                                                   Block         *block  ) ;
extern void GaussianBasisIntegrals_f1Cg2r1 ( const GaussianBasis *iBasis ,
                                             const Real          *rI     ,
                                             const GaussianBasis *jBasis ,
                                             const Real          *rJ     ,
                                             const Real          *rIJ    ,
                                             const Real           rIJ2   ,
                                             const GaussianBasis *fBasis ,
                                             const Real          *rF     ,
                                             const Integer        s3     ,
                                                   Integer       *iWork  ,
                                                   Real          *rWork  ,
                                                   Block         *block  ) ;
extern void GaussianBasisIntegrals_f1Og2i  ( const GaussianBasis *iBasis ,
                                             const Real          *rI     ,
                                             const GaussianBasis *jBasis ,
                                             const Real          *rJ     ,
                                             const GaussianBasis *fBasis ,
                                             const Real          *rF     ,
                                             const Integer        s3     ,
                                                   Real          *rWork  ,
                                                   Block         *block  ) ;
extern void GaussianBasisIntegrals_f1Og2r1 ( const GaussianBasis *iBasis ,
                                             const Real          *rI     ,
                                             const GaussianBasis *jBasis ,
                                             const Real          *rJ     ,
                                             const GaussianBasis *fBasis ,
                                             const Real          *rF     ,
                                             const Integer        s3     ,
                                                   Real          *rWork  ,
                                                   Block         *block  ) ;
# endif
