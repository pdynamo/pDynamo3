# ifndef _GAUSSIANBASISORTHONORMALIZE
# define _GAUSSIANBASISORTHONORMALIZE

# include "GaussianBasis.h"
# include "Integer.h"
# include "OrthogonalizingTransformation.h"
# include "Real.h"
# include "RealArray2D.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void GaussianBasis_Orthonormalize (       GaussianBasis          *self         ,
                                           const GaussianBasisOperator   operator     ,
                                           const OrthogonalizationMethod method       ,
                                                 Integer                *nIndependent ,
                                                 Real                   *deviation    ,
                                                 RealArray2D            *MOut         ,
                                                 RealArray2D            *XOut         ,
                                                 RealArray2D            *YOut         ,
                                                 Status                 *status       ) ;
# endif
