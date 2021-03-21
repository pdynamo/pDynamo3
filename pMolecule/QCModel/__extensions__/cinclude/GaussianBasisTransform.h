# ifndef _GAUSSIANBASISTRANSFORM
# define _GAUSSIANBASISTRANSFORM

# include "RealArray2D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void GaussianBasis_TransformIntegrals2 (       RealArray2D **integrals ,
                                                const RealArray2D  *ic2o      ,
                                                const RealArray2D  *jc2o      ) ;

# endif
