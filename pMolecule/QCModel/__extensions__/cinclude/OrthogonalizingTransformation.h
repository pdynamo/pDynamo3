# ifndef _ORTHOGONALIZINGTRANSFORMATION
# define _ORTHOGONALIZINGTRANSFORMATION

# include "Boolean.h"
# include "Real.h"
# include "RealArray2D.h"
# include "Status.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Real    CheckOrthogonalization        (       RealArray2D     *transformation      ,
                                                     RealArray2D     *inverse             ,
                                                     Status          *status              ) ;
extern Integer OrthogonalizingTransformation (       SymmetricMatrix *S                   ,
                                               const Boolean          doCanonical         ,
                                               const Boolean          preserveInput       ,
                                               const Real            *eigenValueTolerance ,
                                                     RealArray1D     *eigenValues         ,
                                                     RealArray2D     *eigenVectors        ,
                                                     RealArray2D     *X                   ,
                                                     Status          *status              ) ;
# endif
