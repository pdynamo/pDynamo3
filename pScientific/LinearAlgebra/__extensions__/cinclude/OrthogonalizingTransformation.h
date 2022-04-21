# ifndef _ORTHOGONALIZINGTRANSFORMATION
# define _ORTHOGONALIZINGTRANSFORMATION

# include "Boolean.h"
# include "Real.h"
# include "RealArray2D.h"
# include "Status.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Enumerations.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Orthogonalization methods. */
typedef enum {
    OrthogonalizationMethod_Canonical = 1 ,
    OrthogonalizationMethod_Diagonal  = 2 ,
    OrthogonalizationMethod_Symmetric = 3
} OrthogonalizationMethod ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Real    CheckOrthogonalization        (       RealArray2D            *transformation      ,
                                                     RealArray2D            *inverse             ,
                                                     Status                 *status              ) ;
extern Integer OrthogonalizingTransformation (       SymmetricMatrix        *S                   ,
                                               const OrthogonalizationMethod method              ,
                                               const Boolean                 preserveInput       ,
                                               const Real                   *diagonalTolerance   ,
                                               const Real                   *eigenValueTolerance ,
                                                     RealArray1D            *eigenValues         ,
                                                     RealArray2D            *eigenVectors        ,
                                                     RealArray2D            *X                   ,
                                                     Status                 *status              ) ;
# endif
