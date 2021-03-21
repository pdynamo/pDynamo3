# ifndef _GAUSSIANBASISCONTAINER
# define _GAUSSIANBASISCONTAINER

# include "Boolean.h"
# include "GaussianBasis.h"
# include "Integer.h"
# include "IntegerArray1D.h"
# include "RealArray2D.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The container type. */
typedef struct {
    Boolean         isOwner  ;
    Integer         capacity ;
    GaussianBasis **entries  ;
} GaussianBasisContainer ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern GaussianBasisContainer *GaussianBasisContainer_Allocate                    ( const Integer                  capacity ,
                                                                                          Status                  *status   ) ;
extern GaussianBasisContainer *GaussianBasisContainer_Clone                       ( const GaussianBasisContainer  *self     ,
                                                                                          Status                  *status   ) ;
extern void                    GaussianBasisContainer_Deallocate                  (       GaussianBasisContainer **self     ) ;
extern Integer                 GaussianBasisContainer_LargestBasis                ( const GaussianBasisContainer  *self     ,
                                                                                    const Integer                  doWork   ) ;
extern void                    GaussianBasisContainer_MakeBasisAtomIndices        ( const GaussianBasisContainer  *self     ,
                                                                                    const Boolean                  doWork   ,
                                                                                          IntegerArray1D          *indices  ,
                                                                                          Status                  *status   ) ;
extern void                    GaussianBasisContainer_MakeBasisIndices            ( const GaussianBasisContainer  *self     ,
                                                                                    const Boolean                  doWork   ,
                                                                                          IntegerArray1D          *indices  ,
                                                                                          Status                  *status   ) ;
extern void                    GaussianBasisContainer_MakeFunctionTransformations ( const GaussianBasisContainer  *self     ,
                                                                                          RealArray2D             *c2o      ,
                                                                                          RealArray2D             *o2c      ,
                                                                                          Status                  *status   ) ;
extern Integer                 GaussianBasisContainer_NumberOfBasisFunctions      ( const GaussianBasisContainer  *self     ,
                                                                                    const Boolean                  doWork   ) ;

# endif
