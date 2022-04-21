# ifndef _GAUSSIANBASISCONTAINER
# define _GAUSSIANBASISCONTAINER

# include "Boolean.h"
# include "GaussianBasis.h"
# include "Integer.h"
# include "IntegerArray1D.h"
# include "RealArray2D.h"
# include "Status.h"

/*# define _MakeC2S_*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The container type. */
typedef struct {
    Boolean         isOwner                ;
    Integer         capacity               ;
    IntegerArray1D *centerFunctionPointers ;
    IntegerArray1D *functionCenters        ;
    GaussianBasis **entries                ;
} GaussianBasisContainer ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern GaussianBasisContainer *GaussianBasisContainer_Allocate              ( const Integer                  capacity               ,
                                                                                    Status                  *status                 ) ;
extern GaussianBasisContainer *GaussianBasisContainer_Clone                 ( const GaussianBasisContainer  *self                   ,
                                                                                    Status                  *status                 ) ;
extern void                    GaussianBasisContainer_Deallocate            (       GaussianBasisContainer **self                   ) ;
extern Integer                 GaussianBasisContainer_LargestBasis          ( const GaussianBasisContainer  *self                   ,
                                                                              const Boolean                  forC                   ) ;
extern Integer                 GaussianBasisContainer_LargestShell          ( const GaussianBasisContainer  *self                   ,
                                                                              const Boolean                  forC                   ) ;
# ifdef _MakeC2S_
extern void                    GaussianBasisContainer_MakeC2S               ( const GaussianBasisContainer  *self                   ,
                                                                                    RealArray2D             *T                      ,
                                                                                    Status                  *status                 ) ;
# endif
extern void                    GaussianBasisContainer_MakeIndexArrays       (       GaussianBasisContainer  *self                   ,
                                                                                    IntegerArray1D          *centerFunctionPointers ,
                                                                                    IntegerArray1D          *functionCenters        ,
                                                                                    Status                  *status                 ) ;
extern Integer                 GaussianBasisContainer_NumberOfFunctions     ( const GaussianBasisContainer  *self                   ) ;
extern Integer                 GaussianBasisContainer_NumberOfWorkFunctions ( const GaussianBasisContainer  *self                   ) ;

# endif
