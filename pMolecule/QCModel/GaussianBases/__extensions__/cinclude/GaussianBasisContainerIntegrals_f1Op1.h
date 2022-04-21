# ifndef _GAUSSIANBASISCONTAINERINTEGRALS_F1OP1
# define _GAUSSIANBASISCONTAINERINTEGRALS_F1OP1

# include "Boolean.h"
# include "Coordinates3.h"
# include "GaussianBasisContainer.h"
# include "GridFunctionDataBlock.h"
# include "IntegerArray1D.h"
# include "Real.h"
# include "RealArray2D.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void GaussianBasisContainerIntegrals_f1Op1i     ( const GaussianBasisContainer *self         ,
                                                         const Coordinates3           *coordinates3 ,
                                                         const Coordinates3           *rG           ,
                                                               RealArray2D            *values       ,
                                                               Status                 *status       ) ;
extern void GaussianBasisContainerIntegrals_f1Op1ir123 ( const GaussianBasisContainer *self         ,
                                                         const Coordinates3           *coordinates3 ,
                                                         const Coordinates3           *rG           ,
                                                         const Boolean                 resize       ,
                                                         const Real                   *tolerance    ,
                                                               GridFunctionDataBlock  *data         ,
                                                               Status                 *status       ) ;
# endif
