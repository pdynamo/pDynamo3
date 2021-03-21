# ifndef _GAUSSIANBASISCONTAINERINTEGRALS_B1E0N1
# define _GAUSSIANBASISCONTAINERINTEGRALS_B1E0N1

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
extern void GaussianBasisContainerIntegrals_Grid                  ( const GaussianBasisContainer *self         ,
                                                                    const IntegerArray1D         *basisIndices ,
                                                                    const Coordinates3           *coordinates3 ,
                                                                    const Coordinates3           *rG           ,
                                                                          RealArray2D            *values       ,
                                                                          Status                 *status       ) ;
extern void GaussianBasisContainerIntegrals_GridFunctionDataBlock ( const GaussianBasisContainer *self         ,
                                                                    const IntegerArray1D         *basisIndices ,
                                                                    const Coordinates3           *coordinates3 ,
                                                                    const Coordinates3           *rG           ,
                                                                    const Boolean                 resize       ,
                                                                    const Real                   *tolerance    ,
                                                                          GridFunctionDataBlock  *data         ,
                                                                          Status                 *status       ) ;
# endif
