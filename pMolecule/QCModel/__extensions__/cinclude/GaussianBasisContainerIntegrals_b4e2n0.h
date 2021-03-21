# ifndef _GAUSSIANBASISCONTAINERINTEGRALS_B4E2N0
# define _GAUSSIANBASISCONTAINERINTEGRALS_B4E2N0

# include "BlockStorage.h"
# include "Coordinates3.h"
# include "GaussianBasisContainer.h"
# include "Status.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void GaussianBasisContainerIntegrals_TEIs  ( const GaussianBasisContainer *self            ,
                                                    const IntegerArray1D         *basisIndices    ,
                                                    const Coordinates3           *coordinates3    ,
                                                          BlockStorage           *teis            ,
                                                          Status                 *status          ) ;
extern void GaussianBasisContainerIntegrals_TEIsD ( const GaussianBasisContainer *self            ,
                                                    const IntegerArray1D         *basisIndices    ,
                                                    const Coordinates3           *coordinates3    ,
                                                    const SymmetricMatrix        *dTotal          ,
                                                    const SymmetricMatrix        *dSpin           ,
                                                    const Boolean                 doCoulomb       ,
                                                    const Real                    exchangeScaling ,
                                                          Coordinates3           *gradients3      ,
                                                          Status                 *status          ) ;
# endif
