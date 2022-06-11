# ifndef _GAUSSIANBASISCONTAINERINTEGRALS_F2XF2
# define _GAUSSIANBASISCONTAINERINTEGRALS_F2XF2

# include "BlockStorage.h"
# include "Coordinates3.h"
# include "GaussianBasisContainer.h"
# include "Status.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void GaussianBasisContainerIntegrals_f2Cf2i  ( const GaussianBasisContainer *self            ,
                                                      const Coordinates3           *coordinates3    ,
                                                            BlockStorage           *teis            ,
                                                            Status                 *status          ) ;
extern void GaussianBasisContainerIntegrals_f2Cf2R1 ( const GaussianBasisContainer *self            ,
                                                      const Coordinates3           *coordinates3    ,
                                                      const SymmetricMatrix        *dTotal          ,
                                                      const SymmetricMatrix        *dSpin           ,
                                                      const Boolean                 doCoulomb       ,
                                                      const Real                    exchangeScaling ,
                                                            Coordinates3           *gradients3      ,
                                                            Status                 *status          ) ;
extern void GaussianBasisContainerIntegrals_f2Xf2i  ( const GaussianBasisContainer *self            ,
                                                      const Coordinates3           *coordinates3    ,
                                                      const GaussianBasisOperator   operator        ,
                                                            BlockStorage           *teis            ,
                                                            Status                 *status          ) ;
# endif
