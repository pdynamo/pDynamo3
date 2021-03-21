# ifndef _GAUSSIANBASISCONTAINERINTEGRALS_B2E1N1
# define _GAUSSIANBASISCONTAINERINTEGRALS_B2E1N1

# include "Coordinates3.h"
# include "GaussianBasisContainer.h"
# include "IntegerArray1D.h"
# include "RealArray1D.h"
# include "Selection.h"
# include "Status.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void GaussianBasisContainerIntegrals_ElectronNuclear           ( const GaussianBasisContainer *self              ,
                                                                        const IntegerArray1D         *basisIndices      ,
                                                                        const RealArray1D            *charges           ,
                                                                        const RealArray1D            *widthsE           ,
                                                                        const RealArray1D            *widthsN           ,
                                                                        const Coordinates3           *coordinates3      ,
                                                                        const Coordinates3           *coordinates3G     ,
                                                                              Selection              *selectionG        ,
                                                                              SymmetricMatrix        *oneElectronMatrix ,
                                                                              Status                 *status            ) ;
extern void GaussianBasisContainerIntegrals_ElectronNuclearD          ( const GaussianBasisContainer *self              ,
                                                                        const IntegerArray1D         *basisIndices      ,
                                                                        const RealArray1D            *charges           ,
                                                                        const RealArray1D            *widthsE           ,
                                                                        const RealArray1D            *widthsN           ,
                                                                        const Coordinates3           *coordinates3      ,
                                                                        const Coordinates3           *coordinates3G     ,
                                                                              Selection              *selectionG        ,
                                                                        const SymmetricMatrix        *density           ,
                                                                              Coordinates3           *gradients3        ,
                                                                              Coordinates3           *gradients3G       ,
                                                                              Status                 *status            ) ;
extern void GaussianBasisContainerIntegrals_ElectronNuclearPotentials ( const GaussianBasisContainer *self              ,
                                                                        const IntegerArray1D         *basisIndices      ,
                                                                        const RealArray1D            *widthsE           ,
                                                                        const RealArray1D            *widthsN           ,
                                                                        const Coordinates3           *coordinates3      ,
                                                                        const Coordinates3           *coordinates3G     ,
                                                                              Selection              *selectionG        ,
                                                                        const SymmetricMatrix        *density           ,
                                                                              RealArray1D            *potentials        ,
                                                                              Status                 *status            ) ;
# endif
