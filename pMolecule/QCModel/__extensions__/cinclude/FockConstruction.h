# ifndef _FOCKCONSTRUCTION
# define _FOCKCONSTRUCTION

# include "BlockStorage.h"
# include "Real.h"
# include "RealArray1D.h"
# include "Status.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Real Fock_MakeFromFitIntegrals (       BlockStorage    *fitIntegrals         ,
                                              SymmetricMatrix *fitMatrix            ,
                                        const Real             totalCharge          ,
                                              RealArray1D     *fitPotential         ,
                                              SymmetricMatrix *dTotal               ,
                                              SymmetricMatrix *fTotal               ,
                                              Status          *status               ) ;
extern Real Fock_MakeFromTEIs         (       BlockStorage    *twoElectronIntegrals ,
                                        const SymmetricMatrix *dTotal               ,
                                        const SymmetricMatrix *dSpin                ,
                                        const Real             exchangeScaling      ,
                                              SymmetricMatrix *fTotal               ,
                                              SymmetricMatrix *fSpin                ) ;
extern Real Fock_MakeFromTEIsCoulomb  (       BlockStorage    *twoElectronIntegrals ,
                                        const SymmetricMatrix *dTotal               ,
                                              SymmetricMatrix *fTotal               ) ;
extern Real Fock_MakeFromTEIsExchange (       BlockStorage    *twoElectronIntegrals ,
                                        const SymmetricMatrix *dTotal               ,
                                        const SymmetricMatrix *dSpin                ,
                                        const Real             exchangeScaling      ,
                                              SymmetricMatrix *fTotal               ,
                                              SymmetricMatrix *fSpin                ) ;
# endif
