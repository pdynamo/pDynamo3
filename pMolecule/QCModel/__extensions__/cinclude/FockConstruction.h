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
extern void Fock_MakeCoefficientsFromFitIntegrals (       SymmetricMatrix *fitMatrix            ,
                                                          BlockStorage    *fitIntegrals         ,
                                                          SymmetricMatrix *dTotal               ,
                                                    const Real             totalCharge          ,
                                                          RealArray1D     *fitCoefficients      ,
                                                          RealArray1D     *bVector              ,
                                                          Status          *status               ) ;
extern void Fock_MakeFockFromFitIntegrals         (       BlockStorage    *fitIntegrals         ,
                                                    const RealArray1D     *fitVector            ,
                                                          SymmetricMatrix *fTotal               ,
                                                          Status          *status               ) ;
extern Real Fock_MakeFromFitIntegralsCoulomb      (       BlockStorage    *fitIntegrals         ,
                                                          SymmetricMatrix *fitMatrix            ,
                                                    const Real             totalCharge          ,
                                                          RealArray1D     *fitCoefficients      ,
                                                          SymmetricMatrix *dTotal               ,
                                                          SymmetricMatrix *fTotal               ,
                                                          Status          *status               ) ;
extern Real Fock_MakeFromFitIntegralsNonCoulomb   (       BlockStorage    *fitIntegrals         ,
                                                          SymmetricMatrix *fitMatrix            ,
                                                          SymmetricMatrix *fitCoulombMatrix     ,
                                                    const Real             totalCharge          ,
                                                          RealArray1D     *fitCoefficients      ,
                                                          RealArray1D     *fitVectorD           ,
                                                          SymmetricMatrix *dTotal               ,
                                                          SymmetricMatrix *fTotal               ,
                                                          Status          *status               ) ;
extern Real Fock_MakeFromTEIs                     (       BlockStorage    *twoElectronIntegrals ,
                                                    const SymmetricMatrix *dTotal               ,
                                                    const SymmetricMatrix *dSpin                ,
                                                    const Real             exchangeScaling      ,
                                                          SymmetricMatrix *fTotal               ,
                                                          SymmetricMatrix *fSpin                ) ;
extern Real Fock_MakeFromTEIsCoulomb              (       BlockStorage    *twoElectronIntegrals ,
                                                    const SymmetricMatrix *dTotal               ,
                                                          SymmetricMatrix *fTotal               ) ;
extern Real Fock_MakeFromTEIsExchange             (       BlockStorage    *twoElectronIntegrals ,
                                                    const SymmetricMatrix *dTotal               ,
                                                    const SymmetricMatrix *dSpin                ,
                                                    const Real             exchangeScaling      ,
                                                          SymmetricMatrix *fTotal               ,
                                                          SymmetricMatrix *fSpin                ) ;
# endif
