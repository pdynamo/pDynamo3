# ifndef _MULLIKEN
# define _MULLIKEN

# include "IntegerArray1D.h"
# include "RealArray1D.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void Mulliken_AtomicCharges                  ( const IntegerArray1D  *basisIndices   ,
                                                      const SymmetricMatrix *density        ,
                                                      const SymmetricMatrix *overlap        ,
                                                            RealArray1D     *charges        ) ;
extern void Mulliken_BondOrders                     ( const IntegerArray1D  *basisIndices   ,
                                                      const SymmetricMatrix *density        ,
                                                      const SymmetricMatrix *overlap        ,
                                                            SymmetricMatrix *bondOrders     ,
                                                            Status          *status         ) ;
extern void Mulliken_ChargeDensityDerivatives       ( const IntegerArray1D  *basisIndices   ,
                                                      const RealArray1D     *potentials     ,
                                                      const SymmetricMatrix *overlap        ,
                                                            SymmetricMatrix *fock           ) ;
extern Real Mulliken_ChargeRestraintMatrix          ( const IntegerArray1D  *basisIndices   ,
                                                      const RealArray1D     *nuclearCharges ,
                                                      const IntegerArray1D  *crIndices      ,
                                                      const RealArray1D     *crWeights      ,
                                                      const Boolean          isSpin         ,
                                                      const SymmetricMatrix *overlap        ,
                                                            SymmetricMatrix *W              ) ;
extern void Mulliken_ChargeRestraintWeightedDensity ( const IntegerArray1D  *basisIndices   ,
                                                      const IntegerArray1D  *crIndices      ,
                                                      const RealArray1D     *crWeights      ,
                                                      const Boolean          isSpin         ,
                                                      const Real             dRdL           ,
                                                      const SymmetricMatrix *density        ,
                                                            SymmetricMatrix *wdm            ) ;
extern void Mulliken_WeightedDensity                ( const IntegerArray1D  *basisIndices   ,
                                                      const RealArray1D     *potentials     ,
                                                      const SymmetricMatrix *density        ,
                                                             SymmetricMatrix *wDensity       ) ;
# endif
