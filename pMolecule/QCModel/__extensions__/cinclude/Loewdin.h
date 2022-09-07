# ifndef _LOEWDIN
# define _LOEWDIN

# include "IntegerArray1D.h"
# include "Real.h"
# include "RealArray1D.h"
# include "RealArray2D.h"
# include "Status.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void Loewdin_AtomicCharges                  ( const IntegerArray1D  *basisIndices        ,
                                                     const SymmetricMatrix *loewdinT            ,
                                                     const SymmetricMatrix *density             ,
                                                           RealArray1D     *charges             ,
                                                           Status          *status              ) ;
extern void Loewdin_BondOrders                     ( const IntegerArray1D  *basisIndices        ,
                                                     const SymmetricMatrix *loewdinT            ,
                                                     const SymmetricMatrix *density             ,
                                                           SymmetricMatrix *bondOrders          ,
                                                           Status          *status              ) ;
extern void Loewdin_ChargeDensityDerivatives       ( const IntegerArray1D  *basisIndices        ,
                                                     const RealArray1D     *potentials          ,
                                                     const SymmetricMatrix *loewdinT            ,
                                                           SymmetricMatrix *fock                ) ;
extern Real Loewdin_ChargeRestraintMatrix          ( const IntegerArray1D  *basisIndices        ,
                                                     const RealArray1D     *nuclearCharges      ,
                                                     const IntegerArray1D  *crIndices           ,
                                                     const RealArray1D     *crWeights           ,
                                                     const Boolean          isSpin              ,
                                                     const SymmetricMatrix *loewdinT            ,
                                                           SymmetricMatrix *W                   ) ;
extern void Loewdin_ChargeRestraintWeightedDensity ( const IntegerArray1D  *basisIndices        ,
                                                     const IntegerArray1D  *crIndices           ,
                                                     const RealArray1D     *crWeights           ,
                                                     const Boolean          isSpin              ,
                                                     const Real             dRdL                ,
                                                     const RealArray2D     *eigenVectors        ,
                                                     const RealArray2D     *Z                   ,
                                                           SymmetricMatrix *A                   ) ;
extern void Loewdin_WeightedDensity                ( const IntegerArray1D  *basisIndices        ,
                                                     const RealArray1D     *potentials          ,
                                                     const RealArray1D     *eigenValues         ,
                                                     const RealArray2D     *eigenVectors        ,
                                                     const SymmetricMatrix *loewdinT            ,
                                                     const SymmetricMatrix *density             ,
                                                           Real            *eigenValueTolerance ,
                                                           SymmetricMatrix *wDensity            ,
                                                           Status          *status              ) ;
# endif
