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
extern void Loewdin_AtomicCharges            ( const IntegerArray1D  *basisIndices        ,
                                               const RealArray2D     *loewdinT            ,
                                               const SymmetricMatrix *density             ,
                                                     RealArray1D     *charges             ) ;
extern void Loewdin_BondOrders               ( const IntegerArray1D  *basisIndices        ,
                                               const RealArray2D     *loewdinT            ,
                                               const SymmetricMatrix *density             ,
                                                     SymmetricMatrix *bondOrders          ,
                                                     Status          *status              ) ;
extern void Loewdin_ChargeDensityDerivatives ( const IntegerArray1D  *basisIndices        ,
                                               const RealArray1D     *potentials          ,
                                               const RealArray2D     *loewdinT            ,
                                                     SymmetricMatrix *fock                ) ;
extern void Loewdin_WeightedDensity          ( const IntegerArray1D  *basisIndices        ,
                                               const RealArray1D     *potentials          ,
                                               const RealArray1D     *eigenValues         ,
                                               const RealArray2D     *eigenVectors        ,
                                               const RealArray2D     *loewdinT            ,
                                               const RealArray2D     *o2C                 ,
                                               const RealArray2D     *c2O                 ,
                                               const SymmetricMatrix *density             ,
                                                     Real            *eigenValueTolerance ,
                                                     SymmetricMatrix *wDensity            ,
                                                     Status          *status              ) ;
# endif
