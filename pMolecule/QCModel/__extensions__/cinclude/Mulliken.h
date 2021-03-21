# ifndef _MULLIKEN
# define _MULLIKEN

# include "IntegerArray1D.h"
# include "RealArray1D.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void Mulliken_AtomicCharges            ( const IntegerArray1D  *basisIndices ,
                                                const SymmetricMatrix *density      ,
                                                const SymmetricMatrix *overlap      ,
                                                      RealArray1D     *charges      ) ;
extern void Mulliken_BondOrders               ( const IntegerArray1D  *basisIndices ,
                                                const SymmetricMatrix *density      ,
                                                const SymmetricMatrix *overlap      ,
                                                      SymmetricMatrix *bondOrders   ,
                                                      Status          *status       ) ;
extern void Mulliken_ChargeDensityDerivatives ( const IntegerArray1D  *basisIndices ,
                                                const RealArray1D     *potentials   ,
                                                const SymmetricMatrix *overlap      ,
                                                      SymmetricMatrix *fock         ) ;
extern void Mulliken_WeightedDensity          ( const IntegerArray1D  *basisIndices ,
                                                const RealArray1D     *potentials   ,
                                                const SymmetricMatrix *density      ,
                                                      SymmetricMatrix *wDensity     ) ;
# endif
