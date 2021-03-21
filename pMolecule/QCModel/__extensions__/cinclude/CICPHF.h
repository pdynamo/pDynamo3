# ifndef _CICPHF
# define _CICPHF

# include "BlockStorage.h"
# include "Boolean.h"
# include "DoubleSymmetricMatrix.h"
# include "Integer.h"
# include "IntegerArray2D.h"
# include "RealArray1D.h"
# include "RealArray2D.h"
# include "RealArrayND.h"
# include "Status.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void         CICPHF_ApplyCPHFMatrix      ( const Integer                n1                        ,
                                                  const IntegerArray2D        *in1                       ,
                                                  const Integer                n2                        ,
                                                  const IntegerArray2D        *in2                       ,
                                                  const RealArray1D           *aDiagonal                 ,
                                                  const RealArray1D           *b                         ,
                                                  const RealArray2D           *orbitals                  ,
                                                        BlockStorage          *twoElectronIntegrals      ,
                                                        SymmetricMatrix       *work1                     ,
                                                        SymmetricMatrix       *work2                     ,
                                                        RealArray1D           *x                         ) ;
extern void         CICPHF_CalculateCPHFVectors ( const Integer                nActive                   ,
                                                  const Integer                nCore                     ,
                                                  const Integer                nOrbitals                 ,
                                                        BlockStorage          *twoElectronIntegrals      ,
                                                  const DoubleSymmetricMatrix *twoPDM                    ,
                                                  const RealArray1D           *energies                  ,
                                                  const RealArray1D           *occupancies               ,
                                                  const RealArray2D           *orbitals                  ,
                                                  const RealArrayND           *moTEI234                  ,
                                                  const SymmetricMatrix       *fCore                     ,
                                                  const SymmetricMatrix       *onePDM                    ,
                                                  const SymmetricMatrix       *onePDMMO                  ,
                                                        SymmetricMatrix       *work1                     ,
                                                        SymmetricMatrix       *work2                     ,
                                                        Integer               *numberDegenerateRedundant ,
                                                        Integer               *numberNonRedundant        ,
                                                        Integer               *numberRedundant           ,
                                                        IntegerArray2D        *indicesNR                 ,
                                                        IntegerArray2D        *indicesR                  ,
                                                        RealArray1D           *aDiagonal                 ,
                                                        RealArray1D           *qNR                       ,
                                                        RealArray1D           *qR                        ,
                                                        RealArray1D           *preconditioner            ,
                                                        Status                *status                    ) ;
extern RealArray2D *CICPHF_CalculateKPA         ( const Integer                nActive                   ,
                                                  const Integer                nBasis                    ,
                                                  const DoubleSymmetricMatrix *twoPDM                    ,
                                                  const RealArray2D           *orbitals                  ,
                                                  const RealArrayND           *moTEI234                  ,
                                                        Status                *status                    ) ;
extern void         CICPHF_Transform            ( const Integer                n1                        ,
                                                  const IntegerArray2D        *in1                       ,
                                                  const RealArray1D           *x1                        ,
                                                  const Integer                n2                        ,
                                                  const IntegerArray2D        *in2                       ,
                                                  const RealArray1D           *x2                        ,
                                                  const RealArray2D           *orbitals                  ,
                                                  const Boolean                doScale                   ,
                                                        SymmetricMatrix       *work                      ,
                                                        SymmetricMatrix       *z                         ) ;
# endif
