# ifndef _CICONFIGURATIONCONTAINER
# define _CICONFIGURATIONCONTAINER

# include "Boolean.h"
# include "BooleanArray1D.h"
# include "DoubleSymmetricMatrix.h"
# include "Integer.h"
# include "IntegerArray1D.h"
# include "IntegerArray2D.h"
# include "Real.h"
# include "RealArray2D.h"
# include "SparseSymmetricMatrix.h"
# include "Status.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The CI configuration type. */
typedef struct {
    Integer          nAlphas ;
    Integer          nSPQR   ;
    Real             spin    ;
    BooleanArray1D  *parity  ;
    IntegerArray1D  *alphas  ;
    IntegerArray1D  *betas   ;
    IntegerArray1D  *spqr    ;
} CIConfiguration ;

/* . The CI configurations type. */
typedef struct {
    Integer          nActive         ;
    Integer          nConfigurations ;
    Integer          nElectrons      ;
    CIConfiguration *configurations  ;
} CIConfigurationContainer ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern CIConfigurationContainer *CIConfigurationContainer_Allocate                ( const Integer                    nActive               ,
                                                                                    const Integer                    nConfigurations       ,
                                                                                          Status                    *status                ) ;
extern void                      CIConfigurationContainer_Characters              ( const CIConfigurationContainer  *self                  ,
                                                                                    const Boolean                    includeCoreOrbitals   ,
                                                                                    const Integer                    coreOrbitals          ,
                                                                                    const RealArray2D               *orbitalTransformation ,
                                                                                          RealArray2D               *stateTransformation   ,
                                                                                          Status                    *status                ) ;
extern CIConfigurationContainer *CIConfigurationContainer_Clone                   ( const CIConfigurationContainer  *self                  ,
                                                                                          Status                    *status                ) ;
extern void                      CIConfigurationContainer_Deallocate              (       CIConfigurationContainer **self                  ) ;
extern void                      CIConfigurationContainer_GetCIMatrixSparsity     (       CIConfigurationContainer  *self                  ,
                                                                                          Integer                   *nonZero               ,
                                                                                          Real                      *sparsity              ) ;
extern void                      CIConfigurationContainer_MakeCIMatrix            ( const CIConfigurationContainer  *self                  ,
                                                                                    const SymmetricMatrix           *fCoreMO               ,
                                                                                    const DoubleSymmetricMatrix     *moTEIs                ,
                                                                                          SymmetricMatrix           *ciMatrixFull          ,
                                                                                          SparseSymmetricMatrix     *ciMatrixSparse        ,
                                                                                          Status                    *status                ) ;
extern void                      CIConfigurationContainer_MakeDensities           ( const CIConfigurationContainer  *self                  ,
                                                                                    const RealArray1D               *ciVector              ,
                                                                                          SymmetricMatrix           *onePDMMOa             ,
                                                                                          SymmetricMatrix           *onePDMMOb             ,
                                                                                          DoubleSymmetricMatrix     *twoPDM                ,
                                                                                          Status                    *status                ) ;
extern CIConfigurationContainer *CIConfigurationContainer_MakeFull                ( const Integer                    nActive               ,
                                                                                    const Integer                    nUp                   ,
                                                                                    const Integer                    nDown                 ,
                                                                                          Status                    *status                ) ;
extern CIConfigurationContainer *CIConfigurationContainer_MakeSinglesDoubles      ( const Boolean                    doSingles             ,
                                                                                    const Boolean                    doDoubles             ,
                                                                                    const Integer                    nActive               ,
                                                                                    const Integer                    nClosed               ,
                                                                                    const Integer                    nOpen                 ,
                                                                                          Status                    *status                ) ;
extern CIConfigurationContainer *CIConfigurationContainer_MakeUserSpecified       ( const IntegerArray2D            *microStates           ,
                                                                                    const Integer                    activeOrbitals        ,
                                                                                    const Integer                    activeElectrons       ,
                                                                                          Status                    *status                ) ;
extern Integer                   CIConfigurationContainer_NumberOfActiveElectrons ( const CIConfigurationContainer  *self                  ) ;
extern Integer                   CIConfigurationContainer_NumberOfActiveOrbitals  ( const CIConfigurationContainer  *self                  ) ;
extern Integer                   CIConfigurationContainer_NumberOfConfigurations  ( const CIConfigurationContainer  *self                  ) ;
extern void                      CIConfigurationContainer_StateSpins              ( const CIConfigurationContainer  *self                  ,
                                                                                    const RealArray2D               *vectors               ,
                                                                                          RealArray1D               *spins                 ,
                                                                                          Status                    *status                ) ;
extern void                      CIConfigurationContainer_TransitionDipoles       ( const CIConfigurationContainer  *self                  ,
                                                                                          SymmetricMatrix           *tdMOs                 ,
                                                                                          SymmetricMatrix           *tdMatrix              ) ;
# endif
