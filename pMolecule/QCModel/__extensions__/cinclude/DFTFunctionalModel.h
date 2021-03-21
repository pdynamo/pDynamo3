# ifndef _DFTFUNCTIONALMODEL
# define _DFTFUNCTIONALMODEL

# include "Boolean.h"
# include "DFTIntegratorDataBlock.h"
# include "Integer.h"
# include "IntegerArray1D.h"
# include "Status.h"
# include "xc.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The DFT functional model type. */
typedef struct {
    Boolean       hasLaplacian        ;
    Boolean       hasSigma            ;
    Boolean       hasTau              ;
    Boolean       isSpinRestricted    ;
    Integer       numberOfFunctionals ;
    Integer       order               ;
    xc_func_type *functionals         ;
} DFTFunctionalModel ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern DFTFunctionalModel *DFTFunctionalModel_Allocate        ( const Integer                 numberOfFunctionals ,
                                                                      Status                 *status              ) ;
extern DFTFunctionalModel *DFTFunctionalModel_Clone           ( const DFTFunctionalModel     *self                ,
                                                                      Status                 *status              ) ;
extern void                DFTFunctionalModel_Deallocate      (       DFTFunctionalModel    **self                ) ;
extern void                DFTFunctionalModel_Evaluate        ( const DFTFunctionalModel     *self                ,
                                                                      DFTIntegratorDataBlock *data                ) ;
extern Real                DFTFunctionalModel_ExchangeScaling ( const DFTFunctionalModel     *self                ) ;
extern DFTFunctionalModel *DFTFunctionalModel_MakeFromIDs     ( const IntegerArray1D         *ids                 ,
                                                                const Boolean                 isSpinRestricted    ,
                                                                      Status                 *status              ) ;
# endif
