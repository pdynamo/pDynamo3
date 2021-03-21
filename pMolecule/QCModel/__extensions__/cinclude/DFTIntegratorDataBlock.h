# ifndef _DFTINTEGRATORDATABLOCK
# define _DFTINTEGRATORDATABLOCK

# include "Boolean.h"
# include "Integer.h"
# include "RealArray1D.h"
# include "RealArray2D.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The DFT integration block data view type. */
typedef struct {
    RealArray1D dRhoX         ;
    RealArray1D dRhoY         ;
    RealArray1D dRhoZ         ;
    RealArray1D laplacianRho  ;
    RealArray1D rho           ; 
    RealArray1D sigma         ;
    RealArray1D tau           ;
    RealArray1D vLaplacianRho ;
    RealArray1D vRho          ;
    RealArray1D vSigma        ;
    RealArray1D vTau          ;
} DFTIntegratorDataBlockView ;

/* . The DFT integration block data type. */
typedef struct {
    Boolean      hasLocalData        ;
    Integer      numberOfPoints      ;
    RealArray1D *eXC                 ;
    RealArray1D *localEXC            ;
    RealArray2D *dRhoX               ;
    RealArray2D *dRhoY               ;
    RealArray2D *dRhoZ               ;
    RealArray2D *localVLaplacianRho  ;
    RealArray2D *localVRho           ;
    RealArray2D *localVSigma         ;
    RealArray2D *localVTau           ;
    RealArray2D *laplacianRho        ;
    RealArray2D *rho                 ;
    RealArray2D *sigma               ;
    RealArray2D *tau                 ;
    RealArray2D *vLaplacianRho       ;
    RealArray2D *vRho                ;
    RealArray2D *vSigma              ;
    RealArray2D *vTau                ;
    DFTIntegratorDataBlockView viewP ;
    DFTIntegratorDataBlockView viewQ ;
    RealArray1D  sigmaPQ             ;
    RealArray1D  vSigmaPQ            ;
} DFTIntegratorDataBlock ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void                    DFTIntegratorDataBlock_Accumulate     (       DFTIntegratorDataBlock     *self                ) ;
extern DFTIntegratorDataBlock *DFTIntegratorDataBlock_Allocate       ( const Integer                     numberOfFunctionals ,
                                                                       const Integer                     numberOfPoints      ,
                                                                       const Boolean                     hasSigma            ,
                                                                       const Boolean                     hasLaplacian        ,
                                                                       const Boolean                     hasTau              ,
                                                                       const Boolean                     isSpinRestricted    ,
                                                                             Status                     *status              ) ;
extern void                    DFTIntegratorDataBlock_Deallocate     (       DFTIntegratorDataBlock    **self                ) ;
extern void                    DFTIntegratorDataBlock_Initialize     (       DFTIntegratorDataBlock     *self                ) ;
extern void                    DFTIntegratorDataBlock_InitializeView (       DFTIntegratorDataBlock     *self                ,
                                                                       const Integer                     c                   ,
                                                                             DFTIntegratorDataBlockView *view                ) ;
# endif
