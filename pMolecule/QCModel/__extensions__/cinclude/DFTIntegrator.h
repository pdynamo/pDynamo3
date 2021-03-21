# ifndef _DFTINTEGRATOR
# define _DFTINTEGRATOR

# include "Boolean.h"
# include "Coordinates3.h"
# include "DFTFunctionalModel.h"
# include "DFTGrid.h"
# include "GaussianBasisContainer.h"
# include "Real.h"
# include "Status.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void DFTIntegrator_Integrate ( const DFTFunctionalModel     *functionalModel    ,
                                            DFTGrid                *grid               ,
                                      const GaussianBasisContainer *gaussianBases      ,
                                            Coordinates3           *qcCoordinates      ,
                                      const SymmetricMatrix        *densityP           ,
                                      const SymmetricMatrix        *densityQ           ,
                                      const Boolean                 inCore             ,
                                      const Boolean                 isSpinUnrestricted ,
                                            Real                   *eQuad              ,
                                            Real                   *rhoQuad            ,
                                            SymmetricMatrix        *fockA              ,
                                            SymmetricMatrix        *fockB              ,
                                            Coordinates3           *gradients3         ,
                                            Status                 *status             ) ;
# endif
