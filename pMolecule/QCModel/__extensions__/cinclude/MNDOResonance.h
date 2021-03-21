# ifndef _MNDORESONANCE
# define _MNDORESONANCE

# include "Coordinates3.h"
# include "GaussianBasisContainer.h"
# include "IntegerArray1D.h"
# include "MNDOParametersContainer.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void MNDO_ResonanceGradients ( const MNDOParametersContainer *parameters        ,
                                      const GaussianBasisContainer  *bases             ,
                                      const IntegerArray1D          *basisIndices      ,
                                      const Coordinates3            *coordinates3      ,
                                      const SymmetricMatrix         *dTotal            ,
                                            Coordinates3            *gradients3        ) ;
extern void MNDO_ResonanceIntegrals ( const MNDOParametersContainer *parameters        ,
                                      const GaussianBasisContainer  *bases             ,
                                      const IntegerArray1D          *basisIndices      ,
                                      const Coordinates3            *coordinates3      ,
                                            SymmetricMatrix         *oneElectronMatrix ) ;
# endif
