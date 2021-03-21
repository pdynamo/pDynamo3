# ifndef _MNDOELECTRONNUCLEARTEIS
# define _MNDOELECTRONNUCLEARTEIS

# include "BlockStorage.h"
# include "Coordinates3.h"
# include "DoubleSymmetricMatrix.h"
# include "Integer.h"
# include "IntegerArray1D.h"
# include "MNDOParametersContainer.h"
# include "RealArray2D.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void MNDO_ElectronNuclearTEIGradients   ( const MNDOParametersContainer *parameters           ,
                                                 const IntegerArray1D          *basisIndices         ,
                                                 const Coordinates3            *coordinates3         ,
                                                 const SymmetricMatrix         *dTotal               ,
                                                 const SymmetricMatrix         *dSpin                ,
                                                       Coordinates3            *gradients3           ) ;
extern void MNDO_ElectronNuclearTEIGradientsCI ( const Integer                  nActive              ,
                                                 const Integer                  nCore                ,
                                                 const Integer                  nOrbitals            ,
                                                 const MNDOParametersContainer *parameters           ,
                                                 const IntegerArray1D          *basisIndices         ,
                                                 const Coordinates3            *coordinates3         ,
                                                 const DoubleSymmetricMatrix   *twoPDM               ,
                                                 const RealArray2D             *orbitals             ,
                                                 const SymmetricMatrix         *dCore                ,
                                                 const SymmetricMatrix         *dHF                  ,
                                                 const SymmetricMatrix         *dTotalZ              ,
                                                 const SymmetricMatrix         *onePDM               ,
                                                 const SymmetricMatrix         *zMatrix              ,
                                                       Coordinates3            *gradients3           ) ;
extern void MNDO_ElectronNuclearTEIIntegrals   ( const MNDOParametersContainer *parameters           ,
                                                 const IntegerArray1D          *basisIndices         ,
                                                 const Coordinates3            *coordinates3         ,
                                                       SymmetricMatrix         *oneElectronMatrix    ,
                                                       BlockStorage           **twoElectronIntegrals ) ;
# endif
