# ifndef _MNDOINTEGRALS
# define _MNDOINTEGRALS

# include "BlockStorage.h"
# include "Integer.h"
# include "MNDOParameters.h"
# include "RealArray1D.h"
# include "RealArray2D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void MNDOIntegrals_AddInOneCenterTEIs         ( const MNDOParameters *self                 ,
                                                       const Integer         i0                   ,
                                                             BlockStorage   *twoElectronIntegrals ) ;
extern void MNDOIntegrals_MolecularFrame2CIntegrals  ( const MNDOParameters *iData                ,
                                                       const Integer         i0                   ,
                                                       const Real           *xI                   ,
                                                       const MNDOParameters *jData                ,
                                                       const Integer         j0                   ,
                                                       const Real           *xJ                   ,
                                                             RealArray1D    *e1b                  ,
                                                             RealArray1D    *e2a                  ,
                                                             BlockStorage   *twoElectronIntegrals ) ;
extern void MNDOIntegrals_MolecularFrame2CIntegralsD ( const MNDOParameters *iData                ,
                                                       const Integer         i0                   ,
                                                       const Real           *xI                   ,
                                                       const MNDOParameters *jData                ,
                                                       const Integer         j0                   ,
                                                       const Real           *xJ                   ,
                                                       const RealArray1D    *dOneI                ,
                                                       const RealArray1D    *dOneJ                ,
                                                       const RealArray2D    *dTwoIJ               ,
                                                             Real           *gX                   ,
                                                             Real           *gY                   ,
                                                             Real           *gZ                   ) ;
# endif
