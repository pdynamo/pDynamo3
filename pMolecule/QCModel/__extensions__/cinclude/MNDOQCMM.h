# ifndef _MNDOQCMM
# define _MNDOQCMM

# include "BlockStorage.h"
# include "Coordinates3.h"
# include "CubicSplineContainer.h"
# include "IntegerArray1D.h"
# include "MNDOParametersContainer.h"
# include "PairList.h"
# include "Real.h"
# include "RealArray1D.h"
# include "Status.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void MNDO_QCMMGradients ( const IntegerArray1D          *atomIndices         ,
                                 const SymmetricMatrix         *dTotal              ,
                                       BlockStorage            *integrals           ,
                                       Coordinates3            *qcGradients3        ,
                                       Coordinates3            *mmGradients3        ,
                                       Status                  *status              ) ;
extern Real MNDO_QCMMIntegrals ( const MNDOParametersContainer *parameters          ,
                                 const IntegerArray1D          *basisIndices        ,
                                 const CubicSplineContainer    *splines             ,
                                 const Real                     cutOff              ,
                                 const Real                     eScale              ,
                                 const Coordinates3            *qcCoordinates3      ,
                                 const Coordinates3            *mmCoordinates3      ,
                                 const RealArray1D             *mmCharges           ,
                                       PairList                *pairList            ,
                                       SymmetricMatrix         *oneElectronMatrix   ,
                                       Coordinates3            *qcGradients3        ,
                                       Coordinates3            *mmGradients3        ,
                                       BlockStorage           **derivativeIntegrals ,
                                       Status                  *status              ) ;
# endif
