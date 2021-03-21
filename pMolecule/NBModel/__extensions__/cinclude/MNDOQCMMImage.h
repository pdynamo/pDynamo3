# ifndef _MNDOQCMMIMAGE
# define _MNDOQCMMIMAGE

# include "BlockStorageContainer.h"
# include "Coordinates3.h"
# include "CubicSplineContainer.h"
# include "ImagePairListContainer.h"
# include "IntegerArray1D.h"
# include "MNDOParametersContainer.h"
# include "Real.h"
# include "RealArray1D.h"
# include "Status.h"
# include "SymmetricMatrix.h"
# include "SymmetryParameterGradients.h"
# include "SymmetryParameters.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void MNDOQCMMImage_QCMMGradientsImage  ( const IntegerArray1D             *atomIndices                ,
                                                const SymmetricMatrix            *dTotal                     ,
                                                      Coordinates3               *coordinates3A              ,
                                                      Coordinates3               *coordinates3B              ,
                                                      SymmetryParameters         *symmetryParameters         ,
                                                      ImagePairListContainer     *imagePairLists             ,
                                                      BlockStorageContainer      *integralContainer          ,
                                                      Coordinates3               *gradients3A                ,
                                                      Coordinates3               *gradients3B                ,
                                                      SymmetryParameterGradients *symmetryParameterGradients ,
                                                      Status                     *status                     ) ;
extern Real MNDOQCMMImage_QCMMPotentialsImage ( const MNDOParametersContainer    *parameters                 ,
                                                const IntegerArray1D             *basisIndices               ,
                                                const CubicSplineContainer       *splines                    ,
                                                const Real                        cutOff                     ,
                                                const Real                        eScale                     ,
                                                      Coordinates3               *qcCoordinates3             ,
                                                      Coordinates3               *mmCoordinates3             ,
                                                      SymmetryParameters         *symmetryParameters         ,
                                                const RealArray1D                *mmCharges                  ,
                                                      ImagePairListContainer     *imagePairLists             ,
                                                      SymmetricMatrix            *oneElectronMatrix          ,
                                                      Coordinates3               *qcGradients3               ,
                                                      Coordinates3               *mmGradients3               ,
                                                      SymmetryParameterGradients *symmetryParameterGradients ,
                                                      BlockStorageContainer      *integralContainer          ,
                                                      Status                     *status                     ) ;
# endif
