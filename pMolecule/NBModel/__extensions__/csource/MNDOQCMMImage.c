/*==================================================================================================================================
! . MNDO QC/MM interactions for images.
!=================================================================================================================================*/

# include "Boolean.h"
# include "Integer.h"
# include "MNDOQCMM.h"
# include "MNDOQCMMImage.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Image QC/MM gradients.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDOQCMMImage_QCMMGradientsImage ( const IntegerArray1D             *atomIndices                ,
                                        const SymmetricMatrix            *dTotal                     ,
                                              Coordinates3               *coordinates3A              ,
                                              Coordinates3               *coordinates3B              ,
                                              SymmetryParameters         *symmetryParameters         ,
                                              ImagePairListContainer     *imagePairLists             ,
                                              BlockStorageContainer      *integralContainer          ,
                                              Coordinates3               *gradients3A                ,
                                              Coordinates3               *gradients3B                ,
                                              SymmetryParameterGradients *symmetryParameterGradients ,
                                              Status                     *status                     )
{
    if ( ( atomIndices                != NULL ) &&
         ( dTotal                     != NULL ) &&
         ( coordinates3A              != NULL ) &&
         ( coordinates3B              != NULL ) &&
         ( symmetryParameters         != NULL ) &&
         ( imagePairLists             != NULL ) &&
         ( integralContainer          != NULL ) &&
         ( gradients3A                != NULL ) &&
         ( gradients3B                != NULL ) &&
         ( symmetryParameterGradients != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto ImagePairListIterator iterator ;
        ImagePairListIterator_Initialize ( &iterator                  ,
                                           imagePairLists             ,
                                           coordinates3B              ,
                                           symmetryParameters         ,
                                           gradients3B                ,
                                           symmetryParameterGradients ,
                                           status                     ) ;
        if ( Status_IsOK ( status ) )
        {
            auto BlockStorage *integrals ;
            auto Integer       image = 0 ;
            while ( ImagePairListIterator_Next ( &iterator ) )
            {
                integrals = integralContainer->entries[image] ;
                MNDO_QCMMGradients ( atomIndices          ,
                                     dTotal               ,
                                     integrals            ,
                                     gradients3A          ,
                                     iterator.iGradients3 ,
                                     status               ) ;
                ImagePairListIterator_Gradients ( &iterator ) ;
                image ++ ;
            }
        }
        ImagePairListIterator_Finalize ( &iterator ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Image QC/MM potentials.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real MNDOQCMMImage_QCMMPotentialsImage ( const MNDOParametersContainer    *parameters                 ,
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
                                               Status                     *status                     )
{
    Real eCore = 0.0e+00 ;
    if ( ( parameters         != NULL    ) &&
         ( basisIndices       != NULL    ) &&
         ( splines            != NULL    ) &&
         ( eScale             != 0.0e+00 ) &&
         ( qcCoordinates3     != NULL    ) &&
         ( mmCoordinates3     != NULL    ) &&
         ( symmetryParameters != NULL    ) &&
         ( mmCharges          != NULL    ) &&
         ( imagePairLists     != NULL    ) &&
         ( oneElectronMatrix  != NULL    ) &&
         Status_IsOK ( status ) )
    {
        auto ImagePairListIterator iterator ;
        ImagePairListIterator_Initialize ( &iterator                  ,
                                           imagePairLists             ,
                                           mmCoordinates3             ,
                                           symmetryParameters         ,
                                           mmGradients3               ,
                                           symmetryParameterGradients ,
                                           status                     ) ;
        if ( Status_IsOK ( status ) )
        {
            auto BlockStorage *integrals = NULL, **pIntegrals = NULL ;
            auto Boolean       doGradients ;
            auto Integer       image     = 0 ;
            doGradients = ( integralContainer          != NULL ) &&
                          ( mmGradients3               != NULL ) &&
                          ( qcGradients3               != NULL ) &&
                          ( symmetryParameterGradients != NULL ) ;
            if ( doGradients ) pIntegrals = &integrals ;
            while ( ImagePairListIterator_Next ( &iterator ) )
            {
                eCore += MNDO_QCMMIntegrals ( parameters              ,
                                              basisIndices            ,
                                              splines                 ,
                                              cutOff                  ,
                                              eScale * iterator.scale ,
                                              qcCoordinates3          ,
                                              iterator.iCoordinates3  ,
                                              mmCharges               ,
                                              iterator.pairList       ,
                                              oneElectronMatrix       ,
                                              qcGradients3            ,
                                              iterator.iGradients3    ,
                                              pIntegrals              ,
                                              status                  ) ;
                if ( doGradients )
                {
                    integralContainer->entries[image] = integrals ;
                    integrals = NULL ;
                    ImagePairListIterator_Gradients ( &iterator ) ;
                }
                image ++ ;
            }
        }
        ImagePairListIterator_Finalize ( &iterator ) ;
    }
    return eCore ;
}
