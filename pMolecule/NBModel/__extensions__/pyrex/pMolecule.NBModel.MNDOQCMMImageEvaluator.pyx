"""MNDO QC/MM density evaluator."""

from .NBModelError import NBModelError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class MNDOQCMMImageEvaluator ( MNDOQCMMEvaluator ):

    def GradientsImage ( self                                                           ,
                         IntegerArray1D             atomIndices                not None ,
                         SymmetricMatrix            dTotal                     not None ,
                         RealArray2D                qcCoordinates3             not None ,
                         RealArray2D                mmCoordinates3             not None ,
                         SymmetryParameters         symmetryParameters         not None ,
                         ImagePairListContainer     imagePairLists             not None ,
                         BlockStorageContainer      integralContainer          not None ,
                         RealArray2D                qcGradients3               not None ,
                         RealArray2D                mmGradients3               not None ,
                         SymmetryParameterGradients symmetryParameterGradients not None ):
        """Image gradients."""
        cdef CStatus cStatus = CStatus_OK
        MNDOQCMMImage_QCMMGradientsImage  ( atomIndices.cObject                ,
                                            dTotal.cObject                     ,
                                            qcCoordinates3.cObject             ,
                                            mmCoordinates3.cObject             ,
                                            symmetryParameters.cObject         ,
                                            imagePairLists.cObject             ,
                                            integralContainer.cObject          ,
                                            qcGradients3.cObject               ,
                                            mmGradients3.cObject               ,
                                            symmetryParameterGradients.cObject ,
                                            &cStatus                           )
        if cStatus != CStatus_OK: raise NBModelError ( "Error evaluating QC/MM image gradients." )

    def IntegralsImage ( self                                                   ,
                         MNDOParametersContainer    parameters         not None ,
                         IntegerArray1D             basisIndices       not None ,
                         CubicSplineContainer       splines            not None ,
                         CReal                      cutOff                      ,
                         CReal                      electrostaticScale          ,
                         RealArray2D                qcCoordinates3     not None ,
                         RealArray2D                mmCoordinates3     not None ,
                         SymmetryParameters         symmetryParameters not None ,
                         RealArray1D                mmCharges          not None ,
                         ImagePairListContainer     imagePairLists     not None ,
                         SymmetricMatrix            oneElectronMatrix  not None ,
                         Coordinates3               qcGradients3                ,
                         Coordinates3               mmGradients3                ,
                         SymmetryParameterGradients symmetryParameterGradients  ):
        """Image integrals."""
        cdef BlockStorageContainer        integralContainer
        cdef CBlockStorageContainer      *cIntegralContainer          = NULL 
        cdef CReal                        energy                             
        cdef CRealArray2D                *cMMGradients3               = NULL 
        cdef CRealArray2D                *cQCGradients3               = NULL 
        cdef CSymmetryParameterGradients *cSymmetryParameterGradients = NULL
        cdef CStatus                      cStatus                     = CStatus_OK
        if ( mmGradients3               is not None ) and \
           ( qcGradients3               is not None ) and \
           ( symmetryParameterGradients is not None ):
            integralContainer           = BlockStorageContainer ( imagePairLists.numberOfImages )
            cIntegralContainer          = integralContainer.cObject
            cMMGradients3               = mmGradients3.cObject
            cQCGradients3               = qcGradients3.cObject
            cSymmetryParameterGradients = symmetryParameterGradients.cObject
        else:
            integralContainer = None
        energy = MNDOQCMMImage_QCMMPotentialsImage ( parameters.cObject          ,
                                                     basisIndices.cObject        ,
                                                     splines.cObject             ,
                                                     cutOff                      ,
                                                     electrostaticScale          ,
                                                     qcCoordinates3.cObject      ,
                                                     mmCoordinates3.cObject      ,
                                                     symmetryParameters.cObject  ,
                                                     mmCharges.cObject           ,
                                                     imagePairLists.cObject      ,
                                                     oneElectronMatrix.cObject   ,
                                                     cQCGradients3               ,
                                                     cMMGradients3               ,
                                                     cSymmetryParameterGradients ,
                                                     cIntegralContainer          ,
                                                     &cStatus                    )
        if cStatus != CStatus_OK: raise NBModelError ( "Error evaluating QC/MM image integrals." )
        return ( energy, integralContainer )
