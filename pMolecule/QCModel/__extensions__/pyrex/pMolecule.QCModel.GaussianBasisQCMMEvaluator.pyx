"""Gaussian basis QC/MM density evaluator."""

from .QCDefinitions import BasisRepresentation
from .QCModelError  import QCModelError

# . Everything in atomic units.

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class GaussianBasisQCMMEvaluator:

    def Core ( self                                ,
               RealArray1D  chargesI      not None ,
               RealArray1D  chargesJ      not None ,
               Coordinates3 coordinates3I not None ,
               Coordinates3 coordinates3J not None ,
               Selection    selectionJ             ,
               Coordinates3 gradients3I            ,
               Coordinates3 gradients3J            ):
        """The core-MM energy."""
        cdef CReal         energy        = 0.0
        cdef CRealArray2D *cGradients3I = NULL
        cdef CRealArray2D *cGradients3J = NULL
        cdef CSelection   *cSelectionJ  = NULL
        if ( gradients3I is not None ) and \
           ( gradients3J is not None ):
            cGradients3I = gradients3I.cObject
            cGradients3J = gradients3J.cObject
        if selectionJ is not None: cSelectionJ = selectionJ.cObject
        energy = GaussianBasisContainer_NuclearNuclearEnergy ( chargesI.cObject      ,
                                                               chargesJ.cObject      ,
                                                               coordinates3I.cObject ,
                                                               coordinates3J.cObject ,
                                                               NULL                  ,
                                                               cSelectionJ           ,
                                                               NULL                  ,
                                                               NULL                  ,
                                                               NULL                  ,
                                                               NULL                  ,
                                                               cGradients3I          ,
                                                               cGradients3J          )
        return energy

    def Gradients ( self                                          ,
                    GaussianBasisContainer bases         not None ,
                    RealArray1D            charges       not None ,
                    RealArray2D            coordinates3  not None ,
                    RealArray2D            coordinates3G not None ,
                    Selection              selectionG             ,
                    SymmetricMatrix        density       not None ,
                    RealArray2D            gradients3    not None ,
                    RealArray2D            gradients3G   not None ):
        """Gradients of electron-MM energy."""
        cdef IntegerArray1D  basisIndices
        cdef CSelection     *cSelectionG = NULL
        cdef CStatus         cStatus     = CStatus_OK
        basisIndices = bases._centerFunctionPointers[self.basisRepresentation]
        if selectionG is not None: cSelectionG = selectionG.cObject
        GaussianBasisContainerIntegrals_ElectronNuclearD ( bases.cObject         ,
                                                           basisIndices.cObject  ,
                                                           charges.cObject       ,
                                                           NULL                  ,
                                                           NULL                  ,
                                                           coordinates3.cObject  ,
                                                           coordinates3G.cObject ,
                                                           cSelectionG           ,
                                                           density.cObject       ,
                                                           gradients3.cObject    ,
                                                           gradients3G.cObject   ,
                                                           &cStatus              )
        if cStatus != CStatus_OK: raise QCModelError ( "Error evaluating QC/MM gradients." )

    def Integrals ( self                                              ,
                    GaussianBasisContainer bases             not None ,
                    RealArray1D            charges           not None ,
                    RealArray2D            coordinates3      not None ,
                    RealArray2D            coordinates3G     not None ,
                    Selection              selectionG                 ,
                    SymmetricMatrix        oneElectronMatrix not None ):
        """Integrals of electron-MM interaction."""
        cdef IntegerArray1D  basisIndices
        cdef CSelection     *cSelectionG = NULL
        cdef CStatus         cStatus     = CStatus_OK
        basisIndices = bases._centerFunctionPointers[self.basisRepresentation]
        if selectionG is not None: cSelectionG = selectionG.cObject
        GaussianBasisContainerIntegrals_ElectronNuclear  ( bases.cObject             ,
                                                           basisIndices.cObject      ,
                                                           charges.cObject           ,
                                                           NULL                      ,
                                                           NULL                      ,
                                                           coordinates3.cObject      ,
                                                           coordinates3G.cObject     ,
                                                           cSelectionG               ,
                                                           oneElectronMatrix.cObject ,
                                                           &cStatus                  )
        if cStatus != CStatus_OK: raise QCModelError ( "Error evaluating QC/MM integrals." )

    @property
    def basisRepresentation ( self ): return BasisRepresentation.Work
