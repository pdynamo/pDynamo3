"""Gaussian basis QC/MM density evaluator."""

from .GaussianBasisError import GaussianBasisError

# . Everything in atomic units.

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class GaussianBasisQCMMEvaluator:

    def f1Cm1R1 ( self                                            ,
                  GaussianBasisContainer bases           not None ,
                  RealArray1D            charges         not None ,
                  RealArray2D            coordinates3    not None ,
                  RealArray2D            coordinates3G   not None ,
                  Selection              selectionG               ,
                  RealArray1D            fitCoefficients not None ,
                  RealArray2D            gradients3      not None ,
                  RealArray2D            gradients3G     not None ):
        """Gradients of fit-MM energy."""
        cdef CSelection *cSelectionG = NULL
        cdef CStatus     cStatus     = CStatus_OK
        if selectionG is not None: cSelectionG = selectionG.cObject
        GaussianBasisContainerIntegrals_f1Cm1R1 ( bases.cObject           , 
                                                  charges.cObject         , 
                                                  NULL                    , 
                                                  NULL                    , 
                                                  coordinates3.cObject    , 
                                                  coordinates3G.cObject   , 
                                                  cSelectionG             , 
                                                  fitCoefficients.cObject ,
                                                  gradients3.cObject      ,
                                                  gradients3G.cObject     ,
                                                  &cStatus                )
        if cStatus != CStatus_OK: raise GaussianBasisError ( "Error evaluating QC/MM fit gradients." )

    def f1Cm1V ( self                                          ,
                 GaussianBasisContainer bases         not None ,
                 RealArray1D            charges       not None ,
                 RealArray2D            coordinates3  not None ,
                 RealArray2D            coordinates3G not None ,
                 Selection              selectionG             ,
                 RealArray1D            integrals     not None ):
        """Integrals of fit-MM interaction."""
        cdef CSelection *cSelectionG = NULL
        cdef CStatus     cStatus     = CStatus_OK
        if selectionG is not None: cSelectionG = selectionG.cObject
        GaussianBasisContainerIntegrals_f1Cm1V ( bases.cObject             ,
                                                 charges.cObject           ,
                                                 NULL                      ,
                                                 NULL                      ,
                                                 coordinates3.cObject      ,
                                                 coordinates3G.cObject     ,
                                                 cSelectionG               ,
                                                 integrals.cObject         ,
                                                 &cStatus                  )
        if cStatus != CStatus_OK: raise GaussianBasisError ( "Error evaluating QC/MM fit integrals." )

    def f2Cm1R1 ( self                                          ,
                  GaussianBasisContainer bases         not None ,
                  RealArray1D            charges       not None ,
                  RealArray2D            coordinates3  not None ,
                  RealArray2D            coordinates3G not None ,
                  Selection              selectionG             ,
                  SymmetricMatrix        density       not None ,
                  RealArray2D            gradients3    not None ,
                  RealArray2D            gradients3G   not None ):
        """Gradients of electron-MM energy."""
        cdef CSelection *cSelectionG = NULL
        cdef CStatus     cStatus     = CStatus_OK
        if selectionG is not None: cSelectionG = selectionG.cObject
        GaussianBasisContainerIntegrals_f2Cm1R1 ( bases.cObject         ,
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
        if cStatus != CStatus_OK: raise GaussianBasisError ( "Error evaluating QC/MM gradients." )

    def f2Cm1V ( self                                              ,
                 GaussianBasisContainer bases             not None ,
                 RealArray1D            charges           not None ,
                 RealArray2D            coordinates3      not None ,
                 RealArray2D            coordinates3G     not None ,
                 Selection              selectionG                 ,
                 SymmetricMatrix        oneElectronMatrix not None ):
        """Integrals of electron-MM interaction."""
        cdef CSelection *cSelectionG = NULL
        cdef CStatus     cStatus     = CStatus_OK
        if selectionG is not None: cSelectionG = selectionG.cObject
        GaussianBasisContainerIntegrals_f2Cm1V ( bases.cObject             ,
                                                 charges.cObject           ,
                                                 NULL                      ,
                                                 NULL                      ,
                                                 coordinates3.cObject      ,
                                                 coordinates3G.cObject     ,
                                                 cSelectionG               ,
                                                 oneElectronMatrix.cObject ,
                                                 &cStatus                  )
        if cStatus != CStatus_OK: raise GaussianBasisError ( "Error evaluating QC/MM integrals." )

    def m1Cn1ER1 ( self                                ,
                   RealArray1D  chargesI      not None ,
                   RealArray1D  chargesJ      not None ,
                   Coordinates3 coordinates3I not None ,
                   Coordinates3 coordinates3J not None ,
                   Selection    selectionJ             ,
                   Coordinates3 gradients3I            ,
                   Coordinates3 gradients3J            ):
        """The QC core-MM point charge energy and gradients."""
        cdef CReal         energy        = 0.0
        cdef CRealArray2D *cGradients3I = NULL
        cdef CRealArray2D *cGradients3J = NULL
        cdef CSelection   *cSelectionJ  = NULL
        if ( gradients3I is not None ) and \
           ( gradients3J is not None ):
            cGradients3I = gradients3I.cObject
            cGradients3J = gradients3J.cObject
        if selectionJ is not None: cSelectionJ = selectionJ.cObject
        energy = GaussianBasisContainer_m1Cn1ER1 ( chargesI.cObject      ,
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
