"""Defines a cut-off QC/QC electrostatic multipole model."""

from   pCore                         import logFile                                   , \
                                            LogFileActive
from   pScientific.Arrays            import Array                                     , \
                                            StorageType
from  .ImagePairListContainer        import ImagePairListContainer
from  .NBDefaults                    import _CheckCutOffs                             , \
                                            _DefaultGeneratorCutOff                   , \
                                            _DefaultPairwiseInteractionSplineABFSQCQC , \
                                            _ImageScanData                            , \
                                            _QCGrid                                   , \
                                            _QCOccupancy                              , \
                                            _PairListStatistics                       , \
                                            _UpdatablePairLists
from  .PairwiseInteractionABFS       import PairwiseInteractionABFS
from  .PairwiseInteractionSplineABFS import PairwiseInteractionSplineABFS
from  .QCMMElectrostaticModel        import QCMMElectrostaticModel
from ..EnergyModel                   import EnergyClosurePriority                     , \
                                            EnergyModel

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCQCElectrostaticModelMultipoleCutOff ( EnergyModel ):
    """Define a cut-off QC/QC electrostatic model."""

    _attributable             = dict ( EnergyModel._attributable )
    _classLabel               = "CutOff Multipole QC/QC Electrostatic Model"
    _pairwiseInteractionClass = ( PairwiseInteractionABFS, PairwiseInteractionSplineABFS )
    _summarizable             = dict ( EnergyModel._summarizable )
    _attributable.update ( { "dielectric"          : 1.0  ,
                             "generator"           : None ,
                             "pairwiseInteraction" : None } )
    _summarizable.update ( { "dielectric"          : ( "Dielectric", "{:.3f}" ) ,
                             "generator"           : None                       ,
                             "pairwiseInteraction" : None                     } )


    def _CheckOptions ( self ):
        """Check options."""
        if self.generator is None:
            self.generator = _DefaultGeneratorCutOff ( )
        elif not isinstance ( self.generator, PairListGenerator ):
            raise TypeError ( "Invalid pairlist generator attribute." )
        if self.pairwiseInteraction is None:
            self.pairwiseInteraction = _DefaultPairwiseInteractionSplineABFSQCQC ( )
        elif not isinstance ( pairwiseInteraction, self.__class__._pairwiseInteractionClass ):
            raise TypeError ( "Invalid pairwise interaction attribute." )              
        _CheckCutOffs ( self )
        return self

    def EnergyClosures ( self, target ):
        """Return energy closures."""
        def a ( ): self.QCQCPotentials ( target )
        def b ( ): self.QCQCGradients  ( target )
        return [ ( EnergyClosurePriority.QCIntegrals , a, "QC/QC Image Electrostatic Potentials" ) ,
                 ( EnergyClosurePriority.QCGradients , b, "QC/QC Image Electrostatic Gradients"  ) ]

    def QCQCGradients ( self, target ):
        """QC/QC gradients."""
        scratch  = target.scratch
        pNode    = scratch.Get ( _UpdatablePairLists )
        pairList = pNode.Get ( "qcqcImageElectrostatic", None )
        if scratch.doGradients and ( pairList is not None ) and ( len ( pairList ) > 0 ):
            self.pairwiseInteraction.QCQCGradientsImage ( scratch.qcmmMultipoles             ,
                                                          1.0 / self.dielectric              ,
                                                          scratch.qcCoordinates3QCMM         ,
                                                          target.symmetryParameters          ,
                                                          pairList                           ,
                                                          scratch.qcGradients3QCMM           ,
                                                          scratch.symmetryParameterGradients )

    def QCQCPotentials ( self, target ):
        """QC/QC potentials."""
        symmetryParameters = target.symmetryParameters
        if symmetryParameters is not None:
            scratch = target.scratch
            pNode   = scratch.Get ( _UpdatablePairLists )
            if symmetryParameters.IsMinimumImageConventionSatisfied ( self.pairwiseInteraction.range ):
                pNode.qcqcImageElectrostatic = None
            else:
                pairList       = pNode.Get ( "qcqcImageElectrostatic", None )
                qcCoordinates3 = scratch.qcCoordinates3QCMM
                if pairList is None:
                    pairList = ImagePairListContainer.Constructor ( self.generator                      ,
                                                                    None                                ,
                                                                    None                                ,
                                                                    None                                ,
                                                                    qcCoordinates3                      ,
                                                                    qcCoordinates3                      ,
                                                                    symmetryParameters                  ,
                                                                    target.symmetry.transformations     ,
                                                                    pNode.Get ( _ImageScanData , None ) ,
                                                                    pNode.Get ( _QCGrid        , None ) ,
                                                                    pNode.Get ( _QCOccupancy   , None ) ,
                                                                    True                                )
                    pNode.qcqcImageElectrostatic = pairList
                    sNode = scratch.Get ( _PairListStatistics )
                    nI    = float ( pairList.numberOfImages )
                    nP    = float ( pairList.numberOfPairs  )
                    sNode[ "QC/QC Image Electrostatic Images" ]  = nI
                    sNode[ "QC/QC Image Electrostatic Pairs"  ]  = nP
                    sNode["<QC/QC Image Electrostatic Images>"] += nI
                    sNode["<QC/QC Image Electrostatic Pairs>" ] += nP
                if len ( pairList ) > 0:
                    potentials = scratch.Get ( "qcqcPotentials", None )
                    if potentials is None:
                        potentials = Array.WithExtent ( qcCoordinates3.rows, storageType = StorageType.Symmetric )
                        scratch.qcqcPotentials = potentials
                    potentials.Set ( 0.0 )
                    self.pairwiseInteraction.QCQCPotentialsImage ( 1.0 / self.dielectric ,
                                                                   qcCoordinates3        ,
                                                                   symmetryParameters    ,
                                                                   pairList              ,
                                                                   potentials            )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
