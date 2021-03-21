"""Defines a QC/MM Lennard-Jones model that employs cut-offs."""

from   pCore                         import CrossPairList                   , \
                                            logFile                         , \
                                            LogFileActive
from  .ImagePairListContainer        import ImagePairListContainer
from  .NBDefaults                    import _CheckCutOffs                   , \
                                            _DefaultGeneratorCutOff         , \
                                            _DefaultPairwiseInteractionABFS , \
                                            _ImageScanData                  , \
                                            _MMGrid                         , \
                                            _MMOccupancy                    , \
                                            _NonUpdatablePairLists          , \
                                            _PairListStatistics             , \
                                            _UpdatablePairLists
from  .PairwiseInteractionABFS       import PairwiseInteractionABFS
from  .PairwiseInteractionSplineABFS import PairwiseInteractionSplineABFS
from  .QCMMLennardJonesModel         import QCMMLennardJonesModel
from ..EnergyModel                   import EnergyClosurePriority

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCMMLennardJonesModelCutOff ( QCMMLennardJonesModel ):
    """Define a cut-off QC/MM Lennard-Jones model."""

    _attributable             = dict ( QCMMLennardJonesModel._attributable )
    _classLabel               = "CutOff QC/MM Lennard-Jones Model"
    _pairwiseInteractionClass = ( PairwiseInteractionABFS, PairwiseInteractionSplineABFS )
    _summarizable             = dict ( QCMMLennardJonesModel._summarizable )
    _attributable.update ( { "generator" : None } )
    _summarizable.update ( { "generator" : None } )

    def _CheckOptions ( self ):
        """Check options."""
        if self.generator is None:
            self.generator = _DefaultGeneratorCutOff ( )
        elif not isinstance ( self.generator, PairListGenerator ):
            raise TypeError ( "Invalid pairlist generator attribute." )
        if self.pairwiseInteraction is None:
            self.pairwiseInteraction = _DefaultPairwiseInteractionABFS ( )
        elif not isinstance ( pairwiseInteraction, self.__class__._pairwiseInteractionClass ):
            raise TypeError ( "Invalid pairwise interaction attribute." )              
        _CheckCutOffs ( self )
        return self

    def Energy ( self, target ):
        """Energy 1-5+."""
        energies     = {}
        scratch      = target.scratch
        coordinates3 = scratch.Get ( "coordinates3NB", target.coordinates3 )
        pNode        = scratch.GetSetNode ( _UpdatablePairLists )
        pairList     = pNode.Get ( "qcmmLJ", None )
        if pairList is None:
            pairList = self.generator.CrossPairListFromSingleCoordinates3 ( coordinates3                     ,
                                                                            None                             ,
                                                                            target.qcState.pureQCAtoms       ,
                                                                            target.mmState.mmAtoms           ,
                                                                            target.freeAtoms                 ,
                                                                            target.mmState.exclusions        ,
                                                                            False                            ,
                                                                            pNode.Get ( _MMGrid     , None ) ,
                                                                            pNode.Get ( _MMOccupancy, None ) )
            pNode.qcmmLJ = pairList
            sNode = scratch.Get ( _PairListStatistics )
            n     = float ( len ( pairList ) )
            sNode[ "QC/MM LJ Pairs" ]  = n
            sNode["<QC/MM LJ Pairs>"] += n
        if len ( pairList ) > 0:
            gradients3 = scratch.Get ( "gradients3", None )
            ( eElectrostatic, eLennardJones ) = self.pairwiseInteraction.MMMMEnergy ( None                         ,
                                                                                      None                         ,
                                                                                      target.mmState.ljTypeIndices ,
                                                                                      target.mmState.ljTypeIndices ,
                                                                                      target.mmState.ljParameters  ,
                                                                                      0.0                          ,
                                                                                      1.0                          ,
                                                                                      coordinates3                 ,
                                                                                      coordinates3                 ,
                                                                                      pairList                     ,
                                                                                      gradients3                   ,
                                                                                      gradients3                   )
            energies.update ( { "QC/MM Lennard-Jones" : eLennardJones } )
        return energies

    def Energy14 ( self, target ):
        """Energy 1-4."""
        energies = {}
        scratch  = target.scratch
        pNode    = scratch.GetSetNode ( _NonUpdatablePairLists )
        pairList = pNode.Get ( "qcmm14LJ", None )
        if pairList is None:
            pairList = CrossPairList.FromSelfPairList ( target.mmState.interactions14 ,
                                                        target.qcState.pureQCAtoms    ,
                                                        target.mmState.mmAtoms        ,
                                                        target.freeAtoms              )
            pNode.qcmm14LJ = pairList
            sNode          = scratch.Get ( _PairListStatistics )
            sNode["QC/MM 1-4 LJ Pairs"] = float ( len ( pairList ) )
        if len ( pairList ) > 0:
            coordinates3 = scratch.Get ( "coordinates3NB", target.coordinates3 )
            gradients3   = scratch.Get ( "gradients3"    , None                )
            ( eElectrostatic, eLennardJones ) = self.pairwiseInteraction.MMMMEnergy ( None                          ,
                                                                                      None                          ,
                                                                                      target.mmState.ljTypeIndices  ,
                                                                                      target.mmState.ljTypeIndices  ,
                                                                                      target.mmState.ljParameters14 ,
                                                                                      0.0                           ,
                                                                                      1.0                           ,
                                                                                      coordinates3                  ,
                                                                                      coordinates3                  ,
                                                                                      pairList                      ,
                                                                                      gradients3                    ,
                                                                                      gradients3                    )
            energies.update ( { "QC/MM 1-4 Lennard-Jones" : eLennardJones } )
        return energies

    def EnergyClosures ( self, target ):
        """Return energy closures."""
        def a ( ): target.scratch.energyTerms.update ( self.Energy      ( target ) )
        def b ( ): target.scratch.energyTerms.update ( self.Energy14    ( target ) )
        def c ( ): target.scratch.energyTerms.update ( self.EnergyImage ( target ) )
        return [ ( EnergyClosurePriority.IndependentEnergyTerm , a, "QC/MM Lennard-Jones"       ) ,
                 ( EnergyClosurePriority.IndependentEnergyTerm , b, "QC/MM 1-4 Lennard-Jones"   ) ,
                 ( EnergyClosurePriority.IndependentEnergyTerm , c, "QC/MM Image Lennard-Jones" ) ]

    def EnergyImage ( self, target ):
        """Image energy."""
        energies           = {}
        symmetryParameters = target.symmetryParameters
        if symmetryParameters is not None:
            scratch      = target.scratch
            coordinates3 = scratch.Get ( "coordinates3NB", target.coordinates3 )
            pNode        = scratch.GetSetNode ( _UpdatablePairLists )
            pairList     = pNode.Get ( "qcmmImageLJ", None )
            if pairList is None:
                pairList = ImagePairListContainer.Constructor ( self.generator                      ,
                                                                target.qcState.pureQCAtoms          ,
                                                                target.mmState.mmAtoms              ,
                                                                target.freeAtoms                    ,
                                                                coordinates3                        ,
                                                                coordinates3                        ,
                                                                symmetryParameters                  ,
                                                                target.symmetry.transformations     ,
                                                                pNode.Get ( _ImageScanData , None ) ,
                                                                pNode.Get ( _MMGrid        , None ) ,
                                                                pNode.Get ( _MMOccupancy   , None ) ,
                                                                False                               )
                pNode.qcmmImageLJ = pairList
                sNode = scratch.Get ( _PairListStatistics )
                nI    = float ( pairList.numberOfImages )
                nP    = float ( pairList.numberOfPairs  )
                sNode[ "QC/MM Image LJ Images" ]  = nI
                sNode[ "QC/MM Image LJ Pairs"  ]  = nP
                sNode["<QC/MM Image LJ Images>"] += nI
                sNode["<QC/MM Image LJ Pairs>" ] += nP
            if len ( pairList ) > 0:
                ( eElectrostatic, eLennardJones ) = self.pairwiseInteraction.MMMMEnergyImage ( None                         ,
                                                                                               target.mmState.ljTypeIndices ,
                                                                                               target.mmState.ljParameters  ,
                                                                                               0.0                          ,
                                                                                               coordinates3                 ,
                                                                                               symmetryParameters           ,
                                                                                               pairList                     ,
                                                                                               scratch.Get ( "gradients3"                , None ) ,
                                                                                               scratch.Get ( "symmetryParameterGradients", None ) )
                energies.update ( { "QC/MM Image Lennard-Jones" : eLennardJones } )
        return energies

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
