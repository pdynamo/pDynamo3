"""Defines a simple full QC/MM Lennard-Jones model."""

from   pCore                     import CrossPairList
from  .NBDefaults                import _NonUpdatablePairLists  , \
                                        _PairListStatistics
from  .PairwiseInteractionFull   import PairwiseInteractionFull
from  .QCMMLennardJonesModel     import QCMMLennardJonesModel
from ..EnergyModel               import EnergyClosurePriority

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCMMLennardJonesModelFull ( QCMMLennardJonesModel ):
    """Define a full QC/MM Lennard-Jones model."""

    # . Defaults.
    _classLabel               = "Full QC/MM Lennard-Jones Model"
    _pairwiseInteractionClass = PairwiseInteractionFull

    def EnergyClosures ( self, target ):
        """Return energy closures."""
        def a ( ):
            target.scratch.energyTerms.update ( self.Energy ( target ) )
        def b ( ):
            target.scratch.energyTerms.update ( self.Energy14 ( target ) )
        return [ ( EnergyClosurePriority.IndependentEnergyTerm , a, "QC/MM LJ Evaluation"     ) ,
                 ( EnergyClosurePriority.IndependentEnergyTerm , b, "QC/MM 1-4 LJ Evaluation" ) ]

    def Energy ( self, target ):
        """Energy 1-5+."""
        energies = {}
        scratch  = target.scratch
        pNode    = scratch.GetSetNode ( _NonUpdatablePairLists )
        pairList = pNode.Get ( "qcmmLennardJones", None )
        if pairList is None:
            pairList = CrossPairList.FromSelfPairList ( target.mmState.exclusions  ,
                                                        target.qcState.pureQCAtoms ,
                                                        target.mmState.mmAtoms     ,
                                                        target.freeAtoms           , excluded = True )
            pNode.qcmmLennardJones             = pairList
            sNode                              = scratch.Get ( _PairListStatistics )
            sNode["QC/MM Lennard-Jones Pairs"] = float ( len ( pairList ) )
        if len ( pairList ) > 0:
            gradients3 = scratch.Get ( "gradients3", None )
            ( eElectrostatic, eLennardJones ) = self.pairwiseInteraction.MMMMEnergy ( None                         ,
                                                                                      None                         ,
                                                                                      target.mmState.ljTypeIndices ,
                                                                                      target.mmState.ljTypeIndices ,
                                                                                      target.mmState.ljParameters  ,
                                                                                      0.0                          ,
                                                                                      1.0                          ,
                                                                                      target.coordinates3          ,
                                                                                      target.coordinates3          ,
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
        pairList = pNode.Get ( "qcmm14LennardJones", None )
        if pairList is None:
            pairList = CrossPairList.FromSelfPairList ( target.mmState.interactions14 ,
                                                        target.qcState.pureQCAtoms    ,
                                                        target.mmState.mmAtoms        ,
                                                        target.freeAtoms              )
            pNode.qcmm14LennardJones               = pairList
            sNode                                  = scratch.Get ( _PairListStatistics )
            sNode["QC/MM 1-4 Lennard-Jones Pairs"] = float ( len ( pairList ) )
        if len ( pairList ) > 0:
            gradients3 = scratch.Get ( "gradients3", None )
            ( eElectrostatic, eLennardJones ) = self.pairwiseInteraction.MMMMEnergy ( None                          ,
                                                                                      None                          ,
                                                                                      target.mmState.ljTypeIndices  ,
                                                                                      target.mmState.ljTypeIndices  ,
                                                                                      target.mmState.ljParameters14 ,
                                                                                      0.0                           ,
                                                                                      1.0                           ,
                                                                                      target.coordinates3           ,
                                                                                      target.coordinates3           ,
                                                                                      pairList                      ,
                                                                                      gradients3                    ,
                                                                                      gradients3                    )
            energies.update ( { "QC/MM 1-4 Lennard-Jones" : eLennardJones } )
        return energies

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
