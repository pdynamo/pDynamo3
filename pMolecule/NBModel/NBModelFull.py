"""Defines a simple full NB model."""

from   pCore                                          import SelfPairList
from   pMolecule.QCModel                              import QCModelDFT             , \
                                                             QCModelMNDO
from  .NBDefaults                                     import _NonUpdatablePairLists , \
                                                             _PairListStatistics
from  .NBModel                                        import NBModel
from  .PairwiseInteractionFull                        import PairwiseInteractionFull
from  .QCMMElectrostaticModelDensityFullGaussianBasis import QCMMElectrostaticModelDensityFullGaussianBasis
from  .QCMMElectrostaticModelDensityFullMNDO          import QCMMElectrostaticModelDensityFullMNDO
from  .QCMMElectrostaticModelMultipoleFull            import QCMMElectrostaticModelMultipoleFull
from  .QCMMLennardJonesModelFull                      import QCMMLennardJonesModelFull
from ..EnergyModel                                    import EnergyClosurePriority

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class NBModelFull ( NBModel ):
    """Define a full NB model."""

    # . Defaults.
    _classLabel               = "Full NB Model"
    _pairwiseInteractionClass = PairwiseInteractionFull

    def Energy ( self, target ):
        """Energy 1-5+."""
        energies = {}
        pNode    = target.scratch.GetSetNode ( _NonUpdatablePairLists )
        pairList = pNode.Get ( "mmmm", None )
        if pairList is None:
            pairList = SelfPairList.FromSelfPairList ( target.mmState.exclusions , 
                                                       len ( target.atoms )      ,
                                                       target.mmState.mmAtoms    ,
                                                       target.freeAtoms          , excluded = True )
            pNode.mmmm = pairList
            sNode      = target.scratch.Get ( _PairListStatistics )
            sNode["MM/MM Pairs"] = float ( len ( pairList ) )
        if len ( pairList ) > 0:
            gradients3 = target.scratch.Get ( "gradients3", None )
            ( eElectrostatic, eLennardJones ) = self.pairwiseInteraction.MMMMEnergy ( target.mmState.charges       ,
                                                                                      target.mmState.charges       ,
                                                                                      target.mmState.ljTypeIndices ,
                                                                                      target.mmState.ljTypeIndices ,
                                                                                      target.mmState.ljParameters  ,
                                                                                      ( 1.0 / self.dielectric )    ,
                                                                                        1.0                        ,
                                                                                      target.coordinates3          ,
                                                                                      target.coordinates3          ,
                                                                                      pairList                     ,
                                                                                      gradients3                   ,
                                                                                      gradients3                   )
            energies.update ( { "MM/MM Electrostatic" : eElectrostatic ,
                                "MM/MM Lennard-Jones" : eLennardJones  } )
        return energies

    def Energy14 ( self, target ):
        """Energy 1-4."""
        energies = {}
        pNode    = target.scratch.GetSetNode ( _NonUpdatablePairLists )
        pairList = pNode.Get ( "mmmm14", None )
        if pairList is None:
            pairList = SelfPairList.FromSelfPairList ( target.mmState.interactions14 ,
                                                       len ( target.atoms )          ,
                                                       target.mmState.mmAtoms        ,
                                                       target.freeAtoms              )
            pNode.mmmm14 = pairList
            sNode        = target.scratch.Get ( _PairListStatistics )
            sNode["MM/MM 1-4 Pairs"] = float ( len ( pairList ) )
        if pairList is not None:
            gradients3 = target.scratch.Get ( "gradients3", None )
            scale      = ( target.mmModel.electrostaticScale14 / self.dielectric )
            ( eElectrostatic, eLennardJones ) = self.pairwiseInteraction.MMMMEnergy ( target.mmState.charges        ,
                                                                                      target.mmState.charges        ,
                                                                                      target.mmState.ljTypeIndices  ,
                                                                                      target.mmState.ljTypeIndices  ,
                                                                                      target.mmState.ljParameters14 ,
                                                                                      scale                         ,
                                                                                      1.0                           ,
                                                                                      target.coordinates3           ,
                                                                                      target.coordinates3           ,
                                                                                      pairList                      ,
                                                                                      gradients3                    ,
                                                                                      gradients3                    )
            energies.update ( { "MM/MM 1-4 Electrostatic" : eElectrostatic ,
                                "MM/MM 1-4 Lennard-Jones" : eLennardJones  } )
        return energies

    def EnergyClosures ( self, target ):
        """Return energy closures."""
        def a ( ):
            target.scratch.energyTerms.update ( self.Energy ( target ) )
        def b ( ):
            target.scratch.energyTerms.update ( self.Energy14 ( target ) )
        closures = super ( NBModelFull, self ).EnergyClosures ( target )
        closures.extend ( [ ( EnergyClosurePriority.IndependentEnergyTerm, a, "MM/MM NB Evaluation"     ) ,
                            ( EnergyClosurePriority.IndependentEnergyTerm, b, "MM/MM 1-4 NB Evaluation" ) ] )
        return closures

    def QCMMModels ( self, qcModel = None, withSymmetry = False ):
        """Default companion QC/MM models for the model."""
        models = { "qcmmLennardJones" : QCMMLennardJonesModelFull }
        if   isinstance ( qcModel, QCModelDFT  ): model = QCMMElectrostaticModelDensityFullGaussianBasis
        elif isinstance ( qcModel, QCModelMNDO ): model = QCMMElectrostaticModelDensityFullMNDO
        else:                                     model = QCMMElectrostaticModelMultipoleFull
        models["qcmmElectrostatic"] = model
        return models

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
