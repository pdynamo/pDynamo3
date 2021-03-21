"""Defines a simple Monte Carlo NB model."""

from  .NBModel                       import NBModel
from  .PairwiseInteractionMonteCarlo import PairwiseInteractionMonteCarlo
from ..EnergyModel                   import EnergyClosurePriority

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class NBModelMonteCarlo ( NBModel ):
    """Define a Monte Carlo NB model."""

    # . Defaults.
    _classLabel               = "Monte Carlo NB Model"
    _pairwiseInteractionClass = PairwiseInteractionMonteCarlo

    def Energy ( self, target ):
        """Energy."""
        ( eElectrostatic, eLennardJones ) = self.pairwiseInteraction.MMMMEnergy ( target.mmState.charges             ,
                                                                                  target.mmState.ljTypeIndices       ,
                                                                                  target.mmState.ljParameters        ,
                                                                                  ( 1.0 / self.dielectric )          ,
                                                                                  target.connectivity.isolateIndices ,
                                                                                  target.coordinates3                ,
                                                                                  target.symmetryParameters          )
        return { "MM/MM Electrostatic" : eElectrostatic ,
                 "MM/MM Lennard-Jones" : eLennardJones  }

    def EnergyClosures ( self, target ):
        """Return energy closures."""
        def a ( ):
            target.scratch.energyTerms.update ( self.Energy ( target ) )
        closures = super ( NBModelMonteCarlo, self ).EnergyClosures ( target )
        closures.append ( ( EnergyClosurePriority.IndependentEnergyTerm, a, "MM/MM NB Evaluation" ) )
        return closures

    def EnergyOfIsolate ( self, target, isolate ):
        """Energy of a single isolate."""
        ( eElectrostatic, eLennardJones ) = self.pairwiseInteraction.MMMMIsolateEnergy ( isolate                            ,
                                                                                         target.mmState.charges             ,
                                                                                         target.mmState.ljTypeIndices       ,
                                                                                         target.mmState.ljParameters        ,
                                                                                         ( 1.0 / self.dielectric )          ,
                                                                                         target.connectivity.isolateIndices ,
                                                                                         target.coordinates3                ,
                                                                                         target.symmetryParameters          )
        return { "MM/MM Isolate {:d} Electrostatic".format ( isolate ) : eElectrostatic ,
                 "MM/MM Isolate {:d} Lennard-Jones".format ( isolate ) : eLennardJones  }

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
