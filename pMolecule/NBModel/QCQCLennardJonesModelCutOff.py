"""Defines a QC/QC Lennard-Jones model that employs cut-offs."""

from   pCore                       import logFile                , \
                                          LogFileActive
from  .ImagePairListContainer      import ImagePairListContainer
from  .NBDefaults                  import _ImageScanData         , \
                                          _MMGrid                , \
                                          _MMOccupancy           , \
                                          _PairListStatistics    , \
                                          _UpdatablePairLists
from  .QCMMLennardJonesModelCutOff import QCMMLennardJonesModelCutOff
from ..EnergyModel                 import EnergyClosurePriority

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCQCLennardJonesModelCutOff ( QCMMLennardJonesModelCutOff ):
    """Define a cut-off QC/QC Lennard-Jones model."""

    _classLabel = "CutOff QC/QC Lennard-Jones Model"

    def EnergyClosures ( self, target ):
        """Return energy closures."""
        def a ( ): target.scratch.energyTerms.update ( self.EnergyImage ( target ) )
        return [ ( EnergyClosurePriority.IndependentEnergyTerm , a, "QC/QC Image Lennard-Jones" ) ]

    def EnergyImage ( self, target ):
        """Image energy."""
        energies           = {}
        symmetryParameters = target.symmetryParameters
        if ( symmetryParameters is not None ) and \
           ( not symmetryParameters.IsMinimumImageConventionSatisfied ( self.pairwiseInteraction.range ) ):
            scratch      = target.scratch
            coordinates3 = scratch.Get ( "coordinates3NB", target.coordinates3 )
            pNode        = scratch.GetSetNode ( _UpdatablePairLists )
            pairList     = pNode.Get ( "qcqcImageLennardJones", None )
            if pairList is None:
                pairList = ImagePairListContainer.Constructor ( self.generator                      ,
                                                                target.qcState.pureQCAtoms          ,
                                                                target.qcState.pureQCAtoms          ,
                                                                target.freeAtoms                    ,
                                                                coordinates3                        ,
                                                                coordinates3                        ,
                                                                symmetryParameters                  ,
                                                                target.symmetry.transformations     ,
                                                                pNode.Get ( _ImageScanData , None ) ,
                                                                pNode.Get ( _MMGrid        , None ) ,
                                                                pNode.Get ( _MMOccupancy   , None ) ,
                                                                True                                )
                pNode.qcqcImageLennardJones = pairList
                sNode = scratch.Get ( _PairListStatistics )
                nI    = float ( pairList.numberOfImages )
                nP    = float ( pairList.numberOfPairs  )
                sNode[ "QC/QC Image Lennard-Jones Images" ]  = nI
                sNode[ "QC/QC Image Lennard-Jones Pairs"  ]  = nP
                sNode["<QC/QC Image Lennard-Jones Images>"] += nI
                sNode["<QC/QC Image Lennard-Jones Pairs>" ] += nP
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
                energies.update ( { "QC/QC Image LJ" : eLennardJones } )
        return energies

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
