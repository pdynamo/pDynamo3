"""The base class for QC/MM Lennard-Jones models."""


from   pCore        import logFile       , \
                           LogFileActive
from  .NBModelError import NBModelError
from ..EnergyModel  import EnergyModel

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class QCMMLennardJonesModel ( EnergyModel ):
    """The base class for QC/MM Lennard-Jones models."""

    _attributable             = dict ( EnergyModel._attributable )
    _classLabel               = "QC/MM Lennard-Jones Model"
    _pairwiseInteractionClass = None
    _summarizable             = dict ( EnergyModel._summarizable )
    _attributable.update ( { "pairwiseInteraction" : None } )
    _summarizable.update ( { "pairwiseInteraction" : None } )

    def _CheckOptions ( self ):
        """Check options."""
        pwiClass = self.__class__._pairwiseInteractionClass
        if pwiClass is not None:
            if self.pairwiseInteraction is None:
                self.pairwiseInteraction = pwiClass ( )
            elif not isinstance ( self.pairwiseInteraction, pwiClass ):
                raise TypeError ( "Invalid pairwise interaction attribute." ) 

    def ClearScratch ( self, scratch ): scratch.Clear ( )

    def EnergyClosures ( self, target ):
        """Return energy closures."""
        return []

    def Energy ( self, target ):
        """Energy and gradients."""
        pass

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
