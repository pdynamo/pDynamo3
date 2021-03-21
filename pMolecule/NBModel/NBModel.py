"""The base class for MM/MM interaction non-bonding models."""

from   collections  import defaultdict
from   pCore        import logFile               , \
                           LogFileActive
from  .NBDefaults   import _NumberOfCalls        , \
                           _NumberOfUpdates      , \
                           _PairListStatistics
from  .NBModelError import NBModelError
from ..EnergyModel  import EnergyClosurePriority , \
                           EnergyModel           , \
                           EnergyModelPriority

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class NBModel ( EnergyModel ):
    """The base class for MM/MM non-bonding models."""

    _attributable             = dict ( EnergyModel._attributable )
    _classLabel               = "NB Model"
    _pairwiseInteractionClass = None
    _summarizable             = dict ( EnergyModel._summarizable )
    _attributable.update ( { "dielectric"          : 1.0  ,
                             "pairwiseInteraction" : None } )
    _summarizable.update ( { "dielectric"          : ( "Dielectric", "{:.3f}" ) ,
                             "pairwiseInteraction" : None                     } )

    def BuildModel ( self, target, assignQCMMModels = True ):
        """Build the model."""
        state = super ( NBModel, self ).BuildModel ( target )
        if target.mmModel is None: raise NBModelError ( "An MM model is missing." )
        # . Assign QC/MM models.
        # . It is assumed that previous NB and QC/MM models have already been popped.
        if assignQCMMModels               and \
           ( target.qcModel is not None ) and \
           ( len ( target.qcState.qcAtoms ) <= len ( target.atoms ) ):
            withSymmetry = ( target.symmetryParameters is not None )
            for ( key, valueClass ) in self.QCMMModels ( qcModel = target.qcModel, withSymmetry = withSymmetry ).items ( ):
                target._AddEnergyModel ( key, valueClass.WithDefaults ( ), priority = EnergyModelPriority.QCMMModel )

    def _CheckOptions ( self ):
        """Check options."""
        pwiClass = self.__class__._pairwiseInteractionClass
        if pwiClass is not None:
            if self.pairwiseInteraction is None:
                self.pairwiseInteraction = pwiClass ( )
            elif not isinstance ( self.pairwiseInteraction, pwiClass ):
                raise TypeError ( "Invalid pairwise interaction attribute." ) 

    def ClearScratch ( self, scratch ):
        """Clear scratch."""
        scratch.Clear ( )

    def EnergyClosures ( self, target ):
        """Return energy closures."""
        def a ( ): self.EnergyInitialize ( target )
        return [ ( EnergyClosurePriority.NBInitialization, a, "NB Initialization" ) ]

    def EnergyInitialize ( self, target ):
        """Energy initialization."""
        target.scratch.GetSet ( _PairListStatistics, defaultdict, float )

    def QCMMModels ( self, qcModel = None, withSymmetry = False ):
        """Default companion QC/MM models for the model."""
        return {}

    def StatisticsSummary ( self, target, log = logFile ):
        """Summary of pairlist statistics."""
        node = target.scratch.Pop ( _PairListStatistics )
        if LogFileActive ( log ) and ( node is not None ):
            calls   = int ( node.pop ( _NumberOfCalls  , 0.0 ) )
            updates = int ( node.pop ( _NumberOfUpdates, 0.0 ) )
            u       = float ( max ( updates, 1 ) )
            pairs   = []
            for ( key, value ) in node.items ( ):
                if key.startswith ( "<" ):
                    if updates > 1: pairs.append ( ( key, "{:.1f}".format ( value / u ) ) )
                else:
                    pairs.append ( ( key, "{:d}".format  ( int ( value ) ) ) )
            items = []
            if calls   > 0: items.append ( ( _NumberOfCalls  , "{:d}".format ( calls   ) ) )
            if updates > 0: items.append ( ( _NumberOfUpdates, "{:d}".format ( updates ) ) )
            for ( key, value ) in sorted ( pairs ): items.append ( ( key, value ) )
            if updates > 1: items.append ( ( "Calls per Update" , "{:.1f}".format ( float ( calls ) / u ) ) )
            log.SummaryOfItems ( items, order = False, title = "NB Model Pairlist Statistics Summary" )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
