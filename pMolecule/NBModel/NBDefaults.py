"""Various defaults for the NB modules."""

from  pScientific                   import Units
from  pScientific.Geometry3         import PairListGenerator
from .NBModelError                  import NBModelError
from .PairwiseInteractionABFS       import PairwiseInteractionABFS
from .PairwiseInteractionSplineABFS import PairwiseInteractionSplineABFS , \
                                           SplineModel

#===================================================================================================================================
# . Cut-off model defaults.
#===================================================================================================================================
_DampingCutOff   =  0.5
_InnerCutOff     =  8.0
_ListCutOff      = 13.5
_OuterCutOff     = 12.0
_QCDampingCutOff =  0.1
_QCWidth         =  2.0 * Units.Length_Bohrs_To_Angstroms

def _CheckCutOffs ( self ):
    """Check the generator and pairwise interaction for consistency."""
    if ( self.generator is not None ) and ( self.pairwiseInteraction is not None ):
        buffer = self.generator.cutOff - self.pairwiseInteraction.range
        if buffer < 0.0: raise NBModelError ( "Incompatible generator/pairwise interaction cut-offs." )

def _DefaultGeneratorCutOff ( ):
    """The default generator."""
    return PairListGenerator.WithOptions ( cutOff               = _ListCutOff ,
                                           cutOffCellSizeFactor = 0.5         ,
                                           minimumCellExtent    = 2           ,
                                           minimumCellSize      = 3.0         ,
                                           minimumExtentFactor  = 1.5         ,
                                           minimumPoints        = 500         ,
                                           sortIndices          = False       ,
                                           useGridByCell        = True        )

def _DefaultPairwiseInteractionABFS ( ):
    """The default pairwise interaction."""
    return PairwiseInteractionABFS.WithOptions ( dampingCutOff = _DampingCutOff ,
                                                 innerCutOff   = _InnerCutOff   ,
                                                 outerCutOff   = _OuterCutOff   )

def _DefaultPairwiseInteractionSplineABFS ( ):
    """The default pairwise interaction."""
    return PairwiseInteractionSplineABFS.WithOptions ( dampingCutOff = _DampingCutOff ,
                                                       innerCutOff   = _InnerCutOff   ,
                                                       outerCutOff   = _OuterCutOff   )

def _DefaultPairwiseInteractionSplineABFSQCMM ( ):
    """The default pairwise interaction."""
    return PairwiseInteractionSplineABFS.WithOptions ( dampingCutOff      = _QCDampingCutOff           ,
                                                       electrostaticModel = SplineModel.Delta_Gaussian ,
                                                       innerCutOff        = _InnerCutOff               ,
                                                       outerCutOff        = _OuterCutOff               ,
                                                       width1             = _QCWidth                   ,
                                                       width2             = _QCWidth                   )

def _DefaultPairwiseInteractionSplineABFSQCQC ( ):
    """The default pairwise interaction."""
    return PairwiseInteractionSplineABFS.WithOptions ( dampingCutOff      = _QCDampingCutOff              ,
                                                       electrostaticModel = SplineModel.Gaussian_Gaussian ,
                                                       innerCutOff        = _InnerCutOff                  ,
                                                       outerCutOff        = _OuterCutOff                  ,
                                                       width1             = _QCWidth                      ,
                                                       width2             = _QCWidth                      )

#===================================================================================================================================
# . Keys.
#===================================================================================================================================
# . Attribute names in scratch.
_CenteringTranslation3 = "centeringTranslation3"
_ImageScanData         = "imageScanData"
_MMGrid                = "mmGrid"
_MMOccupancy           = "mmOccupancy"
_NonUpdatablePairLists = "nonUpdatablePairLists"
_PairListStatistics    = "pairListStatistics"
_QCGrid                = "qcGrid"
_QCOccupancy           = "qcOccupancy"
_UpdatablePairLists    = "updatablePairLists"
_UpdateChecker         = "updateChecker"

# . Statistics names in pairListStatistics.
_NumberOfCalls         = "Number of Calls"
_NumberOfUpdates       = "Number of Updates"

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
