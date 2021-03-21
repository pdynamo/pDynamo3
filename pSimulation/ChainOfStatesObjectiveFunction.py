"""Chain-of-states objective function."""

# . This is a misnomer as an objective function cannot be defined for many of these methods.
# . Nevertheless, the object serves as an appropriate interface to the system being optimized.

from pMolecule import MultiLayerSystemGeometryObjectiveFunction

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class ChainOfStatesObjectiveFunction ( MultiLayerSystemGeometryObjectiveFunction ):
    """The chain-of-states objective function."""

    _attributable = dict ( MultiLayerSystemGeometryObjectiveFunction._attributable )
    _attributable.update ( { "imageTrajectory" : None ,
                             "numberOfImages"  : 0    } )

    def DumpImage ( self, image ):
        """Dump an image."""
        self.imageTrajectory.WriteOwnerData ( index = image )

    def FinalizeImages ( self ):
        """Finalize image data."""
        self.imageTrajectory.Close ( )

    def GetGradients ( self, gradients ):
        """Get the gradients (without applying linear constaints)."""
        self.system.scratch.gradients3.CopyToToVector ( gradients, selection = self.freeAtoms )
        if self.variableWeights is not None: gradients.Divide ( self.variableWeights )

    def InitializeImages ( self, imageTrajectory ):
        """Initialize image data given a trajectory."""
        self.imageTrajectory = imageTrajectory
        self.numberOfImages  = len ( imageTrajectory )
        if ( self.numberOfImages <= 2 ): raise ValueError ( "Invalid number of images on trajectory: {:d}.".format ( self.numberOfImages ) )

    def LoadImage ( self, image ):
        """Load an image."""
        self.imageTrajectory.RestoreOwnerData ( index = image )
        # . A fudge.
        if self.freeAtoms is None: self.iCoordinates3 = self.system.coordinates3.iterator
        else:                      self.iCoordinates3 = self.system.coordinates3.RowIterator ( selection = self.freeAtoms )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
