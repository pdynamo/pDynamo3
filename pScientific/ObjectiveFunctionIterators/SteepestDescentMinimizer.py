"""Class for performing steepest-descent multidimensional minimization."""

from  .MultiDimensionalMinimizer import MultiDimensionalMinimizer      , \
                                        MultiDimensionalMinimizerState
from ..Arrays                    import Array

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Defaults for various options.
_DefaultMaximumStepSize = 1.0
_DefaultMinimumStepSize = 1.0e-6
_DefaultScaleDown       = 0.5
_DefaultScaleUp         = 1.2
_DefaultStepSize        = 1.0e-2

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class SteepestDescentMinimizerState ( MultiDimensionalMinimizerState ):
    """Steepest descent minimizer state."""

    _attributable = dict ( MultiDimensionalMinimizerState._attributable )
    _attributable.update ( { "fBest"  : 0.0  ,
                             "fOld"   : 0.0  ,
                             "gBest"  : None ,
                             "g2Best" : 0.0  ,
                             "xBest"  : None } )

    def Finalize ( self ):
        """Finalization."""
        # . Make sure the best point is saved.
        self.objectiveFunction.VariablesPut ( self.xBest )
        # . Finish finalizing.
        return super ( SteepestDescentMinimizerState, self ).Finalize ( )

    def SetUp ( self ):
        """Set up the state."""
        super ( SteepestDescentMinimizerState, self ).SetUp ( )
        self.gBest = Array.WithExtent ( self.numberOfVariables ) ; self.gBest.Set ( 0.0 )
        self.xBest = Array.WithExtent ( self.numberOfVariables ) ; self.xBest.Set ( 0.0 )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class SteepestDescentMinimizer ( MultiDimensionalMinimizer ):
    """Steepest descent multidimensional minimizer."""

    _attributable = dict ( MultiDimensionalMinimizer._attributable )
    _classLabel   = "Steepest-Descent Minimizer"
    _stateObject  = SteepestDescentMinimizerState
    _summarizable = dict ( MultiDimensionalMinimizer._summarizable )
    _attributable.update ( { "maximumStepSize" : _DefaultMaximumStepSize ,
                             "minimumStepSize" : _DefaultMinimumStepSize ,
                             "scaleDown"       : _DefaultScaleDown       ,
                             "scaleUp"         : _DefaultScaleUp         ,
                             "stepSize"        : _DefaultStepSize        } )
    _summarizable.update ( { "maximumStepSize" : "Maximum Step Size"     ,
                             "minimumStepSize" : "Minimum Step Size"     , 
                             "scaleDown"       : "Scale Down"            ,
                             "scaleUp"         : "Scale Up"              ,
                             "stepSize"        : "Step Size"             } )

    def Initialize ( self, state ):
        """Initialization before iteration."""
        super ( SteepestDescentMinimizer, self ).Initialize ( state )
        # . Initialization.
        state.alpha  = self.stepSize
        state.fOld   = state.f
        # . Save the best point.
        state.fBest  = state.f
        state.g2Best = state.g.DotSelf ( )
        state.g.CopyTo ( state.gBest )
        state.x.CopyTo ( state.xBest )

    def DetermineStep ( self, state ):
        """Get the step."""
        # . Set the step type.
        state.stepType = "N"
        # . Adjust the step.
        if state.numberOfIterations > 0:
            # . Calculate the change in function value.
            deltaF    = state.f - state.fOld
            state.fOld = state.f
            # . The function value has decreased or remained the same.
            if deltaF <= 0.0:
                if deltaF < 0.0:
                    state.fBest  = state.f
                    state.g2Best = state.g.DotSelf ( )
                    state.g.CopyTo ( state.gBest )
                    state.x.CopyTo ( state.xBest )
            # . Modify the step according to the change in function value.
            if deltaF > 0.0:
                state.alpha   *= self.scaleDown
                state.stepType = "D"
                if state.alpha < self.minimumStepSize: state.error = "Step size too small."
            elif deltaF < 0.0:
                state.alpha   *= self.scaleUp
                state.stepType = "U"
                if state.alpha > self.maximumStepSize: state.alpha = self.maximumStepSize
        # . Take the step.
        state.g.CopyTo ( state.d )
        state.d.Normalize ( )
        state.x.Add ( state.d, scale = -state.alpha )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
