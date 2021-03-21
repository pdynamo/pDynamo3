"""Fast inertial relaxation engine (FIRE) minimizer."""

from  .MultiDimensionalMinimizer import MultiDimensionalMinimizer      , \
                                        MultiDimensionalMinimizerState
from ..Arrays                    import Array

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class FIREMinimizerState ( MultiDimensionalMinimizerState ):
    """FIRE minimizer state."""

    _attributable = dict ( MultiDimensionalMinimizerState._attributable )
    _attributable.update ( { "alpha"          : 0.0  ,
                             "gToAFactor"     : 1.0  ,
                             "numberPositive" : 0    ,
                             "timeStep"       : 0.0  ,
                             "v"              : None } )

    def ExtractSurrogateData ( self, surrogate ):
        """Extract additional data from a surrogate objective function."""
        if surrogate is not None:
            self.gToAFactor = surrogate.AccelerationConversionFactor ( )

    @classmethod
    def FromObjectiveFunction ( selfClass, objectiveFunction ):
        """Constructor given an objective function."""
        self            = super ( FIREMinimizerState, selfClass ).FromObjectiveFunction ( objectiveFunction )
        self.gToAFactor = self.objectiveFunction.AccelerationConversionFactor ( )
        return self

    def SetUp ( self ):
       """Set up the state."""
       super ( FIREMinimizerState, self ).SetUp ( )
       self.v = Array.WithExtent ( self.numberOfVariables ) ; self.v.Set ( 0.0 )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class FIREMinimizer ( MultiDimensionalMinimizer ):
    """FIRE multidimensional minimizer."""

    _attributable = dict ( MultiDimensionalMinimizer._attributable )
    _classLabel   = "FIRE Minimizer"
    _stateObject  = FIREMinimizerState
    _summarizable = dict ( MultiDimensionalMinimizer._summarizable )
    _attributable.update ( { "alpha"              : 0.1                     ,
                             "alphaDownFactor"    : 0.99                    ,
                             "maximumStep"        : 0.01                    ,
                             "maximumTimeStep"    : 1.0                     ,
                             "stepLatencyFactor"  : 5                       ,
                             "timeStep"           : 0.1                     ,
                             "timeStepDownFactor" : 0.5                     ,
                             "timeStepUpFactor"   : 1.1                     } )
    _summarizable.update ( { "alpha"              : "Alpha"                 ,
                             "alphaDownFactor"    : "Alpha Down Factor"     ,
                             "maximumStep"        : "Maximum Step"          ,
                             "maximumTimeStep"    : "Maximum Time Step"     ,
                             "stepLatencyFactor"  : "Step Latency Factor"   ,
                             "timeStep"           : "Time Step"             ,
                             "timeStepDownFactor" : "Time Step Down Factor" ,
                             "timeStepUpFactor"   : "Time Step Up Factor"   } )

    def DetermineStep ( self, state ):
        """Get the step."""
        # . A velocity-Verlet algorithm is used to solve for x and v
        # . as Euler integration does not appear to be sufficient.
        # . Complete the calculation of v from the previous iteration.
        if state.numberOfIterations > 0:
            dt = state.timeStep
            state.v.Add ( state.g, scale = 0.5 * dt * state.gToAFactor )
            # . Perform the FIRE procedure.
            # . Negative or zero power.
            if state.g.Dot ( state.v ) >= 0.0:
                state.alpha          = self.alpha
                state.numberPositive = 0
                state.stepType       = "-"
                state.timeStep      *= self.timeStepDownFactor
                state.v.Set ( 0.0 )
            # . Positive power.
            else:
                # . Modify v.
                gNorm2 = state.g.Norm2 ( )
                vNorm2 = state.v.Norm2 ( )
                state.v.Scale ( 1.0 - state.alpha )
                if gNorm2 != 0.0: state.v.Add ( state.g, scale = - state.alpha * vNorm2 / gNorm2 )
                # . Adjust factors.
                if state.numberPositive > self.stepLatencyFactor:
                    state.alpha   *= self.alphaDownFactor
                    state.timeStep = min ( state.timeStep * self.timeStepUpFactor, self.maximumTimeStep )
                state.numberPositive += 1
                state.stepType        = "+"
        else: state.stepType = "0"
        # . Calculate d and v at the half-step.
        dt = state.timeStep
        state.d.Set ( 0.0 )
        state.d.Add ( state.v, scale =       dt                       )
        state.d.Add ( state.g, scale = 0.5 * dt**2 * state.gToAFactor )
        state.v.Add ( state.g, scale = 0.5 * dt    * state.gToAFactor )
        # . Apply linear constraints - unnecessary if already applied to g.
#        state.objectiveFunction.ApplyLinearConstraints ( state.d )
        # . Check the step length - this appears to be important for this optimizer.
        step = state.d.RootMeanSquare ( )
        if step > self.maximumStep: state.d.Scale ( self.maximumStep / step )
        # . Increment the variables.
        state.x.Add ( state.d )

    def Initialize ( self, state  ):
        """Initialization before iteration."""
        self.Restart ( state )
        super ( FIREMinimizer, self ).Initialize ( state )

    def Restart ( self, state ):
        """Restart the minimization."""
        state.alpha          = self.alpha
        state.numberPositive = 0
        state.stepType       = "R"
        state.timeStep       = self.timeStep
        state.v.Set ( 0.0 )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
