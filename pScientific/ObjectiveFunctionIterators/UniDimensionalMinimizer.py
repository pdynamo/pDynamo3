"""Base classes for unidimensional minimization algorithms."""

# . This may be split as in the multidimensional case.

import sys, traceback

from  pCore                          import AttributableObject              , \
                                            logFile                         , \
                                            LogFileActive                   , \
                                            SummarizableObject
from .ObjectiveFunction              import UniDimensionalObjectiveFunction
from .ObjectiveFunctionIteratorError import ObjectiveFunctionIteratorError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class UniDimensionalMinimizerState ( AttributableObject ):
    """Base class for unidimensional minimization algorithms."""

    _attributable           = dict ( AttributableObject._attributable )
    _objectiveFunctionClass = UniDimensionalObjectiveFunction
    _attributable.update ( { "error"                 : None  ,
                             "f"                     : None  ,
                             "g"                     : None  ,
                             "isConverged"           : False ,
                             "log"                   : None  ,
                             "numberOfFunctionCalls" : 0     ,
                             "numberOfIterations"    : 0     ,
                             "objectiveFunction"     : None  ,
                             "statusMessage"         : None  ,
                             "stepType"              : ""    ,
                             "table"                 : None  ,
                             "x"                     : None  } )

    def ExtractSurrogateData ( self, surrogate ):
        """Extract additional data from a surrogate objective function."""
        pass

    def Finalize ( self ):
        """Finalization."""
        # . Create the report.
        report = { "Converged"      : self.isConverged           ,
                   "Function Calls" : self.numberOfFunctionCalls ,
                   "Function Value" : self.f                     ,
                   "Gradient"       : self.g                     ,
                   "Iterations"     : self.numberOfIterations    ,
                   "Status Message" : self.statusMessage         ,
                   "Variable"       : self.x                     }
        if self.error is not None: report["Error"] = self.error
        return report

    @classmethod
    def FromObjectiveFunction ( selfClass, objectiveFunction ):
        """Constructor given an objective function."""
        # . Create the object.
        self = selfClass ( )
        # . Check the objective function.
        if not isinstance ( objectiveFunction, self.__class__._objectiveFunctionClass ): raise ObjectiveFunctionIteratorError ( "Invalid objective function." )
        self.objectiveFunction = objectiveFunction
        # . Algorithm-specific set up.
        self.SetUp ( )
        return self

    def HandleError ( self, error ):
        """Handle an error."""
        self.error = error.args[0]
        traceback.print_exc ( file = sys.stdout )

    def SetUp ( self ):
        """Set up the state."""
        pass

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class UniDimensionalMinimizer ( SummarizableObject ):
    """Base class for unidimensional minimization algorithms."""

    _attributable = dict ( SummarizableObject._attributable )
    _classLabel   = "Uni-Dimensional Minimizer"
    _stateObject  = UniDimensionalMinimizerState
    _summarizable = dict ( SummarizableObject._summarizable )
    _attributable.update ( { "gradientTolerance" : 1.0e-3               ,
                             "logFrequency"      : 1                    ,
                             "maximumIterations" : 0                    } )
    _summarizable.update ( { "gradientTolerance" : "Gradient Tolerance" ,
                             "logFrequency"      : "Log Frequency"      ,
                             "maximumIterations" : "Maximum Iterations" } )

    def Continue ( self, state ):
        """Check to see if iterations should continue."""
        state.isConverged = ( state.g <= self.gradientTolerance )
        if   state.isConverged:                                  state.statusMessage = "Minimization converged."
        elif state.numberOfIterations >= self.maximumIterations: state.statusMessage = "Too many iterations."
        elif state.error is not None:                            state.statusMessage = "Minimization error: " + state.error
        else: state.statusMessage = None
        return ( state.statusMessage is None )

    def FunctionGradient ( self, state ):
        """Evaluate the function and its gradient."""
        ( state.f, state.g ) = state.objectiveFunction.FunctionGradient ( state.x )
        state.numberOfFunctionCalls += 1

    def Initialize ( self, state ):
        """Initialization before iteration."""
        try:
            self.FunctionGradient ( state )
        except Exception as error:
            state.HandleError ( error )
        state.stepType = "I"

    def Iterate ( self, objectiveFunction, log = logFile ):
        """Apply the algorithm to a function."""
        state = self.StateFromObjectiveFunction ( objectiveFunction )
        self.LogStart     ( state, log = log )
        self.Initialize   ( state )
        self.LogIteration ( state )
        while ( self.Continue ( state ) ):
            self.Iteration    ( state )
            self.LogIteration ( state )
        self.LogStop ( state )
        return state.Finalize ( )

    def Iteration ( self, state ):
        """Perform an iteration."""
        try:
            pass
        except Exception as error:
            state.HandleError ( error )
        state.numberOfIterations += 1

    def LogIteration ( self, state ):
        """Log an iteration."""
        state.objectiveFunction.LogIteration ( state.numberOfIterations )
        if ( state.table is not None ) and ( state.numberOfIterations % self.logFrequency == 0 ):
            state.table.Entry ( "{:d}".format ( state.numberOfIterations ) )
            state.table.Entry ( state.stepType )
            state.table.Entry ( "{:20.8f}".format ( state.f ) )
            state.table.Entry ( "{:20.8f}".format ( state.g ) )
            state.table.Entry ( "{:20.8f}".format ( state.x ) )

    def LogStart ( self, state, log = logFile ):
        """Start logging."""
        state.objectiveFunction.LogStart ( )
        if ( self.logFrequency > 0 ) and LogFileActive ( log ):
            state.log   = log
            state.table = log.GetTable ( columns = [ 6, 6, 20, 20, 20 ] )
            state.table.Start ( )
            state.table.Heading ( "Iteration", columnSpan = 2 )
            state.table.Heading ( "Function" )
            state.table.Heading ( "Gradient" )
            state.table.Heading ( "Variable" )

    def LogStop ( self, state ):
        """Stop logging."""
        state.objectiveFunction.LogStop ( )
        if state.log is not None:
            state.table.Stop ( )
            if state.statusMessage is not None: state.log.Paragraph ( state.statusMessage )
            state.log.SummaryOfItems ( [ ( "Iterations",     "{:d}".format ( state.numberOfIterations    ) ) ,
                                           ( "Function Calls", "{:d}".format ( state.numberOfFunctionCalls ) ) ] ,
                                           title = self.__class__._classLabel + " Statistics" )

    def Restart ( self, state ):
        """Restart the minimization."""
        state.stepType = "R"

    def StateFromObjectiveFunction ( self, objectiveFunction ):
        """Set up the state."""
        return self.__class__._stateObject.FromObjectiveFunction ( objectiveFunction )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
