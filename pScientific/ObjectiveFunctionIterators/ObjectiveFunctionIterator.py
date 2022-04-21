"""Objective function iterator."""

from  pCore                          import AttributableObject , \
                                            logFile            , \
                                            LogFileActive      , \
                                            SummarizableObject
from .ObjectiveFunction              import ObjectiveFunction
from .ObjectiveFunctionIteratorError import ObjectiveFunctionIteratorError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class ObjectiveFunctionIteratorState ( AttributableObject ):
    """Base class for algorithm states that require objective functions."""

    _attributable           = dict ( AttributableObject._attributable )
    _objectiveFunctionClass = ObjectiveFunction
    _attributable.update ( { "error"              : None ,
                             "f"                  : None ,
                             "log"                : None ,
                             "numberOfIterations" : 0    ,
                             "numberOfVariables"  : 0    ,
                             "objectiveFunction"  : None ,
                             "statusMessage"      : None ,
                             "table"              : None ,
                             "x"                  : None } )

    def Finalize ( self ):
        """Finalization."""
        # . Create the report.
        report = { "Function Value" : self.f                  ,
                   "Iterations"     : self.numberOfIterations }
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
        # . Set up variables.
        self.numberOfVariables = objectiveFunction.numberOfVariables
        self.x                 = objectiveFunction.VariablesAllocate ( )
        objectiveFunction.VariablesGet ( self.x )
        # . Finish up.
        self.SetUp ( )
        return self

    def SetUp ( self ):
       """Set up the state."""
       pass

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class ObjectiveFunctionIterator ( SummarizableObject ):
    """Base class for algorithms that require objective functions."""

    _attributable = dict ( SummarizableObject._attributable )
    _classLabel   = "Objective Function Iterator"
    _stateObject  = ObjectiveFunctionIteratorState
    _summarizable = dict ( SummarizableObject._summarizable )
    _attributable.update ( { "logFrequency"      :  0                   ,
                             "maximumIterations" : -1                   } )
    _summarizable.update ( { "logFrequency"      : "Log Frequency"      ,
                             "maximumIterations" : "Maximum Iterations" } )

    def Continue ( self, state ):
        """Check to see if iterations should continue."""
        if   state.error is not None:                            state.statusMessage = "Iterator error: " + state.error
        elif state.numberOfIterations >= self.maximumIterations: state.statusMessage = "Iterator iterations terminated."
        else: state.statusMessage = None
        return ( state.statusMessage is None )

    def Initialize ( self, state ):
        """Initialization before iteration."""
        pass

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

    # . All subclasses should use the same structure for this method so that errors are dealt with properly!
    def Iteration ( self, state ):
        """Perform an iteration."""
        try:
            pass
        except Exception as error:
            state.error = error.args[0]
        state.numberOfIterations += 1

    def LogIteration ( self, state ):
        """Log an iteration."""
        if ( state.log is not None ) and ( state.numberOfIterations % self.logFrequency == 0 ):
            pass

    def LogStart ( self, state, log = logFile ):
        """Start logging."""
        if ( self.logFrequency > 0 ) and LogFileActive ( log ):
            state.log = log

    def LogStop ( self, state ):
        """Stop logging."""
        if ( state.log is not None ) and ( state.statusMessage is not None ):
            state.log.Paragraph ( state.statusMessage )

    def StateFromObjectiveFunction ( self, objectiveFunction ):
        """Set up the state."""
        return self.__class__._stateObject.FromObjectiveFunction ( objectiveFunction )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
