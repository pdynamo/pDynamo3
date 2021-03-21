"""Determine a steepest descent path."""

import math

from   pCore                          import logFile                        , \
                                             LogFileActive
from  .ObjectiveFunctionIterator      import ObjectiveFunctionIterator      , \
                                             ObjectiveFunctionIteratorState
from  .ObjectiveFunctionIteratorError import ObjectiveFunctionIteratorError
from ..Arrays                         import Array
from ..LinearAlgebra                  import EigenPairs

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class SteepestDescentPathFinderState ( ObjectiveFunctionIteratorState ):
    """Class for the state of a steepest descent path calculation."""

    _attributable = dict ( ObjectiveFunctionIteratorState._attributable )
    _attributable.update ( { "direction"             :  0.0 ,
                             "d"                     : None ,
                             "d0"                    : None ,
                             "g"                     : None ,
                             "isFirstBranch"         : True ,
                             "numberOfFunctionCalls" :    0 ,
                             "pathPosition"          :  0.0 ,
                             "saddleStep"            :  0.0 ,
                             "stepLength"            :  0.0 ,
                             "x0"                    : None } )

    def Finalize ( self ):
        """Finalization."""
        report = { "Function Calls" : self.numberOfFunctionCalls }
        report.update ( super ( SteepestDescentPathFinderState, self ).Finalize ( ) )
        return report

    def SetUp ( self ):
       """Set up the state."""
       n = self.numberOfVariables
       if self.d  is None: self.d  = Array.WithExtent ( n ) ; self.d.Set  ( 0.0 )
       if self.d0 is None: self.d0 = Array.WithExtent ( n ) ; self.d0.Set ( 0.0 )
       if self.g  is None: self.g  = Array.WithExtent ( n ) ; self.g.Set  ( 0.0 )
       if self.x0 is None: self.x0 = Array.WithExtent ( n ) ; self.x0.Set ( 0.0 )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class SteepestDescentPathFinder ( ObjectiveFunctionIterator ):
    """Class for determining a steepest descent path using the simplest possible stepping algorithm."""

    _attributable = dict ( ObjectiveFunctionIterator._attributable )
    _classLabel   = "Steepest Descent Path Finder"
    _stateObject  = SteepestDescentPathFinderState
    _summarizable = dict ( ObjectiveFunctionIterator._summarizable )
    _attributable.update ( { "fromSaddle"   : True            ,
                             "functionStep" : 0.2             ,
                             "pathStep"     : 0.025           } )
    _summarizable.update ( { "fromSaddle"   : "From Saddle"   ,
                             "functionStep" : "Function Step" ,
                             "pathStep"     : "Path Step"     } )

    def Continue ( self, state ):
        """Check to see if iterations should continue."""
        toBeContinued = super ( SteepestDescentPathFinder, self ).Continue ( state )
        if ( not toBeContinued ) and self.fromSaddle and state.isFirstBranch:
            state.isFirstBranch      = False
            state.numberOfIterations = 0
            state.pathPosition       = 0.0
            toBeContinued            = True
            self.GetFirstStep ( state )
        return toBeContinued

    def GetFirstStep ( self, state ):
        """Get the initial step."""
        try:
            # . The first point is a saddle point.
            if self.fromSaddle:
                # . First time around.
                if state.isFirstBranch:
                    # . Allocation.
                    n = state.numberOfVariables
                    eigenValues  = Array.WithExtent  ( n    ) ; eigenValues.Set  ( 0.0 )
                    eigenVectors = Array.WithExtents ( n, n ) ; eigenVectors.Set ( 0.0 )
                    # . Get the starting Hessian.
                    if not hasattr ( state.objectiveFunction, "StartingHessian" ): raise ObjectiveFunctionIteratorError ( "Object function has no Hessian evaluator." )
                    hessian = state.objectiveFunction.StartingHessian ( state.x )
                    # . Get the eigenValues and eigenVectors.
                    EigenPairs ( hessian, eigenValues, eigenVectors )
                    # . Get the indices of the negative eigenValues.
                    negative = [ i for ( i, e ) in enumerate ( eigenValues ) if e < 0.0 ]
                    if len ( negative ) != 1: raise ObjectiveFunctionIteratorError ( "The starting point is not a first-order saddle point." )
                    # . Get the direction of the step from the saddle as the negative eigenVector.
                    for i in range ( state.numberOfVariables ): state.d[i] = eigenVectors[i,negative[0]]
                    # . Apply linear constraints.
                    state.objectiveFunction.ApplyLinearConstraints ( state.d )
                    state.d.Normalize ( )
                    # . Get the initial stepsize - the negative direction first.
                    state.direction  = - 1.0
                    state.saddleStep = math.sqrt ( - 2.0 * self.functionStep / eigenValues[negative[0]] )
                    # . Scale the step.
                    state.d.Scale ( state.direction * state.saddleStep )
                    state.stepLength = state.saddleStep
                    # . Save the step and the first point.
                    state.d.CopyTo ( state.d0 )
                    state.x.CopyTo ( state.x0 )
                # . Second time around.
                else:
                    state.direction  = + 1.0
                    state.stepLength = state.saddleStep
                    state.d0.CopyTo ( state.d )
                    state.d.Scale   ( - 1.0   )
                    state.x0.CopyTo ( state.x )
            # . The first point is not a saddle point.
            else:
                state.isFirstBranch = False
                state.direction     = + 1.0
                self.GetStep ( state )
        except Exception as error:
            state.error = error.args[0]

    def GetStep ( self, state ):
        """Get a step."""
        state.g.CopyTo ( state.d )
        state.d.Normalize ( )
        state.d.Scale ( - self.pathStep )
        state.stepLength = self.pathStep

    def Initialize ( self, state ):
        """Initialization before iteration."""
        try:
            # . First function evaluation.
            state.objectiveFunction.VariablesGet ( state.x )
            state.f = state.objectiveFunction.FunctionGradients ( state.x, state.g )
            state.numberOfFunctionCalls += 1
            # . Get the initial step.
            self.GetFirstStep ( state )
        except Exception as error:
            state.error = error.args[0]

    def Iteration ( self, state ):
        """Perform an iteration using the simplest possible step."""
        state.numberOfIterations += 1
        state.pathPosition       += state.direction * state.stepLength
        # . Generate the new variables.
        state.x.Add ( state.d )
        # . Calculate the new function and gradients.
        state.f = state.objectiveFunction.FunctionGradients ( state.x, state.g )
        state.numberOfFunctionCalls += 1
        # . Determine the new step.
        self.GetStep ( state )

    def LogIteration ( self, state ):
        """Log an iteration."""
        state.objectiveFunction.LogIteration ( state.numberOfIterations, identifier = state.pathPosition )
        if ( state.table is not None ) and ( state.numberOfIterations % self.logFrequency == 0 ):
            state.table.Entry ( "{:d}"    .format ( state.numberOfIterations   ) )
            state.table.Entry ( "{:20.8f}".format ( state.pathPosition         ) )
            state.table.Entry ( "{:20.8f}".format ( state.f                    ) )
            state.table.Entry ( "{:20.8f}".format ( state.g.RootMeanSquare ( ) ) )

    def LogStart ( self, state, log = logFile ):
        """Start logging."""
        state.objectiveFunction.LogStart ( )
        if ( self.logFrequency > 0 ) and LogFileActive ( log ):
            state.log   = log
            state.table = log.GetTable ( columns = [ 10, 20, 20, 20 ] )
            state.table.Start ( )
            state.table.Heading ( "Iteration"     )
            state.table.Heading ( "Path Position" )
            state.table.Heading ( "Function"      )
            state.table.Heading ( "RMS Gradient"  )

    def LogStop ( self, state ):
        """Stop logging."""
        state.objectiveFunction.LogStop ( )
        if state.log is not None:
            state.table.Stop ( )
            if state.statusMessage is not None: state.log.Paragraph ( state.statusMessage )
            state.log.SummaryOfItems ( [ ( "Iterations",     "{:d}".format ( state.numberOfIterations    ) ) ,
                                           ( "Function Calls", "{:d}".format ( state.numberOfFunctionCalls ) ) ] ,
                                         title = self.__class__._classLabel + " Statistics" )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
