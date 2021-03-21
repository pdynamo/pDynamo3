"""More-Thuente line searcher."""

from math                                   import fabs, sqrt
from pCore                                  import logFile
from pScientific.ObjectiveFunctionIterators import CubicMinimizerFGFG           , \
                                                   MoreThuenteLineSearcher      , \
                                                   MoreThuenteLineSearcherState , \
                                                   QuadraticExtremumFGF         , \
                                                   QuadraticExtremumGG

# . MJF.
# . Could eventually be merged with original MoreThuenteLineSearcher.
# . Alternatively modify LineSearcherObjectiveFunction (make -ve) when maximizing?

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . A large number for the upper bound of the line search.
_LargeNumber = 1.0e+20

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MoreThuenteLineSearcherWithMaxState ( MoreThuenteLineSearcherState ):
    """More-Thuente line searcher state with maximize option."""

    _attributable = dict ( MoreThuenteLineSearcherState._attributable )
    _attributable.update ( { "maximize" : False } )

    def Finalize ( self ):
        """Finalization."""
        # . Return the best point if abnormal termination.
        if ( not self.isConverged ) and ( self.f > self.fB ):
            self.objectiveFunction.RestoreBestToCurrent ( )
            self.f = self.fB
            self.g = self.gB
            self.x = self.xB        
        if self.maximize:
            self.f *= -1.0
            self.g *= -1.0
        return super ( MoreThuenteLineSearcherWithMaxState, self ).Finalize ( )

    def SetUp ( self ):
        """Set up the state."""
        # . Get the bounds on the function.
        if hasattr ( self.objectiveFunction, "GetBounds" ):
            ( self.lowerBound, self.upperBound ) = self.objectiveFunction.GetBounds ( )
        else:
            self.lowerBound = 0.0
            self.upperBound = _LargeNumber
        # . Get the values at the origin.
        ( self.x0, self.f0, self.g0 ) = self.objectiveFunction.GetValuesAtOrigin ( )
        self.objectiveFunction.StoreCurrentAsBest ( ) # . Assume this is true here.
        if self.maximize:
           self.f0 *= -1.0
           self.g0 *= -1.0
        # . Values at other points.
        self.xB = self.x0 ; self.fB = self.f0 ; self.gB = self.g0
        self.xE = self.x0 ; self.fE = self.f0 ; self.gE = self.g0
        # . Other variables.
        self.isBracketed         = False
        self.useModifiedFunction = True
        self.width               = ( self.upperBound - self.lowerBound )
        self.width2              = 2.0 * self.width

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MoreThuenteLineSearcherWithMax ( MoreThuenteLineSearcher ):
    """More-Thuente line searcher with maximize option."""

    _attributable = dict ( MoreThuenteLineSearcher._attributable )
    _stateObject  = MoreThuenteLineSearcherWithMaxState
    _summarizable = dict ( MoreThuenteLineSearcher._summarizable )
    _attributable.update ( { "maximize"             : False                    ,
                             "rmsdMaximumTolerance" : 1.5                      } )
    _summarizable.update ( { "maximize"             : "Maximize"               ,
                             "rmsdMaximumTolerance" : "RMSD Maximum Tolerance" } )

    def Continue ( self, state ):
        """Check to see if iterations should continue."""
        #if (self.maximize):
        #    state.isConverged = ( state.f <= state.fTest ) and ( fabs ( state.g ) <= ( 0.1) )
        #else:
        state.isConverged = ( state.f <= state.fTest ) and ( fabs ( state.g ) <= ( self.gradientTolerance * fabs ( state.g0 ) ) )
        if   state.isConverged:                                  state.statusMessage = "Converged."
        elif state.numberOfIterations >= self.maximumIterations: state.statusMessage = "Too many iterations."
        elif state.error is not None:                            state.statusMessage = "Error: " + state.error
        elif ( state.x == state.lowerBound ) and ( ( state.f >  state.fTest ) or ( state.g >= state.gTest ) ):
            state.statusMessage = "The line search step is too small."
        elif ( state.x == state.upperBound ) and ( ( state.f <= state.fTest ) or ( state.g <= state.gTest ) ):
            state.statusMessage = "The line search step is too large."
        elif state.isBracketed and ( ( state.xMaximum - state.xMinimum ) <= ( self.variableTolerance * state.xMaximum ) ):
             state.statusMessage = "The relative width of the interval of uncertainty is too small."
        elif state.isBracketed and ( ( state.x <= state.xMinimum ) or ( state.x >= state.xMaximum ) ):
             state.statusMessage = "Rounding errors prevent further progress in the line search."
        elif state.maximize and state.x > self.rmsdMaximumTolerance:
             state.statusMessage = "Too high displacement during maximization! Refinement of existing path point failed."
        else: state.statusMessage = None
        return ( state.statusMessage is None )

    def Initialize ( self, state ):
        """Initialization before iteration."""
        # . Get the initial point.
        state.maximize = self.maximize
        state.SetUp()
        state.x = state.x0 + self.initialStep
        # . Various checks.
        if ( state.x < state.lowerBound ) or ( state.x > state.upperBound ):
            raise ValueError ( "Initial value of variable, {:g}, outside range [{:g},{:g}].".format ( state.x, state.lowerBound, state.upperBound ) )
        isOK = ( ( state.x > state.x0 ) and ( state.g0 < 0.0 ) ) or ( ( state.x < state.x0 ) and ( state.g0 > 0.0 ) )
        # . Function and gradients at the new point.
        super ( MoreThuenteLineSearcher, self ).Initialize ( state )
        if (state.maximize):
           state.f *= -1.0
           state.g *= -1.0
        # . Test quantities.
        state.gTest  = self.functionTolerance * state.g0 # . Curvature.
        state.gTestM = min ( self.functionTolerance, self.gradientTolerance ) * state.g0 # . Test for modified function.
        state.fTest  = state.f0 + state.x * state.gTest # . Sufficient descent. 
        # . Range quantities.
        state.xMinimum = state.x0 
        state.xMaximum = state.x + self.upperStepExtrapolationFactor * ( state.x - state.x0 )

    def Iteration ( self, state ):
        """Perform an iteration."""
        super ( MoreThuenteLineSearcherWithMax, self ).Iteration ( state )
        # The iteration will remain the same, only for the maximization following will be changed
        if state.maximize:
            state.f *= -1.0
            state.g *= -1.0

    @staticmethod
    def TrialSteps ( xP, fP, gP, xQ, fQ, gQ, useQuadraticFGF = False ):
        """Calculate trial steps."""
        try:
           ( cIsOK, cubicStep ) = CubicMinimizerFGFG  ( xP, fP, gP, xQ, fQ, gQ )
        except ZeroDivisionError:
           cIsOK = False
           cubicStep = 0.0
           logFile.Paragraph ( "More-Thuente line searcher: cubicStep too small." )
        if useQuadraticFGF: ( qIsOK, quadraticStep ) = QuadraticExtremumFGF ( xP, fP, gP, xQ, fQ     )
        else:               ( qIsOK, quadraticStep ) = QuadraticExtremumGG  ( xP,     gP, xQ,     gQ )
        if cIsOK and qIsOK: cubicStepIsFarther = ( fabs ( cubicStep - xP ) > fabs ( quadraticStep - xP ) )
        else:               cubicStepIsFarther = False
        return ( cubicStep, quadraticStep, cubicStepIsFarther, cIsOK )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
