"""More-Thuente line searcher."""

# . After the method of J. J. More and D. J. Thuente (1992).

from  math                           import fabs
from .Interpolation                  import CubicMinimizerFGFG             , \
                                            QuadraticExtremumFGF           , \
                                            QuadraticExtremumGG
from .ObjectiveFunctionIteratorError import ObjectiveFunctionIteratorError
from .UniDimensionalMinimizer        import UniDimensionalMinimizer        , \
                                            UniDimensionalMinimizerState

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . A large number for the upper bound of the line search.
_LargeNumber = 1.0e+20

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class MoreThuenteLineSearcherState ( UniDimensionalMinimizerState ):
    """More-Thuente line searcher state."""

    _attributable = dict ( UniDimensionalMinimizerState._attributable )
    _attributable.update ( { "f0"                  : 0.0   ,
                             "fB"                  : 0.0   ,
                             "fE"                  : 0.0   ,
                             "fTest"               : 0.0   ,
                             "g0"                  : 0.0   ,
                             "gB"                  : 0.0   ,
                             "gE"                  : 0.0   ,
                             "gTest"               : 0.0   ,
                             "gTestM"              : 0.0   ,
                             "isBracketed"         : False ,
                             "lowerBound"          : None  ,
                             "upperBound"          : None  ,
                             "useModifiedFunction" : True  ,
                             "width"               : 0.0   ,
                             "width2"              : 0.0   ,
                             "x0"                  : 0.0   ,
                             "xB"                  : 0.0   ,
                             "xE"                  : 0.0   ,
                             "xMaximum"            : 0.0   ,
                             "xMinimum"            : 0.0   } )

    def Finalize ( self ):
        """Finalization."""
        # . Return the best point if abnormal termination.
        if ( not self.isConverged ) and ( self.error is None ) and ( self.f > self.fB ):
            self.objectiveFunction.RestoreBestToCurrent ( )
            self.f = self.fB
            self.g = self.gB
            self.x = self.xB
        return super ( MoreThuenteLineSearcherState, self ).Finalize ( )

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
class MoreThuenteLineSearcher ( UniDimensionalMinimizer ):
    """More-Thuente line searcher."""

    _attributable = dict ( UniDimensionalMinimizer._attributable )
    _classLabel   = "More-Thuente Line Searcher"
    _stateObject  = MoreThuenteLineSearcherState
    _summarizable = dict ( UniDimensionalMinimizer._summarizable )
    _attributable.update ( { "bisectionFactor"              : 0.66                              ,
                             "functionTolerance"            : 1.0e-03                           ,
                             "gradientTolerance"            : 0.9                               , # . Reset this value from superclass default.
                             "initialStep"                  : 1.0                               ,
                             "lowerStepExtrapolationFactor" : 1.1                               , # . Or -1.0?
                             "safeguardFactor"              : 0.66                              ,
                             "upperStepExtrapolationFactor" : 4.0                               ,
                             "variableTolerance"            : 1.0e-15                           } )
    _summarizable.update ( { "bisectionFactor"              : "Bisection Factor"                ,
                             "functionTolerance"            : "Function Tolerance"              ,
                             "lowerStepExtrapolationFactor" : "Lower Step Extrapolation Factor" ,
                             "safeguardFactor"              : "Safeguard Factor"                ,
                             "upperStepExtrapolationFactor" : "Upper Step Extrapolation Factor" ,
                             "variableTolerance"            : "Variable Tolerance"              } )

    def Continue ( self, state ):
        """Check to see if iterations should continue."""
        state.isConverged = ( state.error      is None        ) and \
                            ( state.f          <= state.fTest ) and \
                            ( fabs ( state.g ) <= ( self.gradientTolerance * fabs ( state.g0 ) ) )
        if   state.isConverged:                                  state.statusMessage = "Minimization converged."
        elif state.numberOfIterations >= self.maximumIterations: state.statusMessage = "Too many iterations."
        elif state.error is not None:                            state.statusMessage = "Minimization error: " + state.error
        elif ( state.x == state.lowerBound ) and ( ( state.f >  state.fTest ) or ( state.g >= state.gTest ) ):
            state.statusMessage = "The line search step is too small."
        elif ( state.x == state.upperBound ) and ( ( state.f <= state.fTest ) or ( state.g <= state.gTest ) ):
            state.statusMessage = "The line search step is too large."
        elif state.isBracketed and ( ( state.xMaximum - state.xMinimum ) <= ( self.variableTolerance * state.xMaximum ) ):
             state.statusMessage = "The relative width of the interval of uncertainty is too small."
        elif state.isBracketed and ( ( state.x <= state.xMinimum ) or ( state.x >= state.xMaximum ) ):
             state.statusMessage = "Rounding errors prevent further progress in the line search."
        else: state.statusMessage = None
        return ( state.statusMessage is None )

    def Initialize ( self, state ):
        """Initialization before iteration."""
        # . Get the initial point.
        state.x = state.x0 + self.initialStep
        # . Various checks.
        if ( state.x < state.lowerBound ) or ( state.x > state.upperBound ):
            raise ObjectiveFunctionIteratorError ( "Initial value of variable, {:g}, outside range [{:g},{:g}].".format ( state.x, state.lowerBound, state.upperBound ) )
        isOK = ( ( state.x > state.x0 ) and ( state.g0 < 0.0 ) ) or ( ( state.x < state.x0 ) and ( state.g0 > 0.0 ) )
        if not isOK:
            raise ObjectiveFunctionIteratorError ( "The gradient at the origin does not descend in the direction of the search." )
        # . Function and gradients at the new point.
        super ( MoreThuenteLineSearcher, self ).Initialize ( state )
        # . Test quantities.
        state.gTest  = self.functionTolerance * state.g0 # . Curvature.
        state.gTestM = min ( self.functionTolerance, self.gradientTolerance ) * state.g0 # . Test for modified function.
        state.fTest  = state.f0 + state.x * state.gTest # . Sufficient descent. 
        # . Range quantities.
        state.xMinimum = state.x0
        state.xMaximum = state.x + self.upperStepExtrapolationFactor * ( state.x - state.x0 )

    def Iteration ( self, state ):
        """Perform an iteration."""
        try:
            # . Check whether to use the modified or unmodified function.
            if state.useModifiedFunction and ( state.f < state.fTest ) and ( state.g >= state.gTestM ): state.useModifiedFunction = False
            # . Use the modified function.
            if state.useModifiedFunction and ( state.f < state.fB ) and ( state.f > state.fTest ):
                # . Define the modified function.
                fM  = state.f  - state.gTest * state.x
                fBM = state.fB - state.gTest * state.xB
                fEM = state.fE - state.gTest * state.xE
                gM  = state.g  - state.gTest
                gBM = state.gB - state.gTest
                gEM = state.gE - state.gTest
                # . Compute the new step.
                ( state.stepType, state.x, state.isBracketed, state.xB, fBM, gBM, state.xE, fEM, gEM ) = self.NewStep ( state.xB , fBM , gBM ,
                                                                                                                        state.xE , fEM , gEM ,
                                                                                                                        state.x  , fM  , gM  ,
                                                                                                                        state.isBracketed    ,
                                                                                                                        state.xMinimum       ,
                                                                                                                        state.xMaximum       )
                # . Reset the function.
                state.fB = fBM + state.gTest * state.xB
                state.fE = fEM + state.gTest * state.xE
                state.gB = gBM + state.gTest
                state.gE = gEM + state.gTest
            # . Use the unmodified function.
            else:
                ( state.stepType, state.x, state.isBracketed, state.xB, state.fB, state.gB, state.xE, state.fE, state.gE ) = self.NewStep ( state.xB , state.fB , state.gB ,
                                                                                                                                            state.xE , state.fE , state.gE ,
                                                                                                                                            state.x  , state.f  , state.g  ,
                                                                                                                                            state.isBracketed              ,
                                                                                                                                            state.xMinimum                 ,
                                                                                                                                            state.xMaximum                 )
            # . Ensure that the best point is stored if it is has changed.
            if state.stepType[-1] != "A": state.objectiveFunction.StoreCurrentAsBest ( )
            # . Perform a bisection step.
            if state.isBracketed:
                size = fabs ( state.xE - state.xB )
                if size > ( self.bisectionFactor * state.width2 ): state.x = 0.5 * ( state.xB + state.xE )
                state.width2 = state.width
                state.width  = size
            # . Redetermine the interval of uncertainty.
            if state.isBracketed:
                state.xMinimum = min ( state.xB, state.xE )
                state.xMaximum = max ( state.xB, state.xE )
            else:
                state.xMinimum = state.x + self.lowerStepExtrapolationFactor * ( state.x - state.xB )
                state.xMaximum = state.x + self.upperStepExtrapolationFactor * ( state.x - state.xB )
            # . Ensure that the step is within range.
            state.x = max ( state.x, state.lowerBound )
            state.x = min ( state.x, state.upperBound )
            # . Determine the function and gradient.
            self.FunctionGradient ( state )
            state.fTest = state.f0 + state.x * state.gTest
        except Exception as error:
            state.HandleError ( error )
        state.numberOfIterations += 1

    def NewStep ( self, xB, fB, gB, xE, fE, gE, x, f, g, isBracketed, xMinimum, xMaximum ):
        """Determine a safeguarded step and update the interval of uncertainty."""
        # . Estimate the new step.
        # . Case 1: trial point has a higher function value.
        if f > fB:
            isBracketed = True
            stepType    = "1"
            ( cubicStep, quadraticStep, cubicStepIsFarther, isOK ) = self.TrialSteps ( xB, fB, gB, x, f, g, useQuadraticFGF = True )
            if cubicStepIsFarther: xNew = 0.5 * ( cubicStep + quadraticStep )
            else:                  xNew =         cubicStep
        # . Case 2: trial point has a lower function value and a derivative of opposite sign.
        elif ( ( gB / fabs ( gB ) ) * g ) < 0.0:
            isBracketed = True
            stepType    = "2"
            ( cubicStep, quadraticStep, cubicStepIsFarther, isOK ) = self.TrialSteps ( x, f, g, xB, fB, gB )
            if cubicStepIsFarther: xNew = cubicStep
            else:                  xNew = quadraticStep
        # . Case 3: trial point has a lower function value and a smaller derivative of the same sign.
        elif fabs ( g ) < fabs ( gB ):
            stepType = "3"
            ( cubicStep, quadraticStep, cubicStepIsFarther, isOK ) = self.TrialSteps ( x, f, g, xB, fB, gB )
            if ( not isOK ) or ( ( cubicStep < x ) and ( x > xB ) ) or ( ( cubicStep > x ) and ( x < xB ) ):
                if x > xB: cubicStep = xMaximum
                else:      cubicStep = xMinimum
            if isBracketed:
                if cubicStepIsFarther: xNew = quadraticStep
                else:                  xNew = cubicStep
                # . It is possible that safeguarding should be moved until after the interval
                # . of uncertainty has been updated to ensure that the step is within range.
                xTemp = x + self.safeguardFactor * ( xE - x )
                if x > xB: xNew = min ( xTemp, xNew )
                else:      xNew = max ( xTemp, xNew )
            else:
                if cubicStepIsFarther: xNew = cubicStep
                else:                  xNew = quadraticStep
                if   xNew > xMaximum: xNew = xMaximum
                elif xNew < xMinimum: xNew = xMinimum
        # . Case 4: trial point has a lower function value and a larger derivative of the same sign.
        else:
            stepType = "4"
            if isBracketed: ( isOK, xNew ) = CubicMinimizerFGFG ( x, f, g, xE, fE, gE )
            elif x > xB: xNew = xMaximum
            else:        xNew = xMinimum
        # . Update the interval of uncertainty.
        # . Case A.
        if f > fB:
            xE = x
            fE = f
            gE = g
            stepType += "A"
        else:
            # . Case C.
            if ( g * ( xB - x ) ) < 0.0:
                xE = xB
                fE = fB
                gE = gB
                stepType += "C"
            # . Case B.
            else:
                stepType += "B"
            # . Cases B and C.
            xB = x
            fB = f
            gB = g
        # . Finish up.
        return ( stepType, xNew, isBracketed, xB, fB, gB, xE, fE, gE )

    @staticmethod
    def TrialSteps ( xP, fP, gP, xQ, fQ, gQ, useQuadraticFGF = False ):
        """Calculate trial steps."""
        ( cIsOK, cubicStep ) = CubicMinimizerFGFG  ( xP, fP, gP, xQ, fQ, gQ )
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
