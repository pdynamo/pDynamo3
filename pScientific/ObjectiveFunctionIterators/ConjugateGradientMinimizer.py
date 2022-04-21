"""Conjugate gradient minimizer with a choice of step methods."""

# . The spectral modification follows the method of E. G. Birgin and J. M. Martinez (2001).

import math

from  .MultiDimensionalMinimizer import MultiDimensionalMinimizer      , \
                                        MultiDimensionalMinimizerState
from  .MoreThuenteLineSearch     import MoreThuenteLineSearcher
from  .ObjectiveFunction         import UniDimensionalObjectiveFunction
from ..Arrays                    import Array

#
# . Notes:
#
#   - Additional line search methods.
#   - Methods for determining the initial step for the line search (see NocWri99 p58).
#   - Restart methods (e.g. Pow77 section 4).
#   - Ways of exiting gracefully or recovering if an SCF fails to converge.
#

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Default parameters for the line searcher.
_LineSearcherParameters = { "bisectionFactor"              :  0.66    ,
                            "functionTolerance"            :  1.0e-04 ,
                            "gradientTolerance"            :  0.1     , # . Better for CG methods?
                            "initialStep"                  :  1.0     ,
                            "logFrequency"                 : -1       ,
                            "lowerStepExtrapolationFactor" :  1.1     ,
                            "maximumIterations"            : 20       ,
                            "safeguardFactor"              :  0.66    ,
                            "upperStepExtrapolationFactor" :  4.0     ,
                            "variableTolerance"            :  1.0e-15 }

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class CGObjectiveFunction ( UniDimensionalObjectiveFunction ):
    """The line search objective function."""

    _attributable = dict ( UniDimensionalObjectiveFunction._attributable )
    _attributable.update ( { "d"                 : None ,
                             "f"                 : None ,
                             "f0"                : None ,
                             "fB"                : None ,
                             "g"                 : None ,
                             "g0"                : None ,
                             "gB"                : None ,
                             "objectiveFunction" : None ,
                             "vB"                : None ,
                             "x"                 : None ,
                             "x0"                : None ,
                             "xB"                : None } )

    @classmethod
    def WithOptions ( selfClass, **options ):
        """Constructor from options."""
        self = selfClass ( )
        for ( key, value ) in options.items ( ): setattr ( self, key, value )
        return self

    def FunctionGradient ( self, alpha ):
        """Evaluate the function and gradient."""
        # . Save the explicit variable.
        self.variable = alpha
        # . Set up the implicit variables.
        self.x0.CopyTo ( self.x )
        self.x.Add ( self.d, scale = alpha )
        # . Apply constraints to variables.
        self.objectiveFunction.ApplyConstraintsToVariables ( self.x, referenceVariables = self.x0 )
        # . Function and gradients.
        self.f = self.objectiveFunction.FunctionGradients ( self.x, self.g )
        # . Apply constraints to gradients.
        self.objectiveFunction.ApplyConstraintsToVector ( self.g, referenceVariables = self.x )
        # . Finish up.
        return ( self.f, self.d.Dot ( self.g ) )

#    def GetBounds ( self ): return ( 0.0, _LargeNumber )

    def GetValuesAtOrigin ( self ):
        """Get the values at the origin."""
        return ( 0.0, self.f0, self.g0 )

    def SetValuesAtOrigin ( self, f0, g0 ):
        """Set the values at the origin."""
        self.f0 = f0
        self.g0 = g0

    def RestoreBestToCurrent ( self ):
        """Make the best point the current one."""
        self.f        = self.fB
        self.variable = self.vB
        self.gB.CopyTo ( self.g )
        self.xB.CopyTo ( self.x )

    def StoreCurrentAsBest ( self ):
        """Save the current point as the best point."""
        self.fB = self.f
        self.vB = self.variable
        self.g.CopyTo ( self.gB )
        self.x.CopyTo ( self.xB )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class ConjugateGradientMinimizerState ( MultiDimensionalMinimizerState ):
    """Conjugate gradient minimizer state."""

    _attributable = dict ( MultiDimensionalMinimizerState._attributable )
    _attributable.update ( { "g0"                          : None ,
                             "gDotG"                       : 0.0  ,
                             "lineSearchObjectiveFunction" : None ,
                             "s"                           : None ,
                             "theta"                       : 1.0  ,
                             "x0"                          : None ,
                             "y"                           : None } )

    def SetUp ( self ):
       """Set up the state."""
       # . Arrays.
       super ( ConjugateGradientMinimizerState, self ).SetUp ( )
       self.g0 = Array.WithExtent ( self.numberOfVariables ) ; self.g0.Set ( 0.0 )
       self.s  = Array.WithExtent ( self.numberOfVariables ) ; self.s.Set  ( 0.0 )
       self.x0 = Array.WithExtent ( self.numberOfVariables ) ; self.x0.Set ( 0.0 )
       self.y  = Array.WithExtent ( self.numberOfVariables ) ; self.y.Set  ( 0.0 )
       # . Line-search objective function.
       # . s and y are used here as they are only used for workspace and not for persistent storage.
       self.lineSearchObjectiveFunction = CGObjectiveFunction.WithOptions ( d                 = self.d                 ,
                                                                            g                 = self.g                 ,
                                                                            gB                = self.s                 ,
                                                                            objectiveFunction = self.objectiveFunction ,
                                                                            x                 = self.x                 ,
                                                                            x0                = self.x0                ,
                                                                            xB                = self.y                 )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class ConjugateGradientMinimizer ( MultiDimensionalMinimizer ):
    """Conjugate gradient multidimensional minimizer."""

    # . Possible betaTypes are:
    #   1 and everything else - Perry
    #   2                     - Fletcher-Reeves
    #   3                     - Polak-Ribiere
    #   4                     - Polak-Ribiere Plus

    _attributable = dict ( MultiDimensionalMinimizer._attributable )
    _classLabel   = "Conjugate Gradient Minimizer"
    _stateObject  = ConjugateGradientMinimizerState
    _summarizable = dict ( MultiDimensionalMinimizer._summarizable )
    _attributable.update ( { "betaType"                 : 1                            ,
                             "initialStep"              : 1.0                          ,
                             "lineSearcher"             : None                         ,
                             "maximumTheta"             : 1.0e+10                      ,
                             "minimumTheta"             : 1.0e-10                      ,
                             "steepestDescentTolerance" : 1.0e-03                      ,
                             "useSpectralTheta"         : True                         } )
    _summarizable.update ( { "betaType"                 : "Beta Type"                  ,
                             "initialStep"              : "Initial Step"               ,
                             "lineSearcher"             : "Line Searcher"              ,
                             "maximumTheta"             : "Maximum Theta"              ,
                             "minimumTheta"             : "Minimum Theta"              ,
                             "steepestDescentTolerance" : "Steepest Descent Tolerance" ,
                             "useSpectralTheta"         : "Use Spectral Theta"         } )

    def FunctionGradients ( self, state ):
        """Evaluate the function and its gradients."""
        # . Apply constraints to variables.
        state.objectiveFunction.ApplyConstraintsToVariables ( state.x, referenceVariables = state.x0 ) # . x0 is last iteration.
        # . Function and gradients.
        state.f = state.objectiveFunction.FunctionGradients ( state.x, state.g )
        # . Apply constraints to gradients.
        state.objectiveFunction.ApplyConstraintsToVector ( state.g, referenceVariables = state.x )
        # . Finish up.
        state.rmsGradient = state.g.RootMeanSquare ( )
        state.numberOfFunctionCalls += 1

    def Initialize ( self, state ):
        """Initialization before iteration."""
        super ( ConjugateGradientMinimizer, self ).Initialize ( state )
        state.f0 = state.f
        state.g.CopyTo ( state.g0 )
        state.x.CopyTo ( state.x0 )
        self.InitialSearchDirection ( state )
        parameters = _LineSearcherParameters
        parameters["initialStep"] = self.initialStep
        if self.lineSearcher is None: self.lineSearcher = MoreThuenteLineSearcher.WithOptions ( **parameters )

    def InitialSearchDirection ( self, state ):
        """Determine an initial search direction."""
        state.g.CopyTo ( state.d )
        state.d.Scale  ( -1.0 )
        state.d.Normalize ( )
        state.alpha = self.initialStep
        state.gDotG = state.g.DotSelf ( )

    def Iteration ( self, state ):
        """Perform an iteration."""
        try:
            self.LineSearch         ( state )
            self.NewSearchDirection ( state )
        except Exception as error:
            state.error = error.args[0]
            import traceback, sys
            traceback.print_exc ( file = sys.stdout )
        state.numberOfIterations += 1

    def LineSearch ( self, state ):
        """Do a line search."""
        # . Do the line search and return the best point.
        self.lineSearcher.initialStep = state.alpha
        state.lineSearchObjectiveFunction.SetValuesAtOrigin ( state.f0, state.d.Dot ( state.g0 ) )
        report = self.lineSearcher.Iterate ( state.lineSearchObjectiveFunction, log = None )
        # . Get the results.
        state.alpha                  = report["Variable"      ]
        state.f                      = report["Function Value"]
        state.numberOfFunctionCalls += report.get ( "Function Calls", 0 )
        state.rmsGradient            = state.g.RootMeanSquare ( )
        # . Determine the line search step type (success, partial, error).
        if report["Converged"]:
            statusFlag  = "s"
        elif state.f < state.f0:
            statusFlag  = "p"
        else:
            statusFlag  = "e"
            state.error = "Line search error: " + report["Status Message"] 
        state.stepType = "L{:d}{:s}".format ( report.get ( "Iterations", 0 ), statusFlag )

    def NewSearchDirection ( self, state ):
        """Determine a new search direction."""
        # . Save data.
        state.f0 = state.f
        state.g.CopyTo ( state.y  ) ; state.y.Add ( state.g0, scale = -1.0 )        
        state.x.CopyTo ( state.s  ) ; state.s.Add ( state.x0, scale = -1.0 )        
        state.g.CopyTo ( state.g0 )
        state.x.CopyTo ( state.x0 )
        # . Dot factors.
        gDotG0      = state.gDotG
        state.gDotG = state.g.DotSelf ( )
        gNorm2      = math.sqrt ( state.gDotG )
        sDotS       = state.s.DotSelf ( )
        sDotY       = state.s.Dot ( state.y )
        # . Calculate the theta part of the direction.
        # . Theta is always one if the spectral method is not being used.
        theta0 = state.theta
        if self.useSpectralTheta:
            if sDotY <= 0.0: state.theta = self.maximumTheta
            else:            state.theta = min ( self.maximumTheta, max ( self.minimumTheta, sDotS / sDotY ) )
        state.g.CopyTo (  state.d     )
        state.d.Scale  ( -state.theta )
        # . Calculate the beta part of the direction.
        beta = 0.0
        # . Fletcher-Reeves.
        if self.betaType == 2:
            a = state.theta * state.gDotG
            b = state.alpha * theta0 * gDotG0
        # . Polak-Ribiere (normal and plus).
        elif ( self.betaType == 3 ) or ( self.betaType == 4 ):
            a = state.theta * state.g.Dot ( state.y )
            b = state.alpha * theta0 * gDotG0
            if self.betaType == 4: b = max ( b, 0.0 )
        # . Perry.
        # . y is used as workspace as it is no longer needed (although s is).
        else:
            state.y.Scale ( state.theta )
            state.y.Add ( state.s, scale = -1.0 )
            a = state.y.Dot ( state.g )
            b = sDotY
        # . Do nothing if beta is undefined (i.e. use the steepest descent direction).
        if b != 0.0:
            beta = a / b
            state.d.Add ( state.s, scale = beta )
        # . Choose a steepest descent direction if the composite step is too uphill.
        dNorm2 = state.d.Norm2 ( )
        if state.g.Dot ( state.d ) > ( - self.steepestDescentTolerance * dNorm2 * gNorm2 ):
            state.g.CopyTo ( state.d )
            state.d.Scale  ( -1.0 )
        state.objectiveFunction.ApplyConstraintsToVector ( state.d, referenceVariables = state.x0 ) # . x0 is last iteration.
        state.d.Normalize ( )
        # . Alpha keeps its last value at the moment.
        # . Have an option to reset?

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
