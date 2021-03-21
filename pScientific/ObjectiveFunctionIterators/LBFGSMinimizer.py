"""A very basic L-BFGS minimizer."""

from  .MultiDimensionalMinimizer import MultiDimensionalMinimizer      , \
                                        MultiDimensionalMinimizerState
from ..Arrays                    import Array

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class LBFGSMinimizerState ( MultiDimensionalMinimizerState ):
    """LBFGS minimizer state."""

    _attributable = dict ( MultiDimensionalMinimizerState._attributable )
    _attributable.update ( { "aux"      : None ,
                             "deltaG"   : None ,
                             "deltaX"   : None ,
                             "gOld"     : None ,
                             "history"  :    0 ,
                             "iHistory" :    0 ,
                             "rho"      : None ,
                             "xOld"     : None } )

    @classmethod
    def FromObjectiveFunction ( selfClass, objectiveFunction, history = 0 ):
        """Constructor given an objective function."""
        self         = super ( MultiDimensionalMinimizerState, selfClass ).FromObjectiveFunction ( objectiveFunction )
        self.history = history
        self.SetUp ( )
        return self

    @classmethod
    def FromVariableArray ( selfClass, x, d = None, g = None, history = 0, surrogateObjectiveFunction = None ):
        """Constructor given a variable array."""
        self = selfClass ( )
        self.history           = history
        self.numberOfVariables = len ( x )
        self.d = d
        self.g = g
        self.x = x
        self.SetUp ( )
        self.ExtractSurrogateData ( surrogateObjectiveFunction )
        return self

    def SetUp ( self ):
       """Set up the state."""
       super ( LBFGSMinimizerState, self ).SetUp ( )
       self.aux    =   Array.WithExtent ( self.history           ) ; self.aux.Set  ( 0.0 )
       self.gOld   =   Array.WithExtent ( self.numberOfVariables ) ; self.gOld.Set ( 0.0 )
       self.xOld   =   Array.WithExtent ( self.numberOfVariables ) ; self.xOld.Set ( 0.0 )
       self.deltaG = [ Array.WithExtent ( self.numberOfVariables ) for i in range ( self.history ) ]
       self.deltaX = [ Array.WithExtent ( self.numberOfVariables ) for i in range ( self.history ) ]
       self.rho    = [ 0.0 for i in range ( self.history ) ]
       for ( g, x ) in zip ( self.deltaG, self.deltaX ):
           g.Set ( 0.0 )
           x.Set ( 0.0 )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class LBFGSMinimizer ( MultiDimensionalMinimizer ):
    """L-BFGS minimization."""

    _attributable = dict ( MultiDimensionalMinimizer._attributable )
    _classLabel   = "L-BFGS Minimizer"
    _stateObject  = LBFGSMinimizerState
    _summarizable = dict ( MultiDimensionalMinimizer._summarizable )
    _attributable.update ( { "history"     : 10             ,
                             "maximumStep" : 0.01           } )
    _summarizable.update ( { "history"     : "History"      ,
                             "maximumStep" : "Maximum Step" } )

    def DetermineStep ( self, state ):
        """Determine the step."""
        # . Reshuffle the data.
        if state.iHistory > self.history:
            temp = state.deltaG.pop ( 0 ) ; state.deltaG.append ( temp )
            temp = state.deltaX.pop ( 0 ) ; state.deltaX.append ( temp )
            temp = state.rho.pop    ( 0 ) ; state.rho.append    ( temp )
        # . Calculate new correction data.
        if state.iHistory > 0:
            i = min ( state.iHistory, self.history ) - 1
            state.g.CopyTo ( state.deltaG[i] ) ; state.deltaG[i].Add ( state.gOld, scale = -1.0 )
            state.x.CopyTo ( state.deltaX[i] ) ; state.deltaX[i].Add ( state.xOld, scale = -1.0 )
            hGG = state.deltaG[i].DotSelf ( )
            hGX = state.deltaG[i].Dot ( state.deltaX[i] )
            state.rho[i] = 1.0 / hGX
            hScale       = hGX / hGG
        # . Save the old data.
        state.g.CopyTo ( state.gOld )
        state.x.CopyTo ( state.xOld )
        # . Calculate the step.
        state.g.CopyTo ( state.d )
        state.d.Scale ( -1.0 )
        if state.iHistory == 0:
            state.d.Scale ( 1.0 / state.g.DotSelf ( ) )
        else:
            for i in reversed ( range ( min ( state.iHistory, state.history ) ) ):
                a = state.rho[i] * state.d.Dot ( state.deltaX[i] )
                state.d.Add ( state.deltaG[i], scale = -a )
                state.aux[i] = a
            state.d.Scale ( hScale )
            for i in range ( min ( state.iHistory, state.history ) ):
                a = state.aux[i] - state.rho[i] * state.d.Dot ( state.deltaG[i] )
                state.d.Add ( state.deltaX[i], scale = a )
                state.aux[i] = a
        # . Check the step length.
        step = state.d.RootMeanSquare ( )
        if step > self.maximumStep: state.d.Scale ( self.maximumStep / step )
        # . Increment the variables.
        state.x.Add ( state.d )
        # . Finish up.
        state.stepType  = "S{:d}".format ( min ( state.iHistory, self.history ) )
        state.iHistory += 1

    def Initialize ( self, state  ):
        """Initialization before iteration."""
        self.Restart ( state )
        super ( LBFGSMinimizer, self ).Initialize ( state )

    def Restart ( self, state ):
        """Restart."""
        state.iHistory = 0

    def StateFromObjectiveFunction ( self, objectiveFunction ):
        """Set up the state."""
        return self.__class__._stateObject.FromObjectiveFunction ( objectiveFunction, history = self.history )

    def StateFromVariableArray ( self, x, d = None, g = None, surrogateObjectiveFunction = None ):
        """Set up the state."""
        return self.__class__._stateObject.FromVariableArray ( x, d = d, g = g, history = self.history, surrogateObjectiveFunction = surrogateObjectiveFunction )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
