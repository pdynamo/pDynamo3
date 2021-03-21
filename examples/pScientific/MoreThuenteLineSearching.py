"""Test for the More-Thuente line search algorithm.

The tests are the same as in the original paper."""

from math                                   import cos, sin, pi, sqrt
from pCore                                  import logFile                         , \
                                                   LogFileActive
from pScientific.ObjectiveFunctionIterators import MoreThuenteLineSearcher         , \
                                                   UniDimensionalObjectiveFunction

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . A large number.
_LargeNumber = 1.0e+20

#===================================================================================================================================
# . Test functions.
#===================================================================================================================================
class TestFunction ( UniDimensionalObjectiveFunction ):
    """Base class for test functions."""

    _attributable = dict ( UniDimensionalObjectiveFunction._attributable )

    def FunctionGradient ( self, x ):
        """Evaluate the function and gradient."""
        return ( self.Function ( x ), self.Gradient ( x ) )

    def GetBounds ( self ): return ( 0.0, _LargeNumber )

    def GetValuesAtOrigin ( self ): return ( 0.0, self.Function ( 0.0 ), self.Gradient ( 0.0 ) )

class TestFunction1 ( TestFunction ):
    """Function 1."""

    _attributable = dict ( TestFunction._attributable )
    _attributable.update ( { "beta" : 2.0 } )

    def Function ( self, x ):
        return - x / ( x**2 + self.beta )

    def Gradient ( self, x ):
        a = 2.0 * ( x**2 ) / ( ( x**2 + self.beta )**2 )
        b = 1.0 / ( x**2 + self.beta )
        return ( a - b )

class TestFunction2 ( TestFunction ):
    """Function 2."""


    _attributable = dict ( TestFunction._attributable )
    _attributable.update ( { "beta" : 0.004 } )

    def Function ( self, x ):
        return ( x + self.beta )**5 - 2.0 * ( ( x + self.beta )**4 )

    def Gradient ( self, x ):
        return 5.0 * ( ( x + self.beta )**4 ) - 8.0 * ( ( x + self.beta )**3 )

class TestFunction3 ( TestFunction ):
    """Function 3."""

    _attributable = dict ( TestFunction._attributable )
    _attributable.update ( { "beta" : 0.01 ,
                             "l"    : 39.0 } )

    def Function ( self, x ):
        if x <= 1.0 - self.beta:
            f = 1.0 - x
        elif x >= 1.0 + self.beta:
            f = x - 1.0
        else:
            f = 0.5 / self.beta * ( x - 1.0 )**2 + 0.5 * self.beta
        return ( f + 2.0 * ( 1.0 - self.beta ) / ( self.l * pi ) * sin ( 0.5 * self.l * pi * x ) )

    def Gradient ( self, x ):
        if x <= 1.0 - self.beta:
            g = -1.0
        elif x >= 1.0 + self.beta:
            g = 1.0
        else:
            g = ( x - 1.0 ) / self.beta
        return g + ( 1.0 - self.beta ) * cos ( 0.5 * self.l * pi * x )

class TestFunction456 ( TestFunction ):
    """Function 456."""

    _attributable = dict ( TestFunction._attributable )
    _attributable.update ( { "beta1" : 0.001 ,
                             "beta2" : 0.001 } )

    def Function ( self, x ):
        g1 = sqrt ( 1.0 + self.beta1**2 ) - self.beta1
        g2 = sqrt ( 1.0 + self.beta2**2 ) - self.beta2
        return ( g1 * sqrt ( ( 1.0 - x )**2 + self.beta2**2 ) + g2 * sqrt ( x**2 + self.beta1**2 ) )

    def Gradient ( self, x ):
        g1 = sqrt ( 1.0 + self.beta1**2 ) - self.beta1
        g2 = sqrt ( 1.0 + self.beta2**2 ) - self.beta2
        return ( - g1 * ( 1.0 - x ) / sqrt ( ( 1.0 - x )**2 + self.beta2**2 ) + g2 * x / sqrt ( x**2 + self.beta1**2 ) )

#===================================================================================================================================
# . Other functions.
#===================================================================================================================================
def ReportsSummary ( reports, log = logFile ):
    """Write out a summary of the reports."""
    if LogFileActive ( log ):
        table = log.GetTable ( columns = [ 10, 10, 10, 10, 10, 10, 10, 10 ] )
        table.Start ( )
        table.Title   ( "Minimization Results" )
        table.Heading ( "Function"  )
        table.Heading ( "Step"      )
        table.Heading ( "Info"      )
        table.Heading ( "Calls"     )
        table.Heading ( "Variable"  )
        table.Heading ( "Gradient"  )
        table.Heading ( "Function"  )
        table.Heading ( "Converged" )
        for report in reports:
            table.Entry ( "{:d}"  .format ( report["Function Index"] ) )
            table.Entry ( "{:.1g}".format ( report["Initial Step"  ] ) )
            message = report["Status Message"]
            if   message.startswith ( "Minimization converged" ): info = 1
            elif message.startswith ( "Rounding errors "       ): info = 2
            elif message.startswith ( "The relative width"     ): info = 3
            else: raise ValueError ( "Unknown message: " + message )
            table.Entry ( "{:d}"  .format ( info ) )
            table.Entry ( "{:d}"  .format ( report["Function Calls"] ) )
            table.Entry ( "{:.3f}".format ( report["Variable"      ] ) )
            table.Entry ( "{:.3g}".format ( report["Gradient"      ] ) )
            table.Entry ( "{:.3g}".format ( report["Function Value"] ) )
            table.Entry ( "{!r}"  .format ( report["Converged"     ] ) )
        table.Stop ( )

#===================================================================================================================================
# . Test function data.
#===================================================================================================================================
# . Function data.
_FunctionData = ( ( TestFunction1  , {}, 1.0e-3, 1.0e-1 ) ,
                  ( TestFunction2  , {}, 1.0e-1, 1.0e-1 ) ,
                  ( TestFunction3  , {}, 1.0e-1, 1.0e-1 ) ,
                  ( TestFunction456, { "beta1" : 1.0e-3, "beta2" : 1.0e-3 }, 1.0e-3, 1.0e-3 ) ,
                  ( TestFunction456, { "beta1" : 1.0e-2, "beta2" : 1.0e-3 }, 1.0e-3, 1.0e-3 ) ,
                  ( TestFunction456, { "beta1" : 1.0e-3, "beta2" : 1.0e-2 }, 1.0e-3, 1.0e-3 ) )

# . Initial steps.
_InitialSteps = ( 1.0e-3, 1.0e-1, 1.0e+1, 1.0e+3 )

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Initialization.
isOK = True

# . Run the tests.
reports            = []
numberNotConverged =  0
for ( i, ( functionClass, attributes, functionTolerance, gradientTolerance ) ) in enumerate ( _FunctionData ):
    of = functionClass ( )
    of.__dict__.update ( attributes )
    for initialStep in _InitialSteps:
        mt = MoreThuenteLineSearcher.WithOptions ( functionTolerance            = functionTolerance ,
                                                   gradientTolerance            = gradientTolerance ,
                                                   initialStep                  = initialStep       ,
                                                   lowerStepExtrapolationFactor = 1.1               ,
                                                   maximumIterations            = 500               ,
                                                   variableTolerance            = 1.0e-15           )
        mt.Summary ( )
        report = mt.Iterate ( of )
        report["Function Index"] = i+1
        report["Initial Step"  ] = initialStep
        reports.append ( report )
        if not report.get ( "Converged", False ): numberNotConverged += 1

# . Footer.
ReportsSummary ( reports )
logFile.Footer ( )
if ( numberNotConverged > 0 ): TestScriptExit_Fail ( )
