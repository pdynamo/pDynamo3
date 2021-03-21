"""ABFS integrator class."""

import math

from pCore                  import SummarizableObject
from pScientific.Arrays     import Array
from pScientific.Quadrature import AdaptiveSimpsonsRule

# . Default cutoffs suitable for Angstroms.
# . Scale if want atomic units.

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class ABFSIntegrator ( SummarizableObject ):
    """ABFS integrator."""

    _attributable = dict ( SummarizableObject._attributable )
    _classLabel   = "ABFS Integrator"
    _summarizable = dict ( SummarizableObject._summarizable )
    _attributable.update ( { "dampingCutOff"        :    0.2                  ,   # . Smaller value - good for QC/MM integrals.
                             "innerCutOff"          :    8.0                  , 
                             "outerCutOff"          :   12.0                  , 
                             "integrationTolerance" :    1.0e-15              ,   # . Minimal reasonable value.
                             "maximumIterations"    : 1000                    ,
                             "pointDensity"         :   50                    ,
                             "_x"                   : None                    } ) # . Needs to be removed if options changed!
    _summarizable.update ( { "dampingCutOff"        : "Damping Cutoff"        ,
                             "innerCutOff"          : "Inner Cutoff"          ,
                             "outerCutOff"          : "Outer Cutoff"          ,
                             "integrationTolerance" : "Integration Tolerance" ,
                             "maximumIterations"    : "Maximum Iterations"    ,
                             "pointDensity"         : "Point Density"         } )

    def _Switch ( self ):
        """Default switch function."""
        def S ( r ):
            r2    = r**2
            rOff2 = self.outerCutOff**2
            rOn2  = self.innerCutOff**2
            return ( ( rOff2 - r2 )**2 * ( rOff2 + 2.0 * r2 - 3.0 * rOn2 ) / ( rOff2 - rOn2 )**3 )
        return S

    def _XValues ( self ):
        """The x-values."""
        n = int ( math.ceil ( float ( self.pointDensity ) * self.outerCutOff + 1.0 ) )
        x = Array.WithExtent ( n )
        d = self.outerCutOff / float ( n - 1 )
        x[ 0] = 0.0
        for i in range ( 1, n - 1 ): x[i] = d * float ( i )
        x[-1] = self.outerCutOff
        return x

    def Integrate ( self, F, G ):
        """Integration."""
        # . Function definition.
        S = self._Switch ( )
        def GS ( r ): return ( G ( r ) * S ( r ) )
        # . Initialization.
        n = len ( self.x )
        y = Array.WithExtent ( n )
        # . Calculate y.
        # . Switching region - rOn to rOff.
        integral = 0.0
        iUpper   = n - 2
        y[-1]    = 0.0
        for i in range ( n - 2, -1, -1 ):
            rLower = self.x[i  ]
            rUpper = self.x[i+1]
            if rLower < self.innerCutOff:
                iUpper = i
                if rUpper > self.innerCutOff:
                    ( value, state ) = AdaptiveSimpsonsRule ( GS, self.innerCutOff, rUpper, self.integrationTolerance, self.maximumIterations )
                    integral -= value
                break
            ( value, state ) = AdaptiveSimpsonsRule ( GS, rLower, rUpper, self.integrationTolerance, self.maximumIterations )
            integral -= value
            y[i]      = integral
        # . Normal and damping regions - 0 to rOn.
        # . Determine some constants.
        offSet = integral - F ( self.innerCutOff )
        if self.dampingCutOff > 0.0:
            fDamp = F ( self.dampingCutOff ) + offSet
            gDamp = G ( self.dampingCutOff )
            alpha =   0.5 * gDamp / self.dampingCutOff
            beta  = - 0.5 * gDamp * self.dampingCutOff + fDamp
        else:
            alpha = 0.0
            beta  = 0.0
        # . Remaining points.
        for i in range ( iUpper + 1 ):
            r = self.x[i]
            if r < self.dampingCutOff: y[i] = alpha * r**2 + beta
            else:                      y[i] = F ( r ) + offSet
        # . Finish up.
        return y

    @property
    def x ( self ):
        """The x-values."""
        if self._x is None: self._x = self._XValues ( )
        return self._x

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

