"""A path cubic spline."""

import math

from  .CubicSpline import CubicSpline
from  .SplineError import SplineError
from ..Arrays      import Array

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Gaussian quadrature parameters - 8 point rule (points, weights).
# . Better to put this in Quadrature subpackage?
_GQData = ( (  0.183434642495650, 0.362683783378362 ) ,
            (  0.525532409916329, 0.313706645877887 ) ,
            (  0.796666477413627, 0.222381034453374 ) ,
            (  0.960289856497536, 0.101228536290376 ) ,
            ( -0.960289856497536, 0.101228536290376 ) ,
            ( -0.796666477413627, 0.222381034453374 ) ,
            ( -0.525532409916329, 0.313706645877887 ) ,
            ( -0.183434642495650, 0.362683783378362 ) )

# . Tolerances.
_CurvatureTolerance = 1.0e-7
_RootTolerance      = 1.0e-4

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PathCubicSpline ( CubicSpline ):
    """A multi-cubic spline representing a path."""

    def _Initialize ( self ):
        """Initialization."""
        super ( PathCubicSpline, self )._Initialize ( )
        self.arcLengths = None

    def ArcLength ( self, x0, x1, g = None ):
        """Calculate the arc length between two abscissa values."""
        # . Find the abscissa starting and stopping values (do by segment).
        segmentData = []
        iLower      = self.Locate ( x0 )
        iUpper      = self.Locate ( x1 )
        for iSegment in range ( iLower, iUpper + 1 ):
            if iSegment == iLower: xStart = x0
            else:                  xStart = self.x[iSegment]
            if iSegment == iUpper: xStop  = x1
            else:                  xStop  = self.x[iSegment+1]
            segmentData.append ( ( xStart, xStop ) )
        # . Allocate space.
        if g is None:
            g = Array.WithExtent ( self.numberOfSplines )
            g.Set ( 0.0 )
        # . Do the integration of gNorm2 for each of the segments (using Gaussian quadrature).
        arcLength = 0.0
        for ( xStart, xStop ) in segmentData:
            c = 0.5 * ( xStart + xStop  )
            f = 0.5 * ( xStop  - xStart )
            for ( p, w ) in _GQData:
                arcLength += ( f * w * self.GradientNorm2 ( c + f * p, g = g ) )
        return arcLength

    def Curvature ( self, x, g = None, h = None ):
        """Calculate the curvature."""
        # . Definition?
        if g is None:
            g = Array.WithExtent ( self.numberOfSplines )
            g.Set ( 0.0 )
        if h is None:
            h = Array.WithExtent ( self.numberOfSplines )
            h.Set ( 0.0 )
        self.EvaluateN ( x, g = g, h = h )
        gNorm2 = g.Norm2 ( )
        if gNorm2 > _CurvatureTolerance:
            g.Scale ( 1.0 / gNorm2 )
            ht = h.Dot ( g )
            h.Add ( g, scale = -ht )
            curvature = h.Norm2 ( ) / ( gNorm2**2 )
        else:
            curvature = 0.0
        return curvature

    def Distance ( self, x0, x1 ):
        """Return the distance."""
        # . Where does this come from - like an RMS?
        return self.ArcLength ( x0, x1 ) / math.sqrt ( self.numberOfSplines )

    def DetermineArcLengths ( self ):
        """Calculate the arc lengths for each segment."""
        if self.arcLengths is None:
            self.arcLengths = Array.WithExtent ( self.numberOfPoints - 1 )
            self.arcLengths.Set ( 0.0 )
            for i in range ( self.numberOfPoints - 1 ):
                self.arcLengths[i] = self.ArcLength ( self.x[i], self.x[i+1] )

    def FindAbscissaFromArcLength ( self, arcLength ):
        """Find the value of the spline abscissa given an arc length."""
        # . Make sure the arc lengths exist.
        self.DetermineArcLengths ( )
        # . Check for bad values.
        if ( arcLength < 0.0 ) or ( arcLength > sum ( self.arcLengths ) ): raise SplineError ( "Invalid arc length value." )
        # . Determine the segment within which the arc length falls.
        iSegment = -1
        lower    = 0.0
        upper    = 0.0
        for i in range ( self.numberOfPoints - 1 ):
            lower  = upper
            upper += self.arcLengths[i]
            if   arcLength == lower: return self.x[i]
            elif arcLength == upper: return self.x[i+1]
            elif ( arcLength > lower ) and ( arcLength < upper ):
                iSegment = i
                break
        # . Search for a root with the secant method.
        # . Initialization (the function is arcLength - required arcLength).
        arcLength -= lower
        upper     -= lower
        xLower     = self.x[iSegment]
        xUpper     = self.x[iSegment+1]
        a = xLower ; fa = - arcLength
        b = xUpper ; fb = upper - arcLength
        root = b
        # . Find the root.
        while math.fabs ( fb / upper ) > _RootTolerance:
            root   = b -  fb * ( b - a ) / ( fb - fa )
            length = self.ArcLength ( xLower, root )
            a = b    ; fa = fb
            b = root ; fb = length - arcLength
        return root

    def GradientNorm2 ( self, x, g = None ):
        """Calculate the gradient 2-norm."""
        if g is None:
            g = Array.WithExtent ( self.numberOfSplines )
            g.Set ( 0.0 )
        self.EvaluateN ( x, g = g )
        return g.Norm2 ( )

    def Normal ( self, x, normal = None ):
        """Calculate the normal."""
        # . Second Frenet vector.
        g = Array.WithExtent ( self.numberOfSplines )
        g.Set ( 0.0 )
        if normal is None:
            h = Array.WithExtent ( self.numberOfSplines )
            h.Set ( 0.0 )
        else:
            h = normal
        self.EvaluateN ( x, g = g, h = h )
        g.Normalize ( )
        ht = h.Dot ( g )
        h.Add ( g, scale = -ht )
        h.Normalize ( )
        return h

    def Tangent ( self, x, tangent = None ):
        """Calculate the tangent."""
        # . First Frenet vector.
        if tangent is None:
            tangent = Array.WithExtent ( self.numberOfSplines )
            tangent.Set ( 0.0 )
        self.EvaluateN ( x, g = tangent )
        return tangent.Normalize ( )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
