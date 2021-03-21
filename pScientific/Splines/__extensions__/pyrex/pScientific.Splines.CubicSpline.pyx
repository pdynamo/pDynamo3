"""Handle cubic splines."""

#from Serialization import RawObjectConstructor

from   pCore       import Clone 
from  .SplineError import SplineError
from ..Arrays      import Reshape

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class CubicSpline:

    def __copy__ ( self ):
        """Copying."""
        cdef CubicSpline new
        cdef CStatus     cStatus = CStatus_OK
        new         = self.__class__.Raw ( )
        new.cObject = CubicSpline_Allocate ( NULL )
        CubicSpline_AssignArrays ( new.cObject    ,
                                   self.x.cObject ,
                                   self.y.cObject ,
                                   self.h.cObject ,
                                   &cStatus       )
        if cStatus != CStatus_OK: raise SplineError ( "Error shallow cloning spline." )
        new.x       = self.x
        new.y       = self.y
        new.h       = self.h
        new.isOwner = True
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            CubicSpline_DeassignArrays (  self.cObject )
            CubicSpline_Deallocate     ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        cdef CubicSpline new
        cdef CStatus     cStatus = CStatus_OK
        new         = self.__class__.Raw ( )
        new.x       = Clone ( self.x )
        new.y       = Clone ( self.y )
        new.h       = Clone ( self.h )
        new.cObject = CubicSpline_Allocate ( NULL )
        CubicSpline_AssignArrays ( new.cObject   ,
                                   new.x.cObject ,
                                   new.y.cObject ,
                                   new.h.cObject ,
                                   &cStatus      )
        if cStatus != CStatus_OK: raise SplineError ( "Error deep cloning spline." )
        new.isOwner = True
        return new

    def __init__ ( self ):
        """Constructor."""
        self._Initialize ( )

#    def __reduce_ex__ ( self, protocol ):
#        """Pickling protocol."""
#        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.h       = None
        self.isOwner = False
        self.x       = None
        self.y       = None

    def Evaluate ( self, CReal x, CInteger spline = 0 ):
        """Evaluate a spline given an abscissa value."""
        cdef CReal    f, g, h
        cdef CStatus  cStatus = CStatus_OK
        CubicSpline_Evaluate ( self.cObject, spline, x, &f, &g, &h, &cStatus )
        if cStatus != CStatus_OK: raise SplineError ( "Error evaluating spline." )
        return ( f, g, h )

    def EvaluateN ( self, x, f = None, g = None, h = None ):
        """Evaluate all the splines given an abscissa value."""
        # . f, g and h should be arrays or lists.
        fgh = [ self.Evaluate ( x, spline = s ) for s in range ( self.numberOfSplines ) ]
        if f is not None:
            for ( i, ( fX, _, _ ) ) in enumerate ( fgh ): f[i] = fX
        if g is not None:
            for ( i, ( _, gX, _ ) ) in enumerate ( fgh ): g[i] = gX
        if h is not None:
            for ( i, ( _, _, hX ) ) in enumerate ( fgh ): h[i] = hX

    def FindMaxima ( self, CInteger spline = 0 ):
        """Find the maxima."""
        cdef RealArray1D values
        cdef CInteger    nValues
        cdef CReal       x
        cdef CStatus     cStatus = CStatus_OK
        values = RealArray1D.WithExtent ( len ( self.x ) + 1 )
        CubicSpline_FindExtrema ( self.cObject, spline, values.cObject, NULL, &nValues, NULL, &cStatus )
        if cStatus != CStatus_OK: raise SplineError ( "Error finding spline maxima." )
        return [ ( x, self.Evaluate ( x, spline = spline )[0] ) for x in values[0:min(nValues,len(values))] ]

    @classmethod
    def FromArrays ( selfClass                      ,
                     *args                          ,
                     CInteger lowerDerivative = 2   ,
                     CReal    lowerValue      = 0.0 ,
                     CInteger upperDerivative = 2   ,
                     CReal    upperValue      = 0.0 ):
        """Constructor given arrays and boundary conditions."""
        cdef CubicSpline self
        cdef RealArray1D xLocal = None
        cdef RealArray2D yLocal = None
        cdef RealArray2D hLocal
        cdef CStatus     cStatus = CStatus_OK
        # . Input arrays = y or x, y.
        if len ( args ) == 1:
            x = None
            y = args[0]
        elif len ( args ) == 2:
            x = args[0]
            y = args[1]
        else:
            raise SplineError ( "Invalid constructor arguments." )
        # . y = single spline 1D array or iterable or multispline 2D array.
        if isinstance ( y, RealArray1D ):
            yLocal = Reshape ( y, [ len ( y ), 1 ] )
        elif isinstance ( y, RealArray2D ):
            yLocal = y
        else:
            yLocal = RealArray2D.WithShape ( [ len ( y ), 1 ] )
            for ( i, v ) in enumerate ( y ): yLocal[i,0] = v
        # . x = 1D array or iterable or None in which case integral values are generated.
        if isinstance ( x, RealArray1D ):
            xLocal = x
        else:
            xLocal = RealArray1D.WithExtent ( yLocal.rows )
            if x is None:
                for i in range ( len ( xLocal ) ): xLocal[i] = float ( i )
            else:
                for ( i, v ) in enumerate ( x ): xLocal[i] = v
        hLocal       = RealArray2D.WithShape ( list ( yLocal.shape ) )
        self         = selfClass.Raw ( )
        self.cObject = CubicSpline_FromRealArrays ( xLocal.cObject  ,
                                                    yLocal.cObject  ,
                                                    hLocal.cObject  ,
                                                    lowerDerivative ,
                                                    lowerValue      ,
                                                    upperDerivative ,
                                                    upperValue      ,
                                                    &cStatus        )
        if cStatus  != CStatus_OK: raise SplineError ( "Error creating spline from arrays." )
        self.x       = xLocal
        self.y       = yLocal
        self.h       = hLocal
        self.isOwner = True
        return self

    def Integrate ( self, CReal a, CReal b, CInteger spline = 0 ):
        """Integrate a spline between the limits [a,b]."""
        cdef CReal   value
        cdef CStatus cStatus = CStatus_OK
        value = CubicSpline_Integrate ( self.cObject, spline, a, b, &cStatus )
        if cStatus != CStatus_OK: raise SplineError ( "Error integrating spline within range." )
        return value

    def IntegrateFull ( self, CInteger spline = 0 ):
        """Integrate a spline over its full range."""
        cdef CReal   value
        cdef CStatus cStatus = CStatus_OK
        value = CubicSpline_IntegrateFull ( self.cObject, spline, &cStatus )
        if cStatus != CStatus_OK: raise SplineError ( "Error integrating spline over its full range." )
        return value

    def Locate ( self, x ):
        """Locate the segment within which the abscissa lies."""
        il = 0
        iu = self.numberOfPoints - 1
        while ( iu - il ) > 1:
            im = ( iu + il ) // 2
            if x > self.x[im]: il = im
            else:              iu = im
        return il

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    # . Properties.
    @property
    def lower ( self ):
        return self.x[0]
    @property
    def numberOfPoints ( self ):
        return self.y.rows
    @property
    def numberOfSplines ( self ):
        return self.y.columns
    @property
    def upper ( self ):
        return self.x[-1]
