"""Utilities for handling stick and continuous spectra."""

import math

from pCore     import AttributableObject , \
                      Clone
from ..        import Constants
from ..Arrays  import Array
from ..Splines import CubicSpline

# . Units are not currently treated so any conversions need to be handled explicitly.

# . Surface temperature of the sun in K.
Constants_Sun_Surface_Temperature = 5778.0

#===================================================================================================================================
# . Data ranges and sets.
#===================================================================================================================================
class DataRange ( AttributableObject ):
    """A range over a data set."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "size"   : 0   ,
                             "start"  : 0.0 ,
                             "stop"   : 0.0 ,
                             "stride" : 0.0 } )

    def __iter__ ( self ): return DataRangeIterator ( self )

    def __len__ ( self ): return self.size

    @classmethod
    def FromSizeStartStop ( selfClass, size, start, stop ):
        """Constructor given size, start and stop."""
        self        = selfClass ( )
        self.size   = max ( size, 2 )
        self.start  = start
        self.stop   = stop
        self.stride = ( stop - start ) / float ( self.size - 1 )
        return self

class DataRangeIterator:
    """An iterator over a data range."""

    def __init__ ( self, target ):
        """Constructor."""
        self.size   = target.size
        self.start  = target.start
        self.stride = target.stride
        self.Reset ( )

    def __iter__ ( self ): return self

    def __next__ ( self ): return self.next ( )

    def next ( self ):
        self.counter += 1
        self.value   += self.stride
        if self.counter >= self.size: raise StopIteration
        else: return self.value

    def Reset ( self ):
        """Reset the iterator."""
        self.counter = -1
        self.value   = self.start - self.stride

class DataSet ( AttributableObject ):
    """A set of data with associated information."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "data"  : None ,
                             "label" : None ,
                             "units" : None } )

    def __len__ ( self ): return len ( self.data )

    @classmethod
    def FromDataRange ( selfClass, dataRange ):
        """Constructor given a data range."""
        size = len ( dataRange )
        self = selfClass.WithSize ( size )
        for ( i, d ) in enumerate ( dataRange ): self.data[i] = d
        return self

    @classmethod
    def WithSize ( selfClass, size ):
        """Constructor given a size."""
        self      = selfClass ( )
        self.data = Array.WithExtent ( size )
        return self

#===================================================================================================================================
# . Line shapes.
#===================================================================================================================================
# . FWHM = full band width at half the maximum height.
class LineShape:
    """Base class for line shapes."""

    # . center and fwhm have the same units.

    def __init__ ( self, center, fwhm ):
        """Constructor."""
        self.center = center
        self.fwhm   = fwhm
        self._SetConstants ( )

    def _SetConstants ( self ):
        """Set up some constants for the line shape."""
        pass

    def Height ( self, value ):
        """Determine the line shape height at a given value."""
        return 0.0

    def ResetCenter ( self, center ):
        """Reset the center."""
        self.center = center

    def TailValue ( self, tolerance ):
        """Return (x-c) such that f(x)/fmaximum is equal to tolerance."""
        return 0.0

    @property
    def maximum ( self ): return 0.0

class GaussianLineShape ( LineShape ):
    """Gaussian line shape."""

    def _SetConstants ( self ):
        """Set up some constants for the line shape."""
        self.exponent  = 4.0 * math.log ( 2.0 ) / self.fwhm**2
        self.prefactor = math.sqrt ( self.exponent / math.pi )

    def Height ( self, value ):
        """Determine the line shape height at a given value."""
        return self.prefactor * math.exp ( - self.exponent * ( value - self.center )**2 )

    def TailValue ( self, tolerance ):
        """Return (x-c) such that f(x)/fmaximum is equal to tolerance."""
        return math.sqrt ( max ( 0.0, math.log ( max ( 1.0, 1.0 / tolerance ) ) / self.exponent ) )

    @property
    def maximum ( self ): return self.prefactor

class LorentzianLineShape ( LineShape ):
    """Lorentzian line shape."""

    def _SetConstants ( self ):
        """Set up some constants for the line shape."""
        hwhm = 0.5 * self.fwhm
        self.denominator = hwhm**2
        self.prefactor   = hwhm / math.pi

    def Height ( self, value ):
        """Determine the line shape height at a given value."""
        return self.prefactor / ( self.denominator + ( value - self.center )**2 )

    def TailValue ( self, tolerance ):
        """Return (x-c) such that f(x)/fmaximum is equal to tolerance."""
        return math.sqrt ( max ( 0.0, 1.0 / tolerance - self.denominator ) )

    @property
    def maximum ( self ): return ( self.prefactor / self.denominator )

#===================================================================================================================================
# . Spectra.
#===================================================================================================================================
class Spectrum ( AttributableObject ):
    """Base class for spectra."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "abscissae" : None ,
                             "label"     : None ,
                             "ordinates" : None } )

    def __getstate__ ( self ):
        """Get the state."""
        data = []
        for ( x, y ) in zip ( self.abscissae.data, self.ordinates.data ): data.append ( [ x, y ] )
        return { "Data Fields" : [ "Wavelength", "Amplitude" ] ,
                 "Data Values" : data                          ,
                 "Label"       : self.label                    ,
                 "Units"       : { "Amplitude" : repr ( self.ordinates.units ) , "Wavelength" :  repr ( self.abscissae.units ) } }

    def __iter__ ( self ):
        """Return an iterator over abscissae and ordinates."""
        return SpectrumIterator ( self )

    def __len__ ( self ):
        if self.abscissae is None: return 0
        else:                      return len ( self.abscissae.data )

    def __setstate__ ( self, state ):
        """Set the state."""
        data           = state.get ( "Data Values", []   )
        self.label     = state.get ( "Label"      , None )
        size           = len ( data )
        self.abscissae = DataSet.WithSize ( size )
        self.ordinates = DataSet.WithSize ( size )
        for ( i, ( x, y ) ) in enumerate ( data ):
            self.abscissae.data[i] = x
            self.ordinates.data[i] = y

    @classmethod
    def Raw ( selfClass ): return selfClass ( )

    @classmethod
    def WithAbscissae ( selfClass, abscissae ):
        """Constructor given a set of abscissae."""
        self           = selfClass ( )
        self.abscissae = abscissae
        self.ordinates = DataSet.WithSize ( len ( abscissae ) )
        self.ordinates.data.Set ( 0.0 )
        return self

    @classmethod
    def WithSize ( selfClass, size ):
        """Constructor given a size."""
        self           = selfClass ( )
        self.abscissae = DataSet.WithSize ( size )
        self.ordinates = DataSet.WithSize ( size )
        return self

class SpectrumIterator:
    """An iterator over the abscissae and ordinates of a spectrum."""

    def __init__ ( self, target ):
        """Constructor."""
        self.size  = len ( target )
        self.xData = target.abscissae.data
        self.yData = target.ordinates.data
        self.Reset ( )

    def __iter__ ( self ): return self

    def __next__ ( self ): return self.next ( )

    def next ( self ):
        self.counter += 1
        if self.counter >= self.size: raise StopIteration
        else: return ( self.xData[self.counter], self.yData[self.counter] )

    def Reset ( self ):
        """Reset the iterator."""
        self.counter = -1

class ContinuousSpectrum ( Spectrum ):
    """Continous spectrum."""

    @classmethod
    def FromStickSpectrum ( selfClass, stickSpectrum, abscissae, lineShape, label = None ):
        """Constructor given a stick spectrum, a set of abscissae and a line shape."""
        # . Allocation and initialization.
        self                 = selfClass.WithAbscissae ( abscissae )
        self.label           = label
        self.abscissae.label = stickSpectrum.abscissae.label
        self.abscissae.units = stickSpectrum.abscissae.units
        self.ordinates.label = stickSpectrum.ordinates.label
        self.ordinates.units = stickSpectrum.ordinates.units
        # . Save the line shape center.
        c0 = lineShape.center
        # . Loop over stick peaks and data points.
        xData = self.abscissae.data
        yData = self.ordinates.data
        for ( c, m ) in stickSpectrum:
            lineShape.ResetCenter ( c )
            for ( i, x ) in enumerate ( xData ):
                yData[i] += m * lineShape.Height ( x )
        # . Finish up.
        lineShape.ResetCenter ( c0 )
        return self

    def Normalize ( self ):
        """Normalization so that the integral of the spectrum is one."""
        if self.integral != 0.0: self.ScaleOrdinates ( 1.0 / self.integral )

    def OverlapSpectrum ( self, others ):
        """Create an overlap spectrum from the current spectrum and a set of other spectra."""
        # . Currently self is used as the template for the abscissae.
        result = None
        if len ( others ) == 0:
            result = self
        else:
            # . Limits and splines.
            lower   = self.spline.lower
            upper   = self.spline.upper
            splines = [ self.spline ]
            for other in others:
                lower = max ( lower, other.spline.lower )
                upper = min ( upper, other.spline.upper )
                splines.append ( other.spline )
            # . Abscissae.
            abscissae = []
            for x in self.abscissae.data:
                if ( x >= lower ) and ( x <= upper ): abscissae.append ( x ) 
            abscissae.sort ( )
            # . Construct the overlap.
            if len ( abscissae ) > 0:
                result = self.__class__.WithSize ( len ( abscissae ) )
                for ( i, x ) in enumerate ( abscissae ):
                    y = 1.0
                    for spline in splines: y *= spline.Evaluate ( x )[0]
                    result.abscissae.data[i] = x
                    result.ordinates.data[i] = y
        return result

    def ScaleAbscissae ( self, scale ):
        """Scale the abscissae."""
        self.abscissae.data.Scale ( scale )
        self._integral = None
        self._spline   = None

    def ScaleOrdinates ( self, scale ):
        """Scale the ordinates."""
        self.ordinates.data.Scale ( scale )
        self._integral = None
        self._spline   = None

    def TransmissionSpectrum ( self, factor = 1.0, illuminant = None ):
        """Create a transmission spectrum from the current spectrum given an optional absorption factor and illuminant."""
        f        = max ( factor, 0.0 )
        sMaximum = max ( self.ordinates.data )
        sSpline  = self.spline
        if illuminant is None:
            iMaximum = sMaximum
        else:
            iMaximum = max ( illuminant.ordinates.data )
            iSpline  = illuminant.spline
        result = self.__class__.WithAbscissae ( Clone ( self.abscissae ) )
        for ( i, x ) in enumerate ( result.abscissae.data ):
            if illuminant is None: vI = 1.0
            else:                  vI = iSpline.Evaluate ( x )[0] / iMaximum
            vS = sSpline.Evaluate ( x )[0] / sMaximum
            result.ordinates.data[i] = max ( 0.0, iMaximum * ( vI - f * vS ) )
        return result

    @property
    def integral ( self ):
        result = getattr ( self, "_integral", None )
        if result is None:
            result = self.spline.IntegrateFull ( )
            setattr ( self, "_integral", result )
        return result

    @property
    def spline ( self ):
        result = getattr ( self, "_spline", None )
        if result is None:
            result = CubicSpline.FromArrays ( self.abscissae.data, self.ordinates.data )
            setattr ( self, "_spline", result )
        return result

class StickSpectrum ( Spectrum ):
    """Stick spectrum."""
    pass

#===================================================================================================================================
# . Specific types of continuous spectra.
#===================================================================================================================================
class BlackBodySpectrum ( ContinuousSpectrum ):
    """Black-body spectrum at a given temperature."""

    _attributable = dict ( ContinuousSpectrum._attributable )
    _attributable.update ( { "temperature" : Constants_Sun_Surface_Temperature } )

    def SetOrdinates ( self ):
        """Set the ordinate values."""
        a     = ( Constants.Planck * Constants.Speed_Of_Light) / ( Constants.Boltzmann * self.temperature )
        b     = ( 2.0 * Constants.Planck * Constants.Speed_Of_Light**2 )
        nmToM = 1.0e-9
        for ( i, x ) in enumerate ( self.abscissae.data ):
            l = x * nmToM
            e = a / l
            y = 0.0
            if e <= 500.0:
                y = b / ( math.pow ( l, 5 ) * ( math.exp ( e ) - 1.0 ) )
            self.ordinates.data[i] = y

    @classmethod
    def WithAbscissae ( selfClass, abscissae, temperature = Constants_Sun_Surface_Temperature ):
        """Constructor given abscissae."""
        self             = super ( BlackBodySpectrum, selfClass ).WithAbscissae ( abscissae )
        self.temperature = temperature
        self.SetOrdinates ( )
        return self

class RayleighSpectrum ( ContinuousSpectrum ):
    """Rayleigh spectrum."""

    @classmethod
    def FromIlluminant ( selfClass, illuminant ):
        """Constructor given an illuminant."""
        self                 = selfClass.WithSize ( len ( illuminant ) )
        self.label           = "Rayleigh Spectrum of " + illuminant.label
        self.abscissae.label = illuminant.abscissae.label
        self.abscissae.units = illuminant.abscissae.units
        self.ordinates.label = illuminant.ordinates.label
        self.ordinates.units = illuminant.ordinates.units
        illuminant.abscissae.data.CopyTo ( self.abscissae.data )
        self.SetOrdinates ( )
        self.ordinates.data.Multiply ( illuminant.ordinates.data )
        return self

    def SetOrdinates ( self ):
        """Set the ordinate values."""
        # . The spectrum is scaled so that it has a magnitude of 1 at 555 nm.
        x0 = 555.0
        for ( i, x ) in enumerate ( self.abscissae.data ):
            self.ordinates.data[i] = math.pow ( x / x0, -4 )

    @classmethod
    def WithAbscissae ( selfClass, abscissae ):
        """Constructor given abscissae."""
        self = super ( RayleighSpectrum, selfClass ).WithAbscissae ( abscissae )
        self.SetOrdinates ( )
        return self

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass


