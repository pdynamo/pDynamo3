"""Basic statistical analysis classes and functions.

The Statistics class was inspired by the module from the Python cookbook by Chad J. Schroeder.

Note that the population variance is defined using 1/N whereas the sample variance is defined using 1/(N-1).
"""

import math

from   pCore  import AttributableObject
from ..Arrays import Array

#===============================================================================
# . Class.
#===============================================================================
class Statistics:
    """Class for statistical analysis."""

    # . Public methods.
    def __init__ ( self, data ):
        """Constructor."""
        # . Check that sample is iterable.
        try:    x = iter ( data )
        except: raise TypeError ( "The |data| argument must be iterable." )
        # . Set some attributes.
        self.data             = data
        self.__dict__["size"] = len ( data )
        # . Check for zero size array in constructor.
        if self.size <= 0: raise ValueError ( "There is no data in the array." )

    def __len__ ( self ): return self.size

    def BlockAveraging ( self, maximumBlockSize = None ):
        """Block averaging following Grossfield and Zuckerman 2009."""
        # . The standard errors are calculated as a function of block size.
        # . The plateau (asymptote) value gives the converged estimate of the error.
        # . For data obtained from a trajectory, the correlation time of the data
        # . can be approximated as ( simulation time / number of independent samples )
        # . where the latter is determined from the block size at the plateau value.
        # . Initialization.
        N          = self.size
        if maximumBlockSize is None: maximumBlockSize = N
        nUpper     = min ( maximumBlockSize, N )
        blockMeans = Array.WithExtent ( N          ) ; blockMeans.Set ( 0.0 )
        blockSEs   = Array.WithExtent ( nUpper + 1 ) ; blockSEs.Set   ( 0.0 )
        # . Loop over block sizes.
        for n in range ( 2, nUpper + 1 ):
            M = N // n
            for b in range ( M ):
                blockMeans[b] = Statistics ( self.data[b*n:(b+1)*n] ).mean
            blockSEs[n] = Statistics ( blockMeans[0:M] ).standardDeviation / math.sqrt ( M )
        # . Finish up.
        return blockSEs

    def Count ( self, value, tolerance = 0.0 ):
        """Count the number of elements within |tolerance| of value."""
        n = 0
        for datum in self.data:
            if ( math.fabs ( datum - value ) <= tolerance ): n += 1
        return n

    # . Pass in block sizes directly? E.g. blockSizes = range ( 0, 101 ).
    # . Must then return ( b, s ) instead of just s.
    def StatisticalInefficiency ( self, maximumBlockSize = 100 ):
        """Determine the statistical inefficiency."""
        # . Initialization.
        n = self.size
        u = min ( maximumBlockSize, n )
        m = Array.WithExtent ( n + 1 // 2 ) ; m.Set ( 0.0 )
        s = Array.WithExtent ( u + 1      ) ; s.Set ( 0.0 ) ; s[1] = 1.0
        # . Loop over block sizes.
        v1 = self.variance
        for b in range ( 2, u + 1 ):
            f = ( n // b )
            for i in range ( f ): m[i] = Statistics ( self.data[i*b:(i+1)*b] ).mean
            s[b] = ( float ( b ) * Statistics ( m[0:f] ).variance ) / v1
        # . Finish up.
        return s

    # . Private methods.
    def __GetAbsoluteMaximum ( self ):
        """Absolute maximum."""
        value = self.__dict__.get ( "absoluteMaximum", None )
        if value is None:
            value = max ( self.maximum, math.fabs ( self.minimum ) )
            self.__dict__["absoluteMaximum"] = value
        return value

    def __GetMaximum ( self ):
        """Maximum."""
        value = self.__dict__.get ( "maximum", None )
        if value is None:
            value = max ( self.data )
            self.__dict__["maximum"] = value
        return value

    def __GetMean ( self ):
        """Mean."""
        value = self.__dict__.get ( "mean", None )
        if value is None:
            if self.size > 0: value = self.sum / self.size
            else:             value = 0.0
            self.__dict__["mean"] = value
        return value

    def __GetMinimum ( self ):
        """Minimum."""
        value = self.__dict__.get ( "minimum", None )
        if value is None:
            value = min ( self.data )
            self.__dict__["minimum"] = value
        return value

    def __GetSize ( self ):
        """Size."""
        return self.__dict__["size"]

    def __GetSum ( self ):
        """Sum."""
        value = self.__dict__.get ( "sum", None )
        if value is None:
            value = sum ( self.data )
            self.__dict__["sum"] = value
        return value

    def __GetStandardDeviation ( self ):
        """Population standard deviation."""
        value = self.__dict__.get ( "standarddeviation", None )
        if value is None:
            value = math.sqrt ( self.variance )
            self.__dict__["standarddeviation"] = value
        return value

    def __GetUnsignedMean ( self ):
        """Unsigned Mean."""
        value = self.__dict__.get ( "unsignedMean", None )
        if value is None:
            value = 0.0
            if self.size > 0:
                for datum in self.data:
                    value += math.fabs ( datum )
                value /= self.size
            self.__dict__["unsignedMean"] = value
        return value

    def __GetVariance ( self ):
        """Population variance."""
        value = self.__dict__.get ( "variance", None )
        if value is None:
            value = 0.0
            if self.size > 0:
                for datum in self.data:
                    value += ( datum - self.mean )**2
                value /= self.size
            self.__dict__["variance"] = value
        return value

    # . Properties.
    absoluteMaximum   = property ( __GetAbsoluteMaximum,   None, None, "Absolute Maximum."              )
    maximum           = property ( __GetMaximum,           None, None, "Maximum."                       )
    mean              = property ( __GetMean,              None, None, "Mean."                          )
    minimum           = property ( __GetMinimum,           None, None, "Minimum."                       )
    size              = property ( __GetSize,              None, None, "Size."                          )
    standardDeviation = property ( __GetStandardDeviation, None, None, "Population Standard Deviation." )
    sum               = property ( __GetSum,               None, None, "Sum."                           )
    unsignedMean      = property ( __GetUnsignedMean,      None, None, "Unsigned Mean."                 )
    variance          = property ( __GetVariance,          None, None, "Population Variance."           )

#===============================================================================
# . Class.
#===============================================================================
class StatisticsAccumulator ( AttributableObject ):
    """A class for computing basic statistical properties of data as they are accumulated.

    The data themselves are not stored - only the current value.
    """

    _attributable = dict ( AttributableObject._attributable )
    _converter    = float
    _attributable.update ( { "absoluteMaximum" : None ,
                             "absoluteMinimum" : None ,
                             "current"         : None ,
                             "maximum"         : None ,
                             "minimum"         : None ,
                             "size"            : 0    ,  
                             "sum"             : 0.0  ,
                             "sumSquared"      : 0.0  ,
                             "unsignedMaximum" : None ,
                             "unsignedMinimum" : None ,
                             "unsignedSum"     : 0.0  } )

    def Accumulate ( self, item ):
        """Accumulate a value."""
        try:    value = self.__class__._converter ( item )
        except: raise TypeError ( "Invalid type for array element. Should be " + self.__class__._converter.__name__ + " but is " + type ( item ).__name__ + "." )
        # . Treat the value.
        self.size        += 1
        # . Full value.
        self.current      = value
        self.sum         += value
        self.sumSquared  += value * value
        if self.maximum is None: self.maximum = value
        else:                    self.maximum = max ( self.maximum, value )
        if self.minimum is None: self.minimum = value
        else:                    self.minimum = min ( self.minimum, value )
        # . Absolute values.
        value = math.fabs ( value )
        self.unsignedSum += value
        if self.unsignedMaximum is None: self.unsignedMaximum = value
        else:                            self.unsignedMaximum = max ( self.unsignedMaximum, value )
        if self.unsignedMinimum is None: self.unsignedMinimum = value
        else:                            self.unsignedMinimum = min ( self.unsignedMinimum, value )
        # . Maximum and minimum.
        if self.absoluteMaximum is None: self.absoluteMaximum = value
        else:                            self.absoluteMaximum = max ( self.absoluteMaximum, value )
        if self.absoluteMinimum is None: self.absoluteMinimum = value
        else:                            self.absoluteMinimum = min ( self.absoluteMinimum, value )

    # . Private methods.
#    def __GetAbsoluteMaximum ( self ):
#        """Absolute maximum."""
#        return self.__dict__["absoluteMaximum"]

#    def __GetAbsoluteMinimum ( self ):
#        """Absolute minimum."""
#        return self.__dict__["absoluteMinimum"]

    def __GetMean ( self ):
        """Mean."""
        if self.size > 0: return self.sum / self.size
        else:             return 0.0

#    def __GetSize ( self ):
#        """Size."""
#        return self.__dict__["size"]

    def __GetStandardDeviation ( self ):
        """Population standard deviation."""
        return math.sqrt ( self.variance )

    def __GetUnsignedMean ( self ):
        """Unsigned Mean."""
        if self.size > 0: return self.unsignedSum / self.size
        else:             return 0.0

    def __GetVariance ( self ):
        """Population variance."""
        if self.size > 0: return max ( self.sumSquared / self.size - self.mean**2, 0.0 )
        else:             return 0.0

    # . Properties.
#    absoluteMaximum   = property ( __GetAbsoluteMaximum,   None, None, "Absolute Maximum."              )
#    absoluteMinimum   = property ( __GetAbsoluteMinimum,   None, None, "Absolute Minimum."              )
    mean              = property ( __GetMean,              None, None, "Mean."                          )
#    size              = property ( __GetSize,              None, None, "Size."                          )
    standardDeviation = property ( __GetStandardDeviation, None, None, "Population Standard Deviation." )
    unsignedMean      = property ( __GetUnsignedMean,      None, None, "Unsigned Mean."                 )
    variance          = property ( __GetVariance,          None, None, "Population Variance."           )

#===============================================================================
# . Testing.
#===============================================================================
if __name__ == "__main__" :
    pass
