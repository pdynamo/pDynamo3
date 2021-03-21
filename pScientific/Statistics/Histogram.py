"""Basic histogram classes."""

import math

from   pCore  import AttributableObject , \
                     Clone              , \
                     logFile            , \
                     LogFileActive      , \
                     SummarizableObject
from ..Arrays import Array

#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# . Remarks:
#
# * This will probably be converted to C as it is needed for some low level stuff by the linear algebra modules.
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#===================================================================================================================================
# . Parameter definitions.
#===================================================================================================================================
# . Nudge factor for slightly extending the length of a histogram dimension.
_DEFAULTNUDGE     = 1.0e-6

# . Tolerance for bin size checking.
_DEFAULTTOLERANCE = 1.0e-6

#===================================================================================================================================
# . Classes.
#===================================================================================================================================
# . The initial dimensions change fastest.
# . All data is currently handled as one dimensional.

#-----------------------------------------------------------------------------------------------------------------------------------
# . Histogram class.
#-----------------------------------------------------------------------------------------------------------------------------------
class RegularHistogram ( SummarizableObject ):
    """A class for regular histograms of arbitrary dimension."""

    # . No other _summarizable.
    _attributable = dict ( SummarizableObject._attributable )
    _classLabel   = "Regular Histogram"
    _attributable.update ( { "bins"               : 0    ,
                             "counts"             : list ,
                             "dimensions"         : list ,
                             "integrationElement" : 0.0  } )

    def BinData ( self, data ):
        """Bin data into the histogram."""
        # . The two versions below differ by as much as a factor of 350 (principally due to the function call in the inner loop)!
        # . Cleanest version - but VERY inefficient.
#        for i in range ( 0, len ( data ), len ( self.dimensions ) ):
#            index = 0
#            for ( v, dimension ) in zip ( data[i:], self.dimensions ):
#                i = dimension.GetBinIndex ( v )
#                if i >= 0: index += i * dimension.increment
#                else:
#                    index = -1
#                    break
#            if index >= 0: self.counts[index] += 1
        # . MUCH more efficient version but not as nice. Also uses more space.
        # . Initialization.
        tdata       = len ( data )
        ndata       = tdata // len ( self.dimensions )
        ndimensions = len ( self.dimensions )
        indices     = [ 0 for i in range ( ndata ) ]
        # . Loop over dimensions.
        for ( d, dimension ) in enumerate ( self.dimensions ):
            bins       = dimension.bins
            binSize    = dimension.binSize
            increment  = dimension.increment
            isPeriodic = dimension.isPeriodic
            lower      = dimension.lower
            period     = dimension.period
            for ( i, ( o, v ) ) in enumerate ( zip ( indices, data[d:tdata:ndimensions] ) ):
                if o >= 0:
                    if isPeriodic: v -= period * round ( v / period )
                    index = int ( ( v - lower ) / binSize )
                    if ( index >= 0 ) and ( index < bins ): indices[i] += increment * index
                    else:                                   indices[i]  = -1
        # . Finish up.
        for index in indices:
            if index >= 0: self.counts[index] += 1

    def BinMidPointIterator ( self ):
        """Return an iterator that iterates over bin mid-point values."""
        return RegularHistogramBinMidPointIterator ( self )

    @staticmethod
    def CheckDimensionIndices ( tocheck ):
        """Check a list containing dimension indices."""
        # . Check the dimension indices.
        indices = list ( set ( tocheck ) )
        indices.sort ( )
        if ( len ( indices ) <= 0 ) or ( indices[0] < 0 ) or ( indices[-1] >= len ( self.dimensions ) ): raise TypeError ( "Invalid dimension indices." )
        return indices

    def DataIndexIterator ( self, dimensions = None, indices = None ):
        """Return an iterator over data indices."""
        return RegularHistogramDataIndexIterator ( self, dimension = dimensions, indices = indices )

    @classmethod
    def FromDimensions ( selfClass, dimensions ):
        """Constructor given suitable histogram dimensions."""
        # . Check the dimensions.
        for dimension in dimensions:
            if not isinstance ( dimension, RegularHistogramDimension ) or not dimension.IsValid ( ): raise TypeError ( "Invalid dimension type." )
        # . Create the object.
        self = selfClass ( )
        self.dimensions = dimensions
        # . Determine some remaining attributes - including those in dimensions.
        increment =   1
        volume    = 1.0
        for dimension in self.dimensions:
            dimension.increment = increment
            increment *= dimension.bins
            volume    *= dimension.binSize
        self.bins               = increment
        self.integrationElement = volume
        # . Histogram counts.
        self.InitializeCounts ( )
        return self

    def IndexIterator ( self, dimensions = None, indices = None ):
        """Return an iterator over histogram indices."""
        return RegularHistogramIndexIterator ( self, dimension = dimensions, indices = indices )

    def InitializeCounts ( self ):
        """Initialize the counts."""
        self.counts = self.bins * [ 0 ]

    def MakeReducedHistogram ( self, toReduce ):
        """Create a new histogram by reducing along certain dimensions."""
        # . Get the indices.
        toReduce = self.__class__.CheckDimensionIndices ( toReduce )
        # . Create the new dimensions.
        dimensions = []
        for ( i, dimension ) in enumerate ( self.dimensions ):
            if i not in toReduce: dimensions.append ( Clone ( dimension ) )
        # . Create the object.
        new = self.__class__.FromDimensions ( dimensions )
        # . Finish up.
        if len ( self.counts > 0 ): new.counts = self.ReduceData ( self.counts, toReduce )
        return new

    def Normalize ( self, binData ):
        """Normalize a vector with bin data such that its total integral is one."""
        if len ( binData ) == self.bins:
            volume = sum ( binData ) * self.integrationElement
            binData.Scale ( 1.0 / volume )
        else:
            raise TypeError ( "Invalid bin data argument length {:d}.".format ( len ( binData ) ) )

    def ReduceData ( self, toReduce, binData ):
        """Reduce bin data along certain dimensions."""
        # . Get the indices.
        toReduce = self.__class__.CheckDimensionIndices ( toReduce )
        # . Check bin data.
        if ( len ( binData ) == self.bins ): raise TypeError ( "Invalid length for bin data argument {:d}.".format ( len ( binData ) ) )
        # . Find the indices to keep.
        bins   = 1
        toKeep = []
        for ( i, dimension ) in enumerate ( self.dimensions ):
            if i not in toRemove:
                bins *= dimension.bins
                toKeep.append ( i )
        # . Allocate space.
        new = Array.WithExtent ( bins )
        new.Set ( 0.0 )
        # . Loop over the indices of the two different dimensions.
        for ( i, indices ) in enumerate ( self.IndexIterator ( dimensions = toKeep ) ):
            s = 0
            for j in self.DataIndexIterator ( dimensions = toRemove, indices = indices ): s += binData[j]
            new[i] = s
        # . Finish up.
        return new

    def Summary ( self, log = logFile ):
        """Summary."""
        super ( RegularHistogram, self ).Summary ( log = log )
        self.SummaryDimensions ( log = log )

    def SummaryDimensions ( self, log = logFile ):
        """Summary of dimension data."""
        if LogFileActive ( log ):
            length = max ( max ( [ len ( dimension.label ) for dimension in self.dimensions ] ) + 2, 10 )
            table  = log.GetTable ( columns = [ length, 10 ] + 4 * [ 20 ] )
            table.Start ( )
            table.Title ( "{:s} Dimensions".format ( self.label ) )
            table.Heading ( "Label"    )
            table.Heading ( "Bins"     )
            table.Heading ( "Bin Size" )
            table.Heading ( "Lower"    )
            table.Heading ( "Upper"    )
            table.Heading ( "Period"   )
            for dimension in self.dimensions:
                table.Entry ( dimension.label  )
                table.Entry ( "{:d}"  .format ( dimension.bins    ) )
                table.Entry ( "{:.6g}".format ( dimension.binSize ) )
                table.Entry ( "{:.6g}".format ( dimension.lower   ) )
                table.Entry ( "{:.6g}".format ( dimension.upper   ) )
                if dimension.isPeriodic: table.Entry ( "{:.6g}".format ( dimension.period ) )
                else:                    table.Entry ( "" )
            table.Stop ( )

    def SummaryItems ( self ):
        """Summary items."""
        return [ ( "Bins",   "{:d}".format ( self.bins           ) ) ,
                 ( "Counts", "{:d}".format ( sum ( self.counts ) ) ) ]

    def ToTextFileWithData ( self, path, data, format = None ):
        """Write the histogram to a text file with additional data."""
        # . Data.
        n = len ( data )
        for d in data:
            if len ( d ) != self.bins: raise ValueError ( "Invalid data array dimension." )
        # . Format.
        if format is None: format = " ".join ( ( self.rank + n ) * [ "{:20.3f}" ] ) + "\n"
        # . Output.
        dataFile = open ( path, "w" )
        for ( i, midPoint ) in enumerate ( self.BinMidPointIterator ( ) ):
            items = list ( midPoint ) + [ d[i] for d in data ]
            dataFile.write ( format.format ( *items ) )
        dataFile.close ( )

#-----------------------------------------------------------------------------------------------------------------------------------
# . Iterator class.
#-----------------------------------------------------------------------------------------------------------------------------------
class RegularHistogramBinMidPointIterator:
    """An iterator for a regular histogram that returns the bin mid-point values."""

    def __init__ ( self, histogram ):
        """Constructor."""
        self.histogram = histogram
        self.indices   = [ -1 ] + ( len ( histogram.dimensions ) - 1 ) * [ 0 ]
        self.ntimes    = 0

    def __iter__ ( self ):
        """Iterator interface."""
        return self

    def __next__ ( self ): return self.next ( )

    def next ( self ):
        """Return the next mid-point."""
        if self.ntimes >= self.histogram.bins: raise StopIteration
        # . Find the index of the new point.
        for ( d, dimension ) in enumerate ( self.histogram.dimensions ):
            index = self.indices[d] + 1
            if index >= dimension.bins:
                self.indices[d] = 0
            else:
                self.indices[d] = index
                break
        # . Calculate the midpoint.
        midpoint = []
        for ( index, dimension ) in zip ( self.indices, self.histogram.dimensions ):
            midpoint.append ( float ( index ) * dimension.binSize + dimension.midPointOrigin )
        # . Finish up.
        self.ntimes += 1
        return midpoint

#-----------------------------------------------------------------------------------------------------------------------------------
# . Iterator class.
#-----------------------------------------------------------------------------------------------------------------------------------
class RegularHistogramDataIndexIterator:
    """An iterator over data indices."""

    def __init__ ( self, histogram, dimensions = None, indices = None ):
        """Constructor."""
        self.histogram = histogram
        self.ntimes    = 0
        # . Dimensions.
        if dimensions is None: self.dimensions = range ( len ( histogram.dimensions ) )
        else:                  self.dimensions = dimensions
        # . Indices.
        if indices    is None: self.indices = len ( histogram.dimensions ) * [ -1 ]
        else:                  self.indices = indices

    def __iter__ ( self ):
        """Iterator interface."""
        return self

    def __next__ ( self ): return self.next ( )

    def next ( self ):
        """Return the next data index."""
        if self.ntimes >= self.histogram.bins: raise StopIteration
        # . Find the index of the new point.
        for d in self.dimensions:
            index = self.indices[d] + 1
            if index >= self.histogram.dimensions[d].bins:
                self.indices[d] = 0
            else:
                self.indices[d] = index
                break
        # . Calculate the data index.
        dataindex = 0
        for ( index, dimension ) in zip ( self.indices, self.histogram.dimensions ):
            dataindex += index * dimension.increment
        #. Finish up.
        self.ntimes += 1
        return dataindex

#-----------------------------------------------------------------------------------------------------------------------------------
# . Histogram dimension class.
#-----------------------------------------------------------------------------------------------------------------------------------
class RegularHistogramDimension ( AttributableObject ):
    """A class to represent a regular histogram dimension.

    A periodic dimension is allowed.
    """

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "bins"           :     0 ,
                             "binSize"        :  None ,
                             "increment"      :  None ,
                             "isPeriodic"     : False ,
                             "label"          :  None ,
                             "lower"          :  None ,
                             "midPointOrigin" :  None ,
                             "period"         :  None ,
                             "upper"          :  None } )

    @classmethod
    def FromData ( selfClass, data, binSize = 0.0, bins = 0, label = None, nudge = _DEFAULTNUDGE, period = None, tolerance = _DEFAULTTOLERANCE ):
        """Constructor from data.

        Either the number of bins can be specified or an approximate bin size.
        """
        # . Check for bins or binSize (the former takes precedence).
        if   isinstance ( bins,    int   ) and ( bins    > 0   ): binSize = None
        elif isinstance ( binSize, float ) and ( binSize > 0.0 ): bins    = None
        else: raise TypeError ( "Invalid bin/binSize specification." )
        # . Check for periodicity.
        isPeriodic = ( period is not None )
        # . Find the maximum and minimum of the data.
        lower = data[0]
        upper = data[0]
        for value in data[1:]:
            if isPeriodic: value -= period * round ( value / period )
            lower = min ( lower, value )
            upper = max ( upper, value )
        # . Apply a nudge to ensure all points are covered.
        if nudge > 0.0:
            factor = 0.5 * nudge * ( upper - lower )
            lower -= factor
            upper += factor
        # . An approximate binSize was specified so determine the number of bins
        if bins is None:
            length = upper - lower
            bins   = int ( ( length ) / binSize )
            if ( math.fabs ( length / ( float ( bins ) * binSize ) - 1.0 ) > tolerance ): bins += 1
            delta  = 0.5 * math.fabs ( float ( bins ) * binSize - length )
            lower -= delta
            upper += delta
        # . Adjust for periodic systems.
        if isPeriodic:
            factor = 0.5 * period
            upper = min ( upper,   period )
            lower = max ( lower, - period )
        # . Construct the object.
        self = selfClass ( )
        # . Set attributes.
        self.bins           = bins
        self.binSize        = ( upper - lower ) / float ( bins )
        self.label          = label
        self.lower          = lower
        self.midPointOrigin = lower + 0.5 * self.binSize
        self.upper          = upper
        if isPeriodic:
            self.isPeriodic = True
            self.period     = period
        return self

    def GetBinIndex ( self, value ):
        """Return the bin index given a value."""
        if self.isPeriodic: value -= self.period * round ( value / self.period )
        index = int ( ( value - self.lower ) / self.binSize )
        if ( index < 0 ) or ( index >= self.bins ): index = -1
        return index

    def IsValid ( self, tolerance = _DEFAULTTOLERANCE ):
        """Check for valid dimension attributes."""
        QOK = ( self.bins > 0 ) and ( self.binSize is not None ) and ( self.binSize > 0.0 ) and ( self.label is not None ) and ( ( not self.isPeriodic ) or ( self.isPeriodic and ( self.period is not None ) ) )
        if QOK:
            maxdev = max ( math.fabs ( ( self.upper - self.lower ) / ( float ( self.bins ) * self.binSize ) - 1.0 ), math.fabs ( 2.0 * ( self.midPointOrigin - self.lower ) / self.binSize - 1.0 ) )
            QOK    = ( maxdev < tolerance )
        return QOK

#-----------------------------------------------------------------------------------------------------------------------------------
# . Iterator class.
#-----------------------------------------------------------------------------------------------------------------------------------
class RegularHistogramIndexIterator:
    """An iterator over indices."""

    def __init__ ( self, histogram, dimensions = None, indices = None ):
        """Constructor."""
        self.histogram = histogram
        self.ntimes    = 0
        # . Dimensions.
        if dimensions is None: self.dimensions = range ( len ( histogram.dimensions ) )
        else:                  self.dimensions = dimensions
        # . Indices.
        if indices    is None: self.indices = len ( histogram.dimensions ) * [ -1 ]
        else:                  self.indices = indices

    def __iter__ ( self ):
        """Iterator interface."""
        return self

    def __next__ ( self ): return self.next ( )

    def next ( self ):
        """Return the next index."""
        if self.ntimes >= self.histogram.bins: raise StopIteration
        # . Find the index of the new point.
        for d in self.dimensions:
            index = self.indices[d] + 1
            if index >= self.histogram.dimensions[d].bins:
                self.indices[d] = 0
            else:
                self.indices[d] = index
                break
        #. Finish up.
        self.ntimes += 1
        return self.indices

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
