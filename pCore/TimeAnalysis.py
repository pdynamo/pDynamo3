"""Classes for time analysis."""

from  collections   import defaultdict
from .LogFileWriter import logFile       , \
                           LogFileActive
from .PrintObjects  import Align
from .Time          import CPUTime

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class Timings ( defaultdict ):
    """Timings handler."""

    def __init__ ( self, *arguments, **options ):
        """Constructor."""
        super ( Timings, self ).__init__ ( float, *arguments, **options )

    def Accumulate ( self, other ):
        """Accumulation."""
        for ( key, value ) in other.items ( ): self[key] += value

    def Clear ( self ):
        """Clear all data."""
        self.clear ( )

    def Summary ( self, log = logFile, orderByMagnitude = False, total = None, totalKey = "Total" ):
        """Statistics summary."""
        def ExtractValue ( item ): return -item[1]
        if LogFileActive ( log ):
            # . Find total.
            if total is not None: localTotal = total
            else:                 localTotal = self[totalKey]
            localTotal = self.SummaryValue ( localTotal )
            # . Construct values to output.
            keyLength = 0
            keys      = self.keys ( )
            values    = []
            for ( key, value ) in self.items ( ):
                keyLength  = max ( keyLength, len ( key ) )
                localValue = self.SummaryValue ( value )
                percentage = localValue / localTotal * 100.0
                values.append ( ( key, localValue, percentage ) )
            # . Sort.
            if orderByMagnitude: values.sort ( key = ExtractValue )
            else:                values.sort ( )
            # . Output.
            table = log.GetTable ( columns = [ max ( keyLength + 2, 20 ), 20, 15 ] )
            table.Start  ( )
            table.Title  ( self.SummaryTitle ( ) )
            table.Heading ( "Item"       )
            table.Heading ( "Time"       )
            table.Heading ( "Percentage" )
            for ( key, value, percentage ) in values:
                table.Entry (  key, align = Align.Left  )
                table.Entry ( "{:s}"   .format ( CPUTime.TimeToString ( value ) ) )
                table.Entry ( "{:8.3f}".format ( percentage                     ) )
            table.Stop ( )

    def SummaryTitle ( self ): return "Timings Summary"

    def SummaryValue ( self, value ): return value

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class TimingsAverager ( Timings ):
    """Timings averager."""

    def __init__ ( self, *arguments, **options ):
        """Constructor."""
        super ( TimingsAverager, self ).__init__ ( *arguments, **options )
        self.numberOfSamples = 0

    def Accumulate ( self, other ):
        """Accumulation."""
        super ( TimingsAverager, self ).Accumulate ( other )
        self.numberOfSamples += 1

    def Clear ( self ):
        """Clear all data."""
        super ( TimingsAverager, self ).Clear ( )
        self.numberOfSamples = 0

    def Summary ( self, **options ):
        """Summary."""
        if self.numberOfSamples <= 0: raise ValueError ( "The averager has no samples to summarize." )
        super ( TimingsAverager, self ).Summary ( **options )

    def SummaryTitle ( self ): return ( "Average Timings From {:d} Samples".format ( self.numberOfSamples ) )

    def SummaryValue ( self, value ): return ( value / float ( self.numberOfSamples ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :

    averager = TimingsAverager ( )
    for i in range ( 3 ):
        timings = Timings ( )
        for i in range ( 10 ):
            timings[repr ( i )] += float ( i )
        timings["Total"] = 100.0
        timings.Accumulate ( { repr ( i ) : float ( i ) for i in range ( 21, 25 ) } )
        timings.Summary ( totalKey = "Total" )
        timings.Summary ( totalKey = "Total", orderByMagnitude = True )
        averager.Accumulate ( timings )
    averager.Summary ( orderByMagnitude = True )
