"""Classes and functions for handling time."""

import time

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class CPUTime:
    """A class for measuring CPU time."""

    def __init__ ( self ):
        """Constructor."""
        self.start = time.time ( )

    def Current ( self ):
        """Return the current CPU time as a float."""
        return ( time.time ( ) - self.start )

    def CurrentAsString ( self, compact = True ):
        """Return the current CPU time as a string."""
        value  = self.Current ( )
        return CPUTime.TimeToString ( value, compact = compact )

    @staticmethod
    def TimeToString ( time, compact = True ):
        """Convert a floating point time to a string."""
        value  = time
        fields = []
        for ( size, tag ) in ( ( 60*60*24, "Day" ), ( 60*60, "Hour" ), ( 60, "Minute" ) ):
            if compact: tag = tag[0:1].lower ( )
            ( newf, value ) = divmod ( value, size )
            newi = int ( round ( newf ) )
            if not ( ( newi == 0 ) and ( len ( fields ) == 0 ) ):
                if compact: fields.append ( "{:2d}{:1s}".format ( newi, tag ) )
                else:
                    if newi != 1: tag += "s"
                    fields.append ( "{:2d} {:s}".format ( newi, tag ) )
        tag = "Seconds"
        if compact:
            tag = tag[0:1].lower ( )
            fields.append ( "{:5.3f}{:1s}".format ( value, tag ) )
        else:
            fields.append ( "{:5.3f} {:s}".format ( value, tag ) )
        if compact: return " " .join ( fields )
        else:       return ", ".join ( fields )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :

    # . Print CPU times.
    cputime = CPUTime ( )
    for i in range ( 20 ):
        sum ( range ( 1000000 ) )
        print ( "CPU Times: {:2d}     {:s}".format ( i, cputime.CurrentAsString ( ) ) )
