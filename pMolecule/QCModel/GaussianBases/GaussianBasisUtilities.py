"""Gaussian basis utilities."""

from .GaussianBasisError import GaussianBasisError

#===================================================================================================================================
# . Methods.
#===================================================================================================================================
def AMLabelDecode ( label ):
    """Convert a label to an angular momentum."""
    c = label.upper ( )
    if   c == "S": return 0
    elif c == "P": return 1
    elif c == "D": return 2
    else: return ( ord ( c ) - ord ( "C" ) )

def AMLabelEncode ( l ):
    """Convert an angular momentum to a one character label."""
    if   l == 0: return "S"
    elif l == 1: return "P"
    elif l == 2: return "D"
    else: return chr ( l + ord ( "C" ) )

def ShellLabelDecode ( label ):
    """Convert a shell label to an angular momentum low/high pair."""
    isOK = ( len ( label ) > 0 )
    if isOK:
        for ( i, c ) in enumerate ( label ):
            l = AMLabelDecode ( c )
            if i == 0:
                lHigh = l
                lLow  = l
            elif l == ( lHigh + 1 ):
                lHigh = l
            else:
                isOK = False
                break
    if not isOK: raise GaussianBasisError ( "Invalid shell label." )
    return ( lLow, lHigh )

def ShellLabelEncode ( lLow, lHigh ):
    """Convert an angular momentum low/high pair to a shell label."""
    return "".join ( [ AMLabelEncode ( l ) for l in range ( lLow, lHigh+1 ) ] )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
