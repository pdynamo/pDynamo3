"""Status utilities."""

# . Note that inspect does not work in Cython (at the moment) so cannot use "inspect.stack ( )[1][3]" to get caller name.

from .CoreError import CoreError

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
_Status_Header   = "C library error"
_Status_ToString = { CStatus_AlgorithmError        : "algorithm error"         ,
                     CStatus_IndexOutOfRange       : "index out of range"      ,
                     CStatus_InvalidArgument       : "invalid argument"        ,
                     CStatus_InvalidArrayOperation : "invalid array operation" ,
                     CStatus_MathError             : "math error"              ,
                     CStatus_NonConformableArrays  : "non-conformable arrays"  ,
                     CStatus_OutOfMemory           : "out of memory"           }

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def Status_Check ( CStatus status ):
    """Check the status flag."""
    if status != CStatus_OK:
        string  = _Status_ToString.get ( status, None )
        if string is None: message = "*** {:s} ***".format ( _Status_Header )
        else:              message = "*** {:s}: {:s} ***".format ( _Status_Header, string )
        raise CoreError ( message )
