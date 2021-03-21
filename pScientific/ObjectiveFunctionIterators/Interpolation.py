"""Interpolation functions."""

from math import fabs, sqrt

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def CubicMinimizerFGFG ( xP, fP, gP, xQ, fQ, gQ ):
    """Determine the minimizer of a cubic given two function and gradient values."""
    theta  = 3.0 * ( fP - fQ ) / ( xQ - xP ) + gP + gQ
    s      = max ( fabs ( theta ), fabs ( gP ), fabs ( gQ ) )
    gamma2 = ( theta / s )**2 - ( gP / s ) * ( gQ / s )
    isOK   = ( gamma2 >= 0.0 )
    if isOK:
        gamma = s * sqrt ( gamma2 )
        if xP > xQ: gamma *= -1.0
    else:
        gamma = 0.0
    p    =       gamma + theta - gP
    q    = 2.0 * gamma + gQ    - gP
    isOK = isOK and ( q != 0.0 )
    if isOK: x = xP + ( p / q ) * ( xQ - xP )
    else:    x = 0.0
    return ( isOK, x )

def QuadraticExtremumFGF ( xP, fP, gP, xQ, fQ ):
    """Determine the extremum of a quadratic given two function values and one gradient."""
    delta = ( xP - xQ )
    q     = ( fP - fQ - delta * gP )
    isOK  = ( q != 0.0 )
    if isOK: x = xP + 0.5 * ( gP * delta**2 / q )
    else:    x = 0.0
    return ( isOK, x )

def QuadraticExtremumGG ( xP, gP, xQ, gQ ):
    """Determine the extremum of a quadratic given two gradient values."""
    q    = ( gP - gQ )
    isOK = ( q != 0.0 )
    if isOK: x = xP + ( gP / q ) * ( xQ - xP )
    else:    x = 0.0
    return ( isOK, x )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
