"""Rys quadrature functions."""

from .GaussianBasisError import GaussianBasisError

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def RysQuadrature_MaximumRoots ( ):
    return CRysQuadrature_MaximumRoots ( )

def RysQuadrature_Roots ( CInteger nRoots, CReal x ):
    """Find the quadrature roots and weights."""
    cdef CInteger       r
    cdef CRysQuadrature cRoots
    #cdef CStatus        cStatus = CStatus_OK
    if nRoots > RysQuadrature_MaximumRoots ( ): raise GaussianBasisError ( "Too many Rys roots, {:d}, requested.".format ( nRoots ) )
    CRysQuadrature_Roots ( &cRoots, nRoots, x )
    #if cStatus != CStatus_OK: raise GaussianBasisError ( "Error finding Rys roots." )
    roots   = []
    weights = []
    for r from 0 <= r < nRoots:
        roots.append   ( cRoots.roots  [r] )
        weights.append ( cRoots.weights[r] )
    return ( roots, weights )
