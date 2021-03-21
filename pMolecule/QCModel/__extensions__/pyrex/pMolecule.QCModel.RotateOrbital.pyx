"""Orbital rotation functions."""

from  pScientific.Arrays import Array
from .QCModelError       import QCModelError

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def RotateOrbital ( IntegerArray1D orbitalBasisIndices not None ,
                    Matrix33       rotation            not None ,
                    IntegerArray1D mapping             not None ,
                    RealArray1D    inOrbital           not None ,
                    RealArray1D    outOrbital          not None ):
    """Rotate an orbital."""
    CRotateOrbital ( orbitalBasisIndices.cObject ,
                     rotation.cObject            ,
                     mapping.cObject             ,
                     inOrbital.cObject           ,
                     outOrbital.cObject          )

def RotateOrbital_MakeRotations ( IntegerArray1D orbitalBasisIndices not None ,
                                  Matrix33       rotation            not None ):
    """Make orbital rotations."""
    cdef RealArray2D T
    cdef CInteger    L
    cdef CStatus     cStatus = CStatus_OK
    # . Make basic transformation.
    d = 0
    N = len ( orbitalBasisIndices ) - 1
    for i in range ( N ): d = max ( d, orbitalBasisIndices[i+1] - orbitalBasisIndices[i] )
    if   d == 1: L = 0
    elif d == 4: L = 1
    elif d == 9: L = 2
    else: raise QCModelError ( "Invalid basis set angular momentum." )
    T = Array.WithShape ( [ d, d ] )
    CRotateOrbital_MakeLRotations ( L, rotation.cObject, T.cObject, &cStatus )
    if cStatus != CStatus_OK: raise QCModelError ( "Error making basis set rotation matrices." )
    # . Create transformation matrix array.
    Ts = []
    for i in range ( N ):
        d = orbitalBasisIndices[i+1] - orbitalBasisIndices[i]
        Ts.append ( T[0:d,0:d] )
    return Ts
