"""Real 3x3 matrices."""

import math

from .Geometry3Error import Geometry3Error

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class Matrix33 ( RealArray2D ):

    def _Allocate ( self, rows, columns ):
        """Allocation."""
        if ( rows != 3 ) or ( columns != 3 ): raise Geometry3Error ( "Invalid Matrix33 shape." )
        super ( Matrix33, self )._Allocate ( rows, columns )

    def ApplyTo ( self, Vector3 vector3 ):
        """Apply the matrix to another object - usually a vector."""
        Matrix33_ApplyToVector3 ( self.cObject, vector3.cObject )

    @classmethod
    def Identity ( selfClass ):
        """Constructor."""
        self = selfClass.Null ( )
        self[0,0] = 1.0 ; self[1,1] = 1.0 ; self[2,2] = 1.0
        return self

    def Invert ( self, Matrix33 other ):
        """Invert the matrix and put the result in other."""
        Matrix33_Invert ( other.cObject, self.cObject )

    @classmethod
    def MakeRandomRotation ( selfClass, randomNumberGenerator, tolerance = 1.0e-20 ):
        """Construct a random rotation."""
        self = selfClass.Uninitialized ( )
        self.RandomRotation ( randomNumberGenerator, tolerance = tolerance )
        return self

    @classmethod
    def MakeRotationAboutAxis ( selfClass, CReal angle, Vector3 axis ):
        """Construct a rotation corresponding to a rotation about a normalized axis."""
        self = selfClass.Uninitialized ( )
        self.RotationAboutAxis ( angle, axis )
        return self

    @classmethod
    def Null ( selfClass ):
        """Constructor."""
        self = selfClass.Uninitialized ( )
        self.Set ( 0.0 )
        return self

    def PostMultiplyBy ( self, Matrix33 other ):
        """Postmultiply by another matrix."""
        Matrix33_PostMultiplyBy ( self.cObject, other.cObject )

    def PreMultiplyBy ( self, Matrix33 other ):
        """Premultiply by another matrix."""
        Matrix33_PreMultiplyBy ( self.cObject, other.cObject )

    def RandomRotation ( self, randomNumberGenerator, tolerance = 1.0e-20 ):
        """Fill the matrix with a random rotation."""
        ralpha = 2.0 * math.pi * ( randomNumberGenerator.NextReal ( ) - 0.5 )
        raxis  = Vector3.Null ( )
        norm2  = -1.0
        while norm2 < tolerance:
            for i in range ( 3 ):
                raxis[i] = 2.0 * ( randomNumberGenerator.NextReal ( ) - 0.5 )
            norm2 = raxis.Norm2 ( )
        raxis.Normalize ( )
        self.RotationAboutAxis ( ralpha, raxis )

    def Reflection ( self, Vector3 normal ):
        """Fill the matrix with a reflection through the origin."""
        Matrix33_Reflection ( &self.cObject, normal.cObject )

    def RotationAboutAxis ( self, CReal angle, Vector3 axis ):
        """Fill the matrix with a rotation about a normalized axis."""
        cdef CReal x, y, z
        x = RealArray1D_GetItem ( axis.cObject, 0, NULL )
        y = RealArray1D_GetItem ( axis.cObject, 1, NULL )
        z = RealArray1D_GetItem ( axis.cObject, 2, NULL )
        Matrix33_RotationAboutAxis ( &self.cObject, angle, x, y, z )

    @classmethod
    def Uninitialized ( selfClass ):
        """Constructor."""
        return selfClass ( 3, 3 )

    # . Properties.
    @property
    def isIdentity         ( self ): return ( Matrix33_IsIdentity         ( self.cObject ) == CTrue )
    @property
    def isImproperRotation ( self ): return ( Matrix33_IsImproperRotation ( self.cObject ) == CTrue )
    @property
    def isOrthogonal       ( self ): return ( Matrix33_IsOrthogonal       ( self.cObject ) == CTrue )
    @property
    def isProperRotation   ( self ): return ( Matrix33_IsProperRotation   ( self.cObject ) == CTrue )
