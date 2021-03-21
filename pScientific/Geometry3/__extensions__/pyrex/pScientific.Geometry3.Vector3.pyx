"""Real 3-D vector."""

from .Geometry3Error import Geometry3Error

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class Vector3 ( RealArray1D ):

    def _Allocate ( self, extent ):
        """Allocation."""
        if extent != 3: raise Geometry3Error ( "Invalid Vector3 shape." )
        super ( Vector3, self )._Allocate ( extent )

    def _GetRawArray1D ( self, extent ):
        """Get a raw 1-D array."""
        if extent == 3: return Vector3.RawWithCObject     ( )
        else:           return RealArray1D.RawWithCObject ( )

    def Cross ( Vector3 self, Vector3 other ):
        """In-place cross-product."""
        Vector3_CrossProduct ( self.cObject, other.cObject )

    @classmethod
    def Null ( selfClass ):
        """Constructor."""
        cdef Vector3 self
        self = selfClass.Uninitialized ( )
        self.Set ( 0.0 )
        return self

    @classmethod
    def Uninitialized ( selfClass ):
        """Constructor."""
        return selfClass ( 3 )

    @classmethod
    def WithValues ( selfClass, a, b, c ):
        """Constructor."""
        cdef Vector3 self
        self = selfClass.Uninitialized ( )
        self[0] = a ; self[1] = b ; self[2] = c
        return self
