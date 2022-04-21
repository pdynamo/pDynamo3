"""Antisymmetric matrices."""

from  pCore      import RawObjectConstructor
from .ArrayError import ArrayError
from .Slicing    import ProcessIntegerSlice

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class AntisymmetricMatrix:

    # . Methods.
    def __copy__ ( self ):
        """Copying."""
        cdef AntisymmetricMatrix clone
        cdef CStatus             cStatus = CStatus_OK
        clone         = self.__class__.Raw ( )
        clone.block   = self.block
        clone.cObject = AntisymmetricMatrix_CloneShallow ( self.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error cloning array." )
        return clone

    def __dealloc__ ( self ):
        """Finalization."""
        self.block     = None
        self._iterator = None
        AntisymmetricMatrix_Deallocate ( &self.cObject )

    def __deepcopy__ ( self, memo ):
        """Copying."""
        clone = self.__class__.WithExtent ( self.extent )
        self.CopyTo ( clone )
        return clone

    def __getattr__ ( self, name ):
        """Handle unknown attributes."""
        return getattr ( self.iterator, name )

    def __getitem__ ( self, indices ):
        """Get an item."""
        cdef CInteger i, j
        cdef CReal    value
        ( i, j ) = ProcessIntegerSlice ( self.shape, indices )
        value    = AntisymmetricMatrix_GetItem ( self.cObject, i, j, NULL )
        return value

    def __getstate__ ( self ):
        """Return the state."""
        return { "block"  : self.block  ,
                 "extent" : self.extent }

    def __init__ ( self, extent ):
        """Constructor with extent."""
        self._Initialize ( )
        self._Allocate ( extent )

    def __iter__ ( self ):
        return self._MakeIterator ( )

    def __len__ ( self ):
        """Return the size of the symmetricmatrix."""
        return self.size

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setitem__ ( self, indices, CReal value ):
        """Set an item."""
        cdef CInteger i, j
        ( i, j ) = ProcessIntegerSlice ( self.shape, indices )
        AntisymmetricMatrix_SetItem ( self.cObject, i, j, value, NULL )

    def __setstate__ ( self, state ):
        """Set the state."""
        self._Allocate ( state["extent"], block = state["block"] )

    def _Allocate ( self, CInteger extent, RealBlock block = None ):
        """Allocation."""
        cdef CRealBlock *cBlock  = NULL
        cdef CStatus     cStatus = CStatus_OK
        if block is     None:  block = RealBlock.WithCapacity ( AntisymmetricMatrix_ViewSize ( extent ) )
        if block is not None: cBlock = block.cObject
        self.block   = block
        self.cObject = AntisymmetricMatrix_FromExtentBlock ( extent, cBlock, CTrue, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error allocating antisymmetric matrix." )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject   = NULL
        self.block     = None
        self._iterator = None

    def _MakeIterator ( self ):
        """Make the default array iterator."""
        cdef RealIterator iterator
        cdef CIterator   *cIterator = NULL
        cdef CStatus      cStatus   = CStatus_OK
        iterator           = RealIterator.Raw ( )
        iterator.cIterator = AntisymmetricMatrix_MakeIterator ( self.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Unable to make iterator." )
        iterator.block     = self.block
        iterator._SetDataPointer ( )
        iterator.Reset ( )
        return iterator

    def AbsoluteMaximum ( self ):
        """Return the maximum absolute value in the matrix."""
        return AntisymmetricMatrix_AbsoluteMaximum ( self.cObject )

    def Add ( self, AntisymmetricMatrix other, CReal scale = 1.0 ):
        """Add a scaled matrix."""
        cdef CStatus cStatus = CStatus_OK
        AntisymmetricMatrix_Add ( self.cObject, scale, other.cObject, NULL )
        if cStatus != CStatus_OK: raise ArrayError ( "Error adding antisymmetric matrices." )

    def CopyTo ( self, AntisymmetricMatrix other ):
        """Copying."""
        AntisymmetricMatrix_CopyTo ( self.cObject, other.cObject, NULL )

    def MakeCommutatorSS ( self, SymmetricMatrix a not None ,
                                 SymmetricMatrix b not None ,
                                 RealArray2D     mA = None  ,
                                 RealArray2D     mB = None  ,
                                 RealArray2D     mC = None  ):
        """Make the commutator a * b - b * a."""
        cdef CStatus cStatus = CStatus_OK
        if ( mA is not None ) and ( mB is not None ) and ( mC is not None ):
            AntisymmetricMatrix_CommutatorSS_Fast      ( self.cObject ,
                                                         a.cObject    ,
                                                         b.cObject    ,
                                                         mA.cObject   ,
                                                         mB.cObject   ,
                                                         mC.cObject   ,
                                                         &cStatus     )
        else:
            AntisymmetricMatrix_CommutatorSS_Reference ( self.cObject ,
                                                         a.cObject    ,
                                                         b.cObject    ,
                                                         &cStatus     )
        if cStatus != CStatus_OK: raise ArrayError ( "Error making commutator." )

    def MakeCommutatorSSS ( self, SymmetricMatrix a not None ,
                                  SymmetricMatrix b not None ,
                                  SymmetricMatrix c not None ):
        """Make the commutator a * b * c - c * b * a."""
        cdef CStatus  cStatus    = CStatus_OK
        AntisymmetricMatrix_CommutatorSSS ( self.cObject ,
                                            a.cObject    ,
                                            b.cObject    ,
                                            c.cObject    ,
                                            &cStatus     )
        if cStatus != CStatus_OK: raise ArrayError ( "Error making commutator." )

    def MakeCommutatorXSSSX ( self, SymmetricMatrix a not None ,
                                    SymmetricMatrix b not None ,
                                    SymmetricMatrix c not None ,
                                    RealArray2D     x not None ,
                                    RealArray2D     u not None ,
                                    RealArray2D     v not None ,
                                    RealArray2D     w not None ,
                                            xTranspose = False ):
        """Make the commutator x^T * ( a * b * c - c * b * a ) * x with optional x transposition."""
        cdef CBoolean cXTranpose
        cdef CStatus  cStatus    = CStatus_OK
        if xTranspose: cXTranspose = CTrue
        else:          cXTranspose = CFalse
        AntisymmetricMatrix_CommutatorTSSST ( self.cObject ,
                                              a.cObject    ,
                                              b.cObject    ,
                                              c.cObject    ,
                                              x.cObject    ,
                                              cXTranspose  ,
                                              u.cObject    ,
                                              v.cObject    ,
                                              w.cObject    ,
                                              &cStatus     )
        if cStatus != CStatus_OK: raise ArrayError ( "Error making commutator." )

    def MakeCommutatorXSSY ( self, SymmetricMatrix a not None ,
                                   SymmetricMatrix b not None ,
                                   RealArray2D     x not None ,
                                   RealArray2D     y not None ,
                                   RealArray2D     u not None ,
                                   RealArray2D     v not None ,
                                   RealArray2D     w not None ,
                                           xTranspose = False ,
                                           yTranspose = False ):
        """Make the commutator x * a * b * y - y^T * b * a * x^T with optional x and y transposition."""
        cdef CBoolean cXTranpose
        cdef CBoolean cYTranpose
        cdef CStatus  cStatus    = CStatus_OK
        if xTranspose: cXTranspose = CTrue
        else:          cXTranspose = CFalse
        if yTranspose: cYTranspose = CTrue
        else:          cYTranspose = CFalse
        AntisymmetricMatrix_CommutatorXSSY ( self.cObject ,
                                             a.cObject    ,
                                             b.cObject    ,
                                             x.cObject    ,
                                             y.cObject    ,
                                             cXTranspose  ,
                                             cYTranspose  ,
                                             u.cObject    ,
                                             v.cObject    ,
                                             w.cObject    ,
                                             &cStatus     )
        if cStatus != CStatus_OK: raise ArrayError ( "Error making commutator." )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def Scale ( AntisymmetricMatrix self, CReal value ):
        """Scale all the elements of a symmetricmatrix."""
        AntisymmetricMatrix_Scale ( self.cObject, value )

    def Set ( AntisymmetricMatrix self, CReal value ):
        """Set all the elements of a symmetricmatrix."""
        AntisymmetricMatrix_Set ( self.cObject, value )

    def Trace ( self ):
        """Trace."""
        return 0.0

    def TraceOfProduct ( self, AntisymmetricMatrix other not None ):
        """Trace of product."""
        cdef CReal   value
        cdef CStatus cStatus = CStatus_OK
        value = AntisymmetricMatrix_TraceOfProduct ( self.cObject, other.cObject, &cStatus )
        if cStatus != CStatus_OK: raise ArrayError ( "Error determining trace of product." )
        return value

    def Transform ( self, RealArray2D         x not None ,
                          AntisymmetricMatrix b not None ,
                                      xTranspose = False ):
        """Transform by a matrix (b = x^T * a * x) or its transpose (b = x * a * x^T)."""
        cdef CBoolean cXTranpose
        cdef CStatus  cStatus    = CStatus_OK
        if xTranspose: cXTranspose = CTrue
        else:          cXTranspose = CFalse
        AntisymmetricMatrix_Transform ( self.cObject ,
                                        x.cObject    ,
                                        cXTranspose  ,
                                        b.cObject    ,
                                        &cStatus     )
        if cStatus != CStatus_OK: raise ArrayError ( "Error transforming matrix." )

    @classmethod
    def WithExtent ( selfClass, extent ):
        """Constructor with extent."""
        return selfClass ( extent )

    # . Properties.
    @property
    def columns ( self ): return self.extent
    @property
    def extent ( self ):
        if self.cObject == NULL: return 0
        else:                    return self.cObject.extent
    @property
    def isSquare ( self ): return True
    @property
    def iterator ( self ):
        if self._iterator is None:
            self._iterator = self._MakeIterator ( )
        return self._iterator
    @property
    def rank ( self ): return 2
    @property
    def rows ( self ): return self.extent
    @property
    def shape ( self ): return [ self.extent, self.extent ]
    @property
    def size ( self ):
        if self.cObject == NULL: return 0
        else:                    return self.cObject.size
