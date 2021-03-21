"""Interface to the symmetry parameters C type."""

from  pCore         import logFile              , \
                           LogFileActive        , \
                           RawObjectConstructor
from .SymmetryError import SymmetryError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class SymmetryParameters:
    """Define a set of symmetry parameters."""

    def __copy__ ( self ):
        """Copying."""
        state = self.__getstate__ ( )
        new   = self.__class__.Raw ( )
        new.__setstate__ ( state )
        return new

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            SymmetryParameters_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getstate__ ( self ):
        """Return the state."""
        return { "a"     : self.cObject.a     ,
                 "b"     : self.cObject.b     , 
                 "c"     : self.cObject.c     ,
                 "alpha" : self.cObject.alpha ,
                 "beta"  : self.cObject.beta  ,
                 "gamma" : self.cObject.gamma }

    def __init__ ( self ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate   ( )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        self._Allocate ( )
        a     = state["a"    ]
        b     = state["b"    ]
        c     = state["c"    ]
        alpha = state["alpha"]
        beta  = state["beta" ]
        gamma = state["gamma"]
        self.SetCrystalParameters ( a, b, c, alpha, beta, gamma )

    def _Allocate ( self ):
        """Allocation."""
        cdef CStatus cStatus = CStatus_OK
        self.H        = Matrix33.Null ( )
        self.inverseH = Matrix33.Null ( )
        self.cObject  = SymmetryParameters_AllocateWithMatrices ( self.H.cObject, self.inverseH.cObject, &cStatus )
        self.isOwner  = True
        if cStatus != CStatus_OK: raise SymmetryError ( "Error allocating symmetry parameters." )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject  = NULL
        self.isOwner  = True
        self.H        = None
        self.inverseH = None

    def CenterCoordinatesByAtom ( self, Coordinates3 coordinates3, Selection selection = None ):
        """Center coordinates by atom."""
        cdef CSelection *cSelection
        cdef CStatus     cStatus = CStatus_OK
        if selection is None: cSelection = NULL
        else:                 cSelection = selection.cObject
        SymmetryParameters_CenterCoordinates3ByIndex ( self.cObject, cSelection, coordinates3.cObject, &cStatus )
        if cStatus != CStatus_OK: raise SymmetryError ( "Error centering coordinates by atom." )

    def CenterCoordinates3ByFreeIsolate ( self                                     ,
                                          SelectionContainer isolates     not None ,
                                          BooleanBlock       freeIsolates          ,
                                          Coordinates3       coordinates3 not None ):
        """Center coordinates3 by free isolates."""
        cdef CBooleanBlock *cFreeIsolates = NULL
        cdef CStatus        cStatus       = CStatus_OK
        if freeIsolates is not None: cFreeIsolates = freeIsolates.cObject
        SymmetryParameters_CenterCoordinates3ByFreeIsolate ( self.cObject         ,
                                                             isolates.cObject     ,
                                                             cFreeIsolates        ,
                                                             coordinates3.cObject ,
                                                             &cStatus             )
        if cStatus != CStatus_OK: raise SymmetryError ( "Error centering coordinates by free isolates." )

    def CopyTo ( self, SymmetryParameters other not None ):
        """Copying."""
        SymmetryParameters_CopyTo ( self.cObject, other.cObject )

    def CenterCoordinatesByIsolate ( self, Coordinates3 coordinates3, SelectionContainer isolates, Selection selection = None ):
        """Center coordinates by isolate."""
        cdef CSelection *cSelection
        cdef CStatus     cStatus = CStatus_OK
        if selection is None: cSelection = NULL
        else:                 cSelection = selection.cObject
        SymmetryParameters_CenterCoordinates3ByIsolate ( self.cObject, isolates.cObject, cSelection, coordinates3.cObject, &cStatus )
        if cStatus != CStatus_OK: raise SymmetryError ( "Error centering coordinates by isolate." )

    def IsMinimumImageConventionSatisfied ( self, length ):
        """Is the minimum image convention satisfied given a length."""
        return ( SymmetryParameters_IsMinimumImageConventionSatisfied ( self.cObject, length ) == CTrue )

    def MakeMinimumImageVector3 ( self, Vector3 r, Vector3 dR = None ):
        """Apply the minimum image convention to a vector."""
        cdef CReal *cDR = NULL
        if dR is not None: cDR = RealArray1D_PointerToData ( dR.cObject )
        SymmetryParameters_MakeMinimumImageVector ( self.cObject, RealArray1D_PointerToData ( r.cObject ), cDR )

    @classmethod
    def FromCrystalParameters ( selfClass, a, b, c, alpha, beta, gamma ):
        """Constructor from a set of crystal parameters."""
        self = selfClass ( )
        self.SetCrystalParameters ( a, b, c, alpha, beta, gamma )
        return self

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def SetCrystalParameters ( self, a, b, c, alpha, beta, gamma ):
        """Set the crystal parameters."""
        if self.cObject == NULL: self._Allocate ( )
        SymmetryParameters_SetCrystalParameters ( self.cObject, a, b, c, alpha, beta, gamma )

    @classmethod
    def Uninitialized ( selfClass ):
        """Constructor."""
        return selfClass ( )

    # . Properties.
    property a:
        def __get__ ( self ): return self.cObject.a
    property b:
        def __get__ ( self ): return self.cObject.b
    property c:
        def __get__ ( self ): return self.cObject.c
    property alpha:
        def __get__ ( self ): return self.cObject.alpha
    property beta:
        def __get__ ( self ): return self.cObject.beta
    property gamma:
        def __get__ ( self ): return self.cObject.gamma
    property volume:
        def __get__ ( self ): return SymmetryParameters_Volume ( self.cObject )
