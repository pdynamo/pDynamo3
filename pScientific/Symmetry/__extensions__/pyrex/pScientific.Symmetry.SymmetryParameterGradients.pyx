"""Interface to the symmetry parameter gradients C type."""

from  pCore         import RawObjectConstructor
from .SymmetryError import SymmetryError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class SymmetryParameterGradients:
    """Define the gradients corresponding to a set of symmetry parameters."""

    def __dealloc__ ( self ):
        """Finalization."""
        if self.isOwner:
            SymmetryParameterGradients_Deallocate ( &self.cObject )
            self.isOwner = False

    def __getstate__ ( self ):
        """Return the state."""
        return { "dEdA"     : self.cObject.dEda     ,
                 "dEdB"     : self.cObject.dEdb     ,
                 "dEdC"     : self.cObject.dEdc     ,
                 "dEdAlpha" : self.cObject.dEdalpha ,
                 "dEdBeta"  : self.cObject.dEdbeta  ,
                 "dEdGamma" : self.cObject.dEdgamma ,
                 "dEdH"     : self.dEdH             }

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
        self.cObject.dEda     = state["dEdA"    ]
        self.cObject.dEdb     = state["dEdB"    ]
        self.cObject.dEdc     = state["dEdC"    ]
        self.cObject.dEdalpha = state["dEdAlpha"]
        self.cObject.dEdbeta  = state["dEdBeta" ]
        self.cObject.dEdgamma = state["dEdGamma"]
        self.dEdH             = state["dEdH"    ]

    def _Allocate ( self ):
        """Allocation."""
        cdef CStatus cStatus = CStatus_OK
        self.dEdH    = Matrix33.Null ( )
        self.cObject = SymmetryParameterGradients_AllocateWithMatrix ( self.dEdH.cObject, &cStatus )
        self.isOwner = True
        if cStatus != CStatus_OK: raise SymmetryError ( "Error allocating symmetry parameter gradients." )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = True
        self.dEdH    = None

    def Clear ( self ):
        """Clear the data structure."""
        SymmetryParameterGradients_Initialize ( self.cObject )

    def MakeCrystalDerivatives ( self, SymmetryParameters symmetryParameters ):
        """Make the crystal derivatives from the lattice derivatives."""
        SymmetryParameterGradients_CrystalDerivatives ( self.cObject, symmetryParameters.cObject )

    def MakeFractionalDerivatives ( self, SymmetryParameters symmetryParameters, Coordinates3 coordinates3, Coordinates3 gradients3 ):
        """Convert r/H derivatives to f/H derivatives."""
        SymmetryParameterGradients_FractionalDerivatives ( self.cObject, symmetryParameters.cObject, coordinates3.cObject, gradients3.cObject )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    @classmethod
    def Uninitialized ( selfClass ):
        """Constructor."""
        return selfClass ( )

    # . Properties.
    property dEdA:
        def __get__ ( self ): return self.cObject.dEda
    property dEdB:
        def __get__ ( self ): return self.cObject.dEdb
    property dEdC:
        def __get__ ( self ): return self.cObject.dEdc
    property dEdAlpha:
        def __get__ ( self ): return self.cObject.dEdalpha
    property dEdBeta:
        def __get__ ( self ): return self.cObject.dEdbeta
    property dEdGamma:
        def __get__ ( self ): return self.cObject.dEdgamma
