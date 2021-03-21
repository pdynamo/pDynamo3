"""A DFT functional model."""

from  pCore          import Clone                     , \
                            RawObjectConstructor
from .DFTDefinitions import DFTFunctionalsFromOptions
from .QCModelError   import QCModelError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class DFTFunctionalModel:

    def __copy__ ( self ):
        """Copying."""
        cdef CStatus            cStatus = CStatus_OK
        cdef DFTFunctionalModel new
        new         = self.__class__.Raw ( )
        new.cObject = DFTFunctionalModel_Clone ( self.cObject, &cStatus )
        new.ids     = Clone ( self.ids )
        new.isOwner = True
        if cStatus != CStatus_OK: raise QCModelError ( "Error cloning DFT functional model." )
        return new

    def __dealloc__ ( self ):
        """Destructor."""
        if self.isOwner:
            DFTFunctionalModel_Deallocate ( &self.cObject )
            self.ids     = None
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getstate__ ( self ):
        """Return the state."""
        
        isSpinRestricted = ( self.cObject == NULL ) or ( self.cObject.isSpinRestricted == CTrue )
        return { "IDs"                : self.ids         ,
                 "Is Spin Restricted" : isSpinRestricted }
 
    def __init__ ( self, numberOfFunctionals ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( numberOfFunctionals )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        cdef CBoolean cIsSpinRestricted
        cdef CStatus  cStatus = CStatus_OK
        if state.get ( "Is Spin Restricted", True ): cIsSpinRestricted = CTrue
        else:                                        cIsSpinRestricted = CFalse
        self.ids     = state["IDs"]
        self.cObject = DFTFunctionalModel_MakeFromIDs ( self.ids.cObject, cIsSpinRestricted, &cStatus )
        self.isOwner = True
        if cStatus != CStatus_OK: raise QCModelError ( "Error setting the state of DFT functional model." )

    def _Allocate ( self, numberOfFunctionals ):
        """Constructor."""
        cdef CStatus  cStatus = CStatus_OK
        self.cObject = DFTFunctionalModel_Allocate ( numberOfFunctionals, &cStatus )
        self.isOwner = True
        if cStatus != CStatus_OK: raise QCModelError ( "Error allocating DFT functional model." )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.ids     = None
        self.isOwner = False

    # . HF is accepted but will return None.
    @staticmethod
    def FindIDs ( options, separator = "/" ):
        """Find functional IDs given a set of options."""
        ids   = None
        items = DFTFunctionalsFromOptions ( options, separator = separator )
        if len ( items ) > 0:
            ids = IntegerArray1D.WithExtent ( len ( items ) )
            for ( i, item ) in enumerate ( items ): ids[i] = item.value
        return ids

    @classmethod
    def FromIDs ( selfClass, IntegerArray1D ids not None, isSpinRestricted = True ):
        """Constructor given a set of functional IDs."""
        cdef CBoolean           cIsSpinRestricted
        cdef CStatus            cStatus = CStatus_OK
        cdef DFTFunctionalModel self
        if isSpinRestricted: cIsSpinRestricted = CTrue
        else:                cIsSpinRestricted = CFalse
        self         = selfClass.Raw ( )
        self.ids     = ids
        self.cObject = DFTFunctionalModel_MakeFromIDs ( self.ids.cObject, cIsSpinRestricted, &cStatus )
        self.isOwner = True
        if cStatus != CStatus_OK: raise QCModelError ( "Error constructing DFT functional model from IDs." )
        return self

    @classmethod
    def FromOptions ( selfClass, options, isSpinRestricted = True ):
        """Constructor from options."""
        ids = selfClass.FindIDs ( options )
        if ids is None: return None
        else:           return selfClass.FromIDs ( ids, isSpinRestricted = isSpinRestricted )

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    @property
    def exchangeScaling ( self ):
        return DFTFunctionalModel_ExchangeScaling ( self.cObject )

    @property
    def numberOfFunctionals ( self ):
        if self.ids is None: return 0
        else:                return len ( self.ids )

    @property
    def order ( self ):
        if self.cObject == NULL: return -1
        else:                    return self.cObject.order

