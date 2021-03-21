"""A container for Gaussian basis sets."""

import os, os.path

from  pCore              import Clone                    , \
                                DataType                 , \
                                RawObjectConstructor     , \
                                YAMLMappingFile_ToObject , \
                                YAMLPickleFileExtension
from  pScientific        import PeriodicTable
from  pScientific.Arrays import Array
from .QCDefinitions      import BasisRepresentation
from .QCModelError       import QCModelError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class GaussianBasisContainer:

    def __copy__ ( self ):
        """Copying."""
        cdef GaussianBasisContainer new
        new                     = self.__class__.Raw ( )
        new.cObject             = GaussianBasisContainer_Clone ( self.cObject, NULL )
        new.isOwner             = True
        new.label               = self.label
        new.uniqueEntries       = Clone ( self.uniqueEntries )
        new.basisRepresentation = self.basisRepresentation
        new._MakeBasisRepresentations ( )
        return new

    def __dealloc__ ( self ):
        """Destructor."""
        self.uniqueEntries = None
        if self.isOwner:
            GaussianBasisContainer_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getstate__ ( self ):
        """Return the state."""
        state = { "Atomic Numbers"       : self.atomicNumbers        ,
                  "Basis Representation" : self._representation.name ,
                  "Unique Entries"       : self.uniqueEntries        }
        if self.label is not None: state["Label"] = self.label
        return state

    def __init__ ( self, capacity ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( capacity )

    def __len__ ( self ):
        """The number of basis functions in the current representation."""
        return self.numberOfFunctions

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        self._CreateObject ( state["Unique Entries"], state["Atomic Numbers"] )
        self.label = state.get ( "Label", None )
        self._MakeBasisRepresentations ( )
        if "Basis Representation" in state:
            self.basisRepresentation = BasisRepresentation.__dict__[state["Basis Representation"]]

    def _Allocate ( self, capacity ):
        """Constructor."""
        self.cObject = GaussianBasisContainer_Allocate ( capacity, NULL )
        self.isOwner = True

    def _CreateObject ( self, uniqueEntries, atomicNumbers ):
        """Create the object."""
        cdef GaussianBasis entry
        cdef CInteger      i
        capacity = len ( atomicNumbers )
        self._Allocate ( capacity )
        self.uniqueEntries = uniqueEntries
        for i from 0 <= i < self.cObject.capacity:
            entry = uniqueEntries[atomicNumbers[i]]
            self.cObject.entries[i] = entry.cObject

    def _Initialize ( self ):
        """Initialization."""
        self.cObject                 = NULL
        self.isOwner                 = False
        self.label                   = None
        self.uniqueEntries           = None
        self._centerFunctionPointers = {}
        self._functionCenters        = {}
        self._numberOfFunctions      = {}
        self._representation         = None
        self._a2w                    = None
        self._w2a                    = None

    def _MakeBasisRepresentations ( self ):
        """Make the basis representations."""
        cdef IntegerArray1D centerFunctionPointers
        cdef IntegerArray1D functionCenters
        cdef CBoolean       cFlag
        cdef CInteger       n
        cdef CStatus        cStatus = CStatus_OK
        # . Initialization.
        if self.cObject != NULL:
            self._centerFunctionPointers = {}
            self._functionCenters        = {}
            self._numberOfFunctions      = {}
            # . The only representations this container can handle.
            for ( value, cFlag ) in ( ( BasisRepresentation.Actual, CFalse ) ,
                                      ( BasisRepresentation.Work  , CTrue  ) ):
                centerFunctionPointers = Array.WithExtent ( self.cObject.capacity + 1, dataType = DataType.Integer )
                GaussianBasisContainer_MakeBasisIndices ( self.cObject, cFlag, centerFunctionPointers.cObject, &cStatus )
                if cStatus != CStatus_OK: raise QCModelError ( "Error creating center function pointers." )
                n = centerFunctionPointers[-1]
                functionCenters = Array.WithExtent ( n, dataType = DataType.Integer  )
                GaussianBasisContainer_MakeBasisAtomIndices ( self.cObject, cFlag, functionCenters.cObject, &cStatus )
                if cStatus != CStatus_OK: raise QCModelError ( "Error creating function centers." )
                self._centerFunctionPointers[value] = centerFunctionPointers
                self._functionCenters       [value] = functionCenters
                self._numberOfFunctions     [value] = n

    def _MakeFunctionTransformations ( self ):
        """Make the function transformations (actual <-> work)."""
        cdef RealArray2D a2w
        cdef RealArray2D w2a
        cdef CInteger    nA, nW
        cdef CStatus     cStatus = CStatus_OK
        if ( self.cObject != NULL ) and ( ( self._a2w is None ) or ( self._w2a is None ) ):
            nA  = GaussianBasisContainer_NumberOfBasisFunctions ( self.cObject, CFalse )
            nW  = GaussianBasisContainer_NumberOfBasisFunctions ( self.cObject, CTrue  )
            w2a = Array.WithExtents ( nW, nA )
            a2w = Array.WithExtents ( nW, nA )
            GaussianBasisContainer_MakeFunctionTransformations ( self.cObject ,
                                                                 w2a.cObject  ,
                                                                 a2w.cObject  ,
                                                                 &cStatus     )
            if cStatus != CStatus_OK: raise QCModelError ( "Error creating function transformations." )
            self._a2w = a2w
            self._w2a = w2a

    @classmethod
    def FromParameterDirectory ( selfClass, basisLabel, atomicNumbers, path = None ):
        """Constructor given a basis label and a list of atomic numbers."""
        if path is None: path = os.path.join ( os.getenv ( "PDYNAMO3_PARAMETERS" ), "gaussianBasisSets", basisLabel )
        missing       = set ( )
        uniqueEntries = {}
        for atomicNumber in sorted ( set ( atomicNumbers ) ):
            entryPath = os.path.join ( path, PeriodicTable.Symbol ( atomicNumber ) + YAMLPickleFileExtension )
            try:
                basis = YAMLMappingFile_ToObject ( entryPath, GaussianBasis )
                basis.Normalize ( doReport = False ) # . With normalization.
                uniqueEntries[atomicNumber] = basis
            except: missing.add ( "{:d}".format ( atomicNumber ) )
        if len ( missing ) > 0:
            raise QCModelError ( "There are missing {:s} basis sets: {:s}.".format ( basisLabel, ", ".join ( sorted ( missing ) ) ) )
        else:
            return selfClass.FromUniqueEntries ( uniqueEntries, atomicNumbers, label = basisLabel )

    @classmethod
    def FromUniqueEntries ( selfClass, uniqueEntries, atomicNumbers, label = None ):
        """Constructor from unique entries."""
        cdef GaussianBasisContainer self
        self = selfClass.Raw ( )
        self._CreateObject ( uniqueEntries, atomicNumbers )
        if label is not None: self.label = label
        self._MakeBasisRepresentations ( )
        self._representation = BasisRepresentation.Actual # . Set the actual representation as the default.
        return self

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def RotationMatrices ( self, RealArray2D R not None ):
        """Make rotation matrices."""
        cdef GaussianBasis entry
        cdef RealArray2D   T
        cdef RealArray2D   Tc
        cdef CBoolean      cDoW2A
        cdef CStatus       cStatus = CStatus_OK
        doW2A = ( self._representation is BasisRepresentation.Actual )
        if doW2A: cDoW2A = CTrue
        else:     cDoW2A = CFalse
        L  = max ( [ entry.maximumAngularMomentum for entry in self.uniqueEntries.values ( ) ] )
        d  = ( ( L + 1 ) * ( L + 2 ) * ( L + 3 ) ) // 6
        Tc = Array.WithShape ( [ d, d ] )
        GaussianBasis_MakeLRotations ( L, R.cObject, Tc.cObject, &cStatus )
        Tn = {}
        for ( n, entry ) in self.uniqueEntries.items ( ):
            if doW2A: d = entry.numberOfFunctions
            else:     d = entry.numberOfWorkFunctions
            T = Array.WithShape ( [ d, d ] )
            GaussianBasis_MakeRotationMatrix ( entry.cObject, Tc.cObject, cDoW2A, T.cObject, &cStatus )
            Tn[n] = T
        if cStatus != CStatus_OK: raise QCModelError ( "Error making rotation matrices." )
        return [ Tn[n] for n in self.atomicNumbers ]

    @property
    def a2w ( self ):
        if self._a2w is None: self._MakeFunctionTransformations ( )
        return self._a2w

    @property
    def atomicNumbers ( self ):
        cdef CInteger i
        values = []
        if self.cObject != NULL:
            for i from 0 <= i < self.cObject.capacity:
                values.append ( self.cObject.entries[i].atomicNumber )
        return values

    @property
    def basisRepresentation ( self ):
        return self._representation
    @basisRepresentation.setter
    def basisRepresentation ( self, value ):
        """Set the basis representation."""
        if ( value is BasisRepresentation.Actual ) or ( value is BasisRepresentation.Work ):
            self._representation = value
        else:
            raise QCModelError ( "Invalid basis representation." )

    @property
    def centerFunctionPointers ( self ):
        return self._centerFunctionPointers[self._representation]

    @property
    def functionCenters ( self ):
        return self._functionCenters[self._representation]

    @property
    def nuclearCharges ( self ):
        cdef RealArray1D charges = None
        cdef CInteger    i
        if self.cObject != NULL:
            charges = Array.WithExtent ( self.cObject.capacity )
            for i from 0 <= i < self.cObject.capacity:
                charges[i] = float ( self.cObject.entries[i].atomicNumber ) # . Eventually core charge when have ECPs, etc.
        return charges

    @property
    def numberOfCenters ( self ):
        if self.cObject == NULL: return 0
        else:                    return self.cObject.capacity

    @property
    def numberOfFunctions ( self ):
        return self._numberOfFunctions[self._representation]

    @property
    def w2a ( self ):
        if self._w2a is None: self._MakeFunctionTransformations ( )
        return self._w2a
