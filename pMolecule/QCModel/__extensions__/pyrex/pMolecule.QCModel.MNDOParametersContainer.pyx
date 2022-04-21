"""A container for MNDO parameters."""

import os, os.path

from  pCore                     import Clone                    , \
                                       logFile                  , \
                                       LogFileActive            , \
                                       RawObjectConstructor     , \
                                       YAMLMappingFile_ToObject , \
                                       YAMLPickleFileExtension
from  pScientific               import PeriodicTable
from  pScientific.Arrays        import Array
from  pScientific.LinearAlgebra import OrthogonalizationMethod
from .GaussianBases             import GaussianBasis
from .QCModelError              import QCModelError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class MNDOParametersContainer:

    def __copy__ ( self ):
        """Copying."""
        cdef MNDOParametersContainer new
        new               = self.__class__.Raw ( )
        new.cObject       = MNDOParametersContainer_Clone ( self.cObject, NULL )
        new.isOwner       = True
        new.label         = self.label
        new.uniqueEntries = Clone ( self.uniqueEntries  )
        return new

    def __dealloc__ ( self ):
        """Destructor."""
        self.uniqueEntries = None
        if self.isOwner:
            MNDOParametersContainer_Deallocate ( &self.cObject )
            self.isOwner = False

    def __deepcopy__ ( self, memo ):
        """Copying."""
        return self.__copy__ ( )

    def __getstate__ ( self ):
        """Return the state."""
        state = { "Atomic Numbers" : self.atomicNumbers ,
                  "Unique Entries" : self.uniqueEntries }
        if self.label is not None: state["Label"] = self.label
        return state

    def __init__ ( self, capacity ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( capacity )

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        self._CreateObject ( state["Unique Entries"], state["Atomic Numbers"] )
        self.label = state.get ( "Label", None )

    def _Allocate ( self, capacity ):
        """Constructor."""
        self.cObject = MNDOParametersContainer_Allocate ( capacity, NULL )
        self.isOwner = True

    def _CheckDiatomicTerms ( self ):
        """Check for missing diatomic terms."""
        numbers  = set ( self.uniqueEntries.keys ( ) )
        warnings = set ( )
        for n in numbers:
            elements = self.uniqueEntries[n].diatomicTerms
            if elements is not None:
                missing = numbers.difference ( set ( elements ) )
                for m in missing:
                    warnings.add ( ( max ( m, n ), min ( m, n ) ) )
        if ( len ( warnings ) > 0 ) and LogFileActive ( logFile ):
            logFile.Paragraph ( "Possible missing MNDOdiatomic core-core interaction parameters: {:s}.".format ( repr ( warnings ) ) )

    def _CreateObject ( self, uniqueEntries, atomicNumbers ):
        """Create the object."""
        cdef MNDOParameters entry
        cdef CInteger       i
        capacity = len ( atomicNumbers )
        self._Allocate ( capacity )
        self.uniqueEntries = uniqueEntries
        for i from 0 <= i < self.cObject.capacity:
            entry = uniqueEntries[atomicNumbers[i]]
            self.cObject.entries[i] = entry.cObject

    def _Initialize ( self ):
        """Initialization."""
        self.cObject       = NULL
        self.isOwner       = False
        self.label         = None
        self.uniqueEntries = None

    @classmethod
    def FromParameterDirectory ( selfClass, method, atomicNumbers, path = None ):
        """Constructor given a method label and a list of atomic numbers."""
        if path is None: path = os.path.join ( os.getenv ( "PDYNAMO3_PARAMETERS" ), "mndoParameters", method )
        missing       = set ( )
        uniqueEntries = {}
        for atomicNumber in sorted ( set ( atomicNumbers ) ):
            entryPath = os.path.join ( path, PeriodicTable.Symbol ( atomicNumber ) + YAMLPickleFileExtension )
            try   : uniqueEntries[atomicNumber] = YAMLMappingFile_ToObject ( entryPath, MNDOParameters )
            except: missing.add ( "{:d}".format ( atomicNumber ) )
        if len ( missing ) > 0:
            raise QCModelError ( "There are missing {:s} parameter sets: {:s}.".format ( method, ", ".join ( sorted ( missing ) ) ) )
        return selfClass.FromUniqueEntries ( uniqueEntries, atomicNumbers, label = method )

    @classmethod
    def FromUniqueEntries ( selfClass, uniqueEntries, atomicNumbers, label = None ):
        """Constructor from unique entries."""
        cdef MNDOParametersContainer self
        self = selfClass.Raw ( )
        self._CreateObject ( uniqueEntries, atomicNumbers )
        self._CheckDiatomicTerms ( )
        return self

    def MakeOrbitalBasis ( self, atomicNumbers, path = None ):
        """Make the orbital basis."""
        if path is None: path = os.path.join ( os.getenv ( "PDYNAMO3_PARAMETERS" ), "mndoParameters", "mndostong" )
        bases   = {}
        missing = set ( )
        for ( key, entry ) in self.uniqueEntries.items ( ):
            basisLabel = entry.orbitalBasisLabel
            if basisLabel is None:
                missing.add ( "atomic number {:d}".format ( key ) )
            else:
                basisPath = os.path.join ( path, basisLabel + ".yaml" )
                if os.path.exists ( basisPath ):
                    basis              = YAMLMappingFile_ToObject ( basisPath, GaussianBasis )
                    bases[key]         = basis # . No normalization.
                    basis.atomicNumber = key
                else: missing.add ( basisLabel )
        if len ( missing ) > 0:
            raise QCModelError ( "There are missing MNDO basis sets: {:s}.".format ( ", ".join ( sorted ( missing ) ) ) )
        orbitalBases = GaussianBasisContainer.FromUniqueEntries ( bases, atomicNumbers, label = "MNDO Basis" )
        self.ScaleBasesExponents ( orbitalBases ) # . With normalization.
        return orbitalBases

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    def ScaleBasesExponents ( self, bases ):
        """Scale basis container exponents."""
        # . Separated off for use separately after, for example, pickling.
        for ( key, entry ) in self.uniqueEntries.items ( ):
            basis = bases.uniqueEntries[key]
            for ( i, zeta ) in enumerate ( entry.shellExponents ):
                basis.ScaleShellExponents ( i, zeta )
            entry.DetermineNormalization ( basis )

    @property
    def atomicNumbers ( self ):
        cdef CInteger  i
        values = []
        if self.cObject != NULL:
            for i from 0 <= i < self.cObject.capacity:
                values.append ( self.cObject.entries[i].atomicNumber )
        return values

    @property
    def energyBaseLine ( self ):
        """The energy base line."""
        cdef CInteger  i
        cdef CReal     value = 0.0
        if self.cObject != NULL:
            for i from 0 <= i < self.cObject.capacity:
                value += ( self.cObject.entries[i].eheat - self.cObject.entries[i].eisol )
        return value

    @property
    def coreCharges ( self ):
        cdef RealArray1D charges = None
        cdef CInteger    i
        if self.cObject != NULL:
            charges = Array.WithExtent ( self.cObject.capacity )
            for i from 0 <= i < self.cObject.capacity:
                charges[i] = self.cObject.entries[i].zcore
        return charges
