"""A container for Gaussian basis sets."""

import os, os.path

from   pCore              import Clone                    , \
                                 DataType                 , \
                                 RawObjectConstructor     , \
                                 YAMLMappingFile_ToObject , \
                                 YAMLPickleFileExtension
from   pScientific        import PeriodicTable
from   pScientific.Arrays import Array
from  .GaussianBasisError import GaussianBasisError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class GaussianBasisContainer:

    def __copy__ ( self ):
        """Copying."""
        cdef GaussianBasisContainer new
        cdef CStatus cStatus = CStatus_OK
        new               = self.__class__.Raw ( )
        new.cObject       = GaussianBasisContainer_Clone ( self.cObject, &cStatus )
        if cStatus != CStatus_OK: raise GaussianBasisError ( "Error cloning container." )
        new.isOwner       = True
        new.label         = self.label
        new.uniqueEntries = Clone ( self.uniqueEntries )
        new._MakeFunctionData ( )
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
        state = { "Atomic Numbers" : self.atomicNumbers ,
                  "Unique Entries" : self.uniqueEntries }
        if self.label is not None: state["Label"] = self.label
        return state

    def __init__ ( self, capacity ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( capacity )

    def __len__ ( self ):
        """The number of basis functions."""
        return self._numberOfFunctions

    def __reduce_ex__ ( self, protocol ):
        """Pickling protocol."""
        return ( RawObjectConstructor, ( self.__class__, ), self.__getstate__ ( ) )

    def __setstate__ ( self, state ):
        """Set the state."""
        self._CreateObject ( state["Unique Entries"], state["Atomic Numbers"] )
        self.label = state.get ( "Label", None )
        self._MakeFunctionData ( )

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
        self._centerFunctionPointers = None
        self._functionCenters        = None
        self._functionLabels         = None
        self._numberOfFunctions      = -1
        self._numberOfWorkFunctions  = -1

    def _MakeFunctionData ( self ):
        """Make some function data."""
        cdef IntegerArray1D centerFunctionPointers
        cdef IntegerArray1D functionCenters
        cdef CInteger       i, n
        cdef CStatus        cStatus = CStatus_OK
        if self.cObject != NULL:
            # . Number of functions.
            self._numberOfFunctions     = GaussianBasisContainer_NumberOfFunctions     ( self.cObject )
            self._numberOfWorkFunctions = GaussianBasisContainer_NumberOfWorkFunctions ( self.cObject )
            # . Index arrays.
            centerFunctionPointers  = Array.WithExtent ( self.cObject.capacity+1, dataType = DataType.Integer )
            functionCenters         = Array.WithExtent ( self._numberOfFunctions, dataType = DataType.Integer )
            GaussianBasisContainer_MakeIndexArrays ( self.cObject, centerFunctionPointers.cObject, functionCenters.cObject, &cStatus )
            if cStatus != CStatus_OK: raise GaussianBasisError ( "Error creating index arrays." )
            self._centerFunctionPointers = centerFunctionPointers
            self._functionCenters        = functionCenters
            # . Function labels.
            labels = []
            if self.cObject != NULL:
                uniqueLabels = { atomicNumber : basis.MakeFunctionLabels ( ) for ( atomicNumber, basis ) in self.uniqueEntries.items ( ) }
                for i from 0 <= i < self.cObject.capacity:
                    labels.extend ( uniqueLabels[self.cObject.entries[i].atomicNumber] )
            self._functionLabels = labels

    @classmethod
    def FromParameterDirectory ( selfClass            ,  
                                 basisLabel           ,  
                                 atomicNumbers        ,  
                                 basisType     = None ,
                                 path          = None ):
        """Constructor given a basis label and a list of atomic numbers."""
        if path is None: path = os.path.join ( os.getenv ( "PDYNAMO3_PARAMETERS" ), "gaussianBasisSets", basisLabel )
        missing       = set ( )
        uniqueEntries = {}
        for atomicNumber in sorted ( set ( atomicNumbers ) ):
            entryPath = os.path.join ( path, PeriodicTable.Symbol ( atomicNumber ) + YAMLPickleFileExtension )
            #try:
            basis = YAMLMappingFile_ToObject ( entryPath, GaussianBasis )
            basis.Finalize ( )
            uniqueEntries[atomicNumber] = basis
            #except: missing.add ( "{:d}".format ( atomicNumber ) )
        if len ( missing ) > 0:
            raise GaussianBasisError ( "There are missing {:s} basis sets: {:s}.".format ( basisLabel, ", ".join ( sorted ( missing ) ) ) )
        else:
            return selfClass.FromUniqueEntries ( uniqueEntries, atomicNumbers, label = basisLabel )

    @classmethod
    def FromUniqueEntries ( selfClass, uniqueEntries, atomicNumbers, label = None ):
        """Constructor from unique entries."""
        cdef GaussianBasisContainer self
        self = selfClass.Raw ( )
        self._CreateObject ( uniqueEntries, atomicNumbers )
        if label is not None: self.label = label
        self._MakeFunctionData ( )
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
        cdef CStatus       cStatus = CStatus_OK
        L  = max ( [ entry.maximumAngularMomentum for entry in self.uniqueEntries.values ( ) ] )
        d  = ( ( L + 1 ) * ( L + 2 ) * ( L + 3 ) ) // 6
        Tc = Array.WithShape ( [ d, d ] )
        GaussianBasis_MakeLRotations ( L, R.cObject, Tc.cObject, &cStatus )
        Tn = {}
        for ( n, entry ) in self.uniqueEntries.items ( ):
            d = entry.numberOfFunctions
            T = Array.WithShape ( [ d, d ] )
            GaussianBasis_MakeRotationMatrix ( entry.cObject, Tc.cObject, T.cObject, &cStatus )
            Tn[n] = T
        if cStatus != CStatus_OK: raise GaussianBasisError ( "Error making rotation matrices." )
        return [ Tn[n] for n in self.atomicNumbers ]

    @property
    def atomicNumbers ( self ):
        cdef CInteger i
        values = []
        if self.cObject != NULL:
            for i from 0 <= i < self.cObject.capacity:
                values.append ( self.cObject.entries[i].atomicNumber )
        return values

    # . For debugging only.
#   @property
#   def c2s ( self ):
#       cdef CStatus     cStatus = CStatus_OK
#       cdef RealArray2D T
#       T = Array.WithShape ( [ self.numberOfWorkFunctions, self.numberOfFunctions ] )
#       GaussianBasisContainer_MakeC2S ( self.cObject, T.cObject, &cStatus )
#       if cStatus != CStatus_OK: raise GaussianBasisError ( "Error making C->S transformation." )
#       return T

    @property
    def centerFunctionPointers ( self ):
        return self._centerFunctionPointers

    @property
    def functionCenters ( self ):
        return self._functionCenters

    @property
    def functionLabels ( self ):
        return self._functionLabels

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
        return self._numberOfFunctions

    @property
    def numberOfWorkFunctions ( self ):
        return self._numberOfWorkFunctions
