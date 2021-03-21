"""A DFT grid."""

from  pCore          import Clone                , \
                            RawObjectConstructor
from  pScientific    import Magnitude_Adjust
from .DFTDefinitions import DFTGridAccuracy
from .QCModelError   import QCModelError

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Accuracies.
DFTGridAccuracy_FromCEnum = { DFTGridAccuracy_VeryLow  : DFTGridAccuracy.VeryLow  ,
                              DFTGridAccuracy_Low      : DFTGridAccuracy.Low      ,
                              DFTGridAccuracy_Medium   : DFTGridAccuracy.Medium   ,
                              DFTGridAccuracy_High     : DFTGridAccuracy.High     ,
                              DFTGridAccuracy_VeryHigh : DFTGridAccuracy.VeryHigh }
DFTGridAccuracy_ToCEnum   = { DFTGridAccuracy.VeryLow  : DFTGridAccuracy_VeryLow  ,
                              DFTGridAccuracy.Low      : DFTGridAccuracy_Low      ,
                              DFTGridAccuracy.Medium   : DFTGridAccuracy_Medium   ,
                              DFTGridAccuracy.High     : DFTGridAccuracy_High     ,
                              DFTGridAccuracy.VeryHigh : DFTGridAccuracy_VeryHigh }

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class DFTGrid:

    def __dealloc__ ( self ):
        """Destructor."""
        cdef CStatus  cStatus = CStatus_OK
        if self.isOwner:
            DFTGrid_Deallocate ( &self.cObject, &cStatus )
            self.isOwner = False
            if cStatus  != CStatus_OK: raise QCModelError ( "Error deallocating DFT grid." )
 
    def __init__ ( self, accuracy ):
        """Constructor."""
        self._Initialize ( )
        self._Allocate ( accuracy )

    def _Allocate ( self, accuracy ):
        """Constructor."""
        cdef CDFTGridAccuracy cAccuracy
        cdef CStatus          cStatus = CStatus_OK
        cAccuracy    = DFTGridAccuracy_ToCEnum[accuracy]
        self.cObject = DFTGrid_Allocate ( cAccuracy, &cStatus )
        self.isOwner = True
        if cStatus != CStatus_OK: raise QCModelError ( "Error allocating DFT grid." )

    def _Initialize ( self ):
        """Initialization."""
        self.cObject = NULL
        self.isOwner = False

    @classmethod
    def Construct (                selfClass                ,
                                   accuracy                 ,
                    IntegerArray1D atomicNumbers   not None , 
                    Coordinates3   qcCoordinates3  not None ):
        """Constructor of a grid of a given accuracy from atomic numbers and coordinates."""
        cdef CDFTGridAccuracy cAccuracy
        cdef CStatus          cStatus = CStatus_OK
        cdef DFTGrid          self
        cAccuracy    = DFTGridAccuracy_ToCEnum[accuracy]
        self         = selfClass.Raw ( )
        self.cObject = DFTGrid_Construct ( cAccuracy               ,
                                           atomicNumbers.cObject   ,
                                           qcCoordinates3.cObject  ,
                                           &cStatus                )
        self.isOwner = True
        if cStatus != CStatus_OK: raise QCModelError ( "Error constructing DFT grid." )
        return self

    @staticmethod
    def EstimatedPoints (                accuracy               ,
                          IntegerArray1D atomicNumbers not None ):
        """Estimate the number of points in the grid."""
        cdef CDFTGridAccuracy cAccuracy
        cdef CInteger         p
        cdef CStatus          cStatus = CStatus_OK
        cAccuracy = DFTGridAccuracy_ToCEnum[accuracy]
        p = DFTGrid_EstimatedPoints ( cAccuracy, atomicNumbers.cObject, &cStatus )
        if cStatus != CStatus_OK: raise QCModelError ( "Error estimating number of grid points." )
        return p

    @classmethod
    def Raw ( selfClass ):
        """Raw constructor."""
        self = selfClass.__new__ ( selfClass )
        self._Initialize ( )
        return self

    @property
    def functionByteSize ( self ):
        cdef CReal    result
        cdef CStatus  cStatus = CStatus_OK
        result = DFTGrid_FunctionByteSize ( self.cObject, &cStatus )
        if cStatus != CStatus_OK: raise QCModelError ( "Error determining DFT grid function data size." )
        return Magnitude_Adjust ( result )

    @property
    def hasFunctionData ( self ):
        cdef CBoolean result
        cdef CStatus  cStatus = CStatus_OK
        result = DFTGrid_HasFunctionData ( self.cObject, &cStatus )
        if cStatus != CStatus_OK: raise QCModelError ( "Error determining whether DFT grid has function data." )
        return ( result == CTrue )

    @property
    def numberOfFunctionValues ( self ):
        cdef CInteger result
        cdef CStatus  cStatus = CStatus_OK
        result = DFTGrid_NumberOfFunctionValues ( self.cObject, &cStatus )
        if cStatus != CStatus_OK: raise QCModelError ( "Error determining DFT grid number of function values." )
        return result

    @property
    def numberOfPoints ( self ):
        return DFTGrid_NumberOfPoints ( self.cObject )
