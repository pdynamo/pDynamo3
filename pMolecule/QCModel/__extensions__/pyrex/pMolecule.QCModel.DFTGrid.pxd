from pCore.CPrimitiveTypes              cimport CBoolean        , \
                                                CFalse          , \
                                                CTrue           , \
                                                CInteger        , \
                                                CReal
from pCore.Status                       cimport CStatus         , \
                                                CStatus_OK
from pScientific.Arrays.IntegerArray1D  cimport CIntegerArray1D , \
                                                IntegerArray1D        
from pScientific.Arrays.RealArray2D     cimport CRealArray2D
from pScientific.Geometry3.Coordinates3 cimport Coordinates3

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "DFTGridWeights.h":

    ctypedef struct CDFTGridWeights "DFTGridWeights":
        pass

cdef extern from "List.h":

    ctypedef struct CList "List":
        pass

cdef extern from "DFTGrid.h":

    ctypedef enum CDFTGridAccuracy "DFTGridAccuracy":
        DFTGridAccuracy_VeryLow  = 0 ,
        DFTGridAccuracy_Low      = 1 ,
        DFTGridAccuracy_Medium   = 2 ,
        DFTGridAccuracy_High     = 3 ,
        DFTGridAccuracy_VeryHigh = 4

    ctypedef struct CDFTGridWeights "DFTGridWeights":
        pass

    ctypedef struct CDFTGridPointBlock "DFTGridPointBlock":
        pass

    ctypedef struct CDFTGrid "DFTGrid":
        CDFTGridAccuracy     accuracy
        CInteger             blockSize
        CInteger             numberOfPoints
        CInteger             numberOfRecords
        CReal                bfTolerance
        CReal                rhoTolerance
        CDFTGridPointBlock **records
        CDFTGridWeights     *weights
        CList               *points

    cdef CDFTGrid *DFTGrid_Allocate              ( CDFTGridAccuracy  accuracy       ,
                                                   CStatus          *status         )
    cdef CDFTGrid *DFTGrid_Construct             ( CDFTGridAccuracy  accuracy       ,
                                                   CIntegerArray1D  *atomicNumbers  ,
                                                   CRealArray2D     *qcCoordinates3 ,
                                                   CStatus          *status         )
    cdef void     DFTGrid_Deallocate             ( CDFTGrid        **self           ,
                                                   CStatus          *status         )
    cdef void     DFTGrid_DeallocateFunctionData ( CDFTGrid         *self           ,
                                                   CStatus          *status         )
    cdef CInteger DFTGrid_EstimatedPoints        ( CDFTGridAccuracy  accuracy       ,
                                                   CIntegerArray1D  *atomicNumbers  ,
                                                   CStatus          *status         )
    cdef CReal    DFTGrid_FunctionByteSize       ( CDFTGrid         *self           ,
                                                   CStatus          *status         )
    cdef CBoolean DFTGrid_HasFunctionData        ( CDFTGrid         *self           ,
                                                   CStatus          *status         )
    cdef CInteger DFTGrid_NumberOfFunctionValues ( CDFTGrid         *self           ,
                                                   CStatus          *status         )
    cdef CInteger DFTGrid_NumberOfPoints         ( CDFTGrid         *self           )
    cdef CInteger DFTGrid_NumberOfRecords        ( CDFTGrid         *self           )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class DFTGrid:

    cdef CDFTGrid      *cObject
    cdef public object  isOwner
