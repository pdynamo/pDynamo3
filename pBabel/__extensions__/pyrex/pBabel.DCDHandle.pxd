from pCore.CPrimitiveTypes                   cimport CBoolean, CFalse, CInteger, CReal, CTrue
from pCore.Selection                         cimport CSelection
from pScientific.Arrays.RealArray2D          cimport CRealArray2D
from pScientific.Geometry3.Coordinates3      cimport Coordinates3
from pScientific.Symmetry.SymmetryParameters cimport CSymmetryParameters

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "fastio.h":
    pass

cdef extern from "DCDHandle.h":

    ctypedef enum CDCDStatus "DCDStatus":
        CDCDStatus_Normal                  "DCDStatus_Normal"
        CDCDStatus_AtomNumberMismatch      "DCDStatus_AtomNumberMismatch"
        CDCDStatus_BadFormat               "DCDStatus_BadFormat"
        CDCDStatus_BadRead                 "DCDStatus_BadRead"
        CDCDStatus_BadSeek                 "DCDStatus_BadSeek"
        CDCDStatus_BadWrite                "DCDStatus_BadWrite"
        CDCDStatus_FileAccessFailure       "DCDStatus_FileAccessFailure"
        CDCDStatus_InvalidDataObject       "DCDStatus_InvalidDataObject"
        CDCDStatus_InvalidFrameIndex       "DCDStatus_InvalidFrameIndex"
        CDCDStatus_OpenFailed              "DCDStatus_OpenFailed"
        CDCDStatus_OutOfMemory             "DCDStatus_OutOfMemory"

    ctypedef struct CDCDHandle "DCDHandle":
        CBoolean hasUnitCell
        CBoolean isXPLOR
        CBoolean useVelocityHeader
        CInteger currentFrame
        CInteger numberOfAtomIndices
        CInteger numberOfAtoms
        CInteger numberOfFrames
        CInteger saveFrequency
        CInteger startingFrame
        CReal    timeStep

    cdef CDCDHandle *DCDHandle_Allocate              ( )
    cdef CDCDStatus  DCDHandle_AllocateQW            ( CDCDHandle  *self )
    cdef CDCDStatus  DCDHandle_CheckNumberOfAtoms    ( CDCDHandle  *self, CInteger numberOfAtoms )
    cdef CInteger    DCDHandle_CurrentFrame          ( CDCDHandle  *self )
    cdef void        DCDHandle_Deallocate            ( CDCDHandle **self )
    cdef CInteger    DCDHandle_NumberOfFrames        ( CDCDHandle  *self )
    cdef CDCDStatus  DCDHandle_SetAtomIndices        ( CDCDHandle  *self, CSelection *selection )
    cdef CDCDStatus  DCDHandle_SetData3              ( CDCDHandle  *self, CRealArray2D *data3 )
    cdef CDCDStatus  DCDHandle_SetSymmetryParameters ( CDCDHandle  *self, CSymmetryParameters *symmetryParameters )

# . Functions.
cdef DCDStatus_Check ( CDCDStatus status )
