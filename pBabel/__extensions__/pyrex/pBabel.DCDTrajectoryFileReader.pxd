from pBabel.DCDHandle                        cimport CDCDHandle                      , \
                                                     CDCDStatus                      , \
                                                     CDCDStatus_OutOfMemory          , \
                                                     DCDHandle_Allocate              , \
                                                     DCDHandle_AllocateQW            , \
                                                     DCDHandle_CheckNumberOfAtoms    , \
                                                     DCDHandle_CurrentFrame          , \
                                                     DCDHandle_NumberOfFrames        , \
                                                     DCDHandle_SetAtomIndices        , \
                                                     DCDHandle_SetData3              , \
                                                     DCDHandle_SetSymmetryParameters , \
                                                     DCDStatus_Check
from pCore.CPrimitiveTypes                   cimport CBoolean, CFalse, CInteger, CReal, CTrue
from pCore.Selection                         cimport Selection
from pScientific.Geometry3.Coordinates3      cimport Coordinates3
from pScientific.Symmetry.SymmetryParameters cimport SymmetryParameters

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "DCDRead.h":

    cdef void       DCDRead_Close     ( CDCDHandle **self )
    cdef CDCDStatus DCDRead_Frame     ( CDCDHandle  *self )
    cdef CDCDStatus DCDRead_GotoFrame ( CDCDHandle  *self, CInteger f )
    cdef CDCDStatus DCDRead_Header    ( CDCDHandle  *self )
    cdef CDCDStatus DCDRead_Open      ( CDCDHandle  *self, char *path  )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class DCDTrajectoryFileReader:

    cdef CDCDHandle    *cObject
    cdef public object  isOpen
    cdef public object  isTrajectory
    cdef public object  owner
    cdef public object  path
