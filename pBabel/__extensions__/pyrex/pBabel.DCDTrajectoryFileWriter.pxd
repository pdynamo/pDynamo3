from pBabel.DCDHandle                        cimport CDCDHandle                      , \
                                                     CDCDStatus                      , \
                                                     CDCDStatus_OutOfMemory          , \
                                                     DCDHandle_Allocate              , \
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
cdef extern from "DCDWrite.h":

    cdef void       DCDWrite_Close  ( CDCDHandle **self )
    cdef CDCDStatus DCDWrite_Frame  ( CDCDHandle  *self )
    cdef void       DCDWrite_Header ( CDCDHandle  *self, char *title )
    cdef CDCDStatus DCDWrite_Open   ( CDCDHandle  *self, char *path  )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class DCDTrajectoryFileWriter:

    cdef CDCDHandle    *cObject
    cdef public object  isOpen
    cdef public object  isTrajectory
    cdef public object  owner
    cdef public object  path
    cdef public object  title
