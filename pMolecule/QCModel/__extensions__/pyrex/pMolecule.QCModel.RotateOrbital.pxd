from pCore.CPrimitiveTypes              cimport CBoolean         , \
                                                CFalse           , \
                                                CTrue            , \
                                                CInteger         , \
                                                CReal
from pCore.Status                       cimport CStatus          , \
                                                CStatus_OK
from pScientific.Arrays.IntegerArray1D  cimport IntegerArray1D   , \
                                                CIntegerArray1D
from pScientific.Arrays.RealArray1D     cimport CRealArray1D     , \
                                                RealArray1D
from pScientific.Arrays.RealArray2D     cimport CRealArray2D     , \
                                                RealArray2D
from pScientific.Geometry3.Matrix33     cimport Matrix33

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "RotateOrbital.h":

    cdef void CRotateOrbital                "RotateOrbital"                ( CIntegerArray1D  *orbitalBasisIndices ,
                                                                             CRealArray2D     *rotation            ,
                                                                             CIntegerArray1D  *mapping             ,
                                                                             CRealArray1D     *inOrbital           ,
                                                                             CRealArray1D     *outOrbital          )
    cdef void CRotateOrbital_MakeLRotations "RotateOrbital_MakeLRotations" ( CInteger          L                   ,
                                                                             CRealArray2D     *R                   ,
                                                                             CRealArray2D     *T                   ,
                                                                             CStatus          *status              )
