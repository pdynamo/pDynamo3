from pCore.CPrimitiveTypes             cimport CBoolean        , \
                                               CFalse          , \
                                               CInteger        , \
                                               CReal           , \
                                               CTrue
from pCore.Status                      cimport CStatus         , \
                                               CStatus_OK
from pScientific.Arrays.IntegerArray2D cimport CIntegerArray2D , \
                                               IntegerArray2D
from pScientific.Arrays.RealArray1D    cimport CRealArray1D    , \
                                               RealArray1D 
from pScientific.Arrays.RealArray2D    cimport CRealArray2D    , \
                                               RealArray2D 
from pScientific.Arrays.RealArrayND    cimport CRealArrayND    , \
                                               RealArrayND 
from pScientific.Geometry3.RegularGrid cimport CRegularGrid    , \
                                               RegularGrid 

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "MarchingCubes.h":

    # . Functions.
    cdef void CMarchingCubes_Isosurface3D "MarchingCubes_Isosurface3D" ( CRegularGrid    *grid           ,
                                                                         CRealArrayND    *data           ,
                                                                         CReal            isovalue       ,
                                                                         CRealArray2D    *polygonNormals ,
                                                                         CIntegerArray2D *polygons       ,
                                                                         CRealArray2D    *vertexNormals  ,
                                                                         CRealArray2D    *vertices       ,
                                                                         CStatus         *status         )
