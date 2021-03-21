from pCore.CPrimitiveTypes              cimport CBoolean     , \
                                                CFalse       , \
                                                CInteger     , \
                                                CReal        , \
                                                CTrue
from pScientific.Arrays.RealArray1D     cimport CRealArray1D , \
                                                RealArray1D
from pScientific.Arrays.RealArray2D     cimport CRealArray2D
from pScientific.Geometry3.Coordinates3 cimport Coordinates3

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "Lebedev.h":

    cdef void      LebedevLaikov_GridPointsWeights ( CInteger      numberOfPoints ,
                                                     CRealArray2D *gridPoints     ,
                                                     CRealArray1D *weights        )
    cdef CInteger  LebedevLaikov_Number_Of_Points  ( CInteger      lValue         )
