from pCore.CPrimitiveTypes              cimport CBoolean        , \
                                                CInteger        , \
                                                CFalse          , \
                                                CTrue           , \
                                                CInteger        , \
                                                CReal              
from pScientific.Arrays.RealArray1D     cimport CRealArray1D    , \
                                                RealArray1D        
from pScientific.Arrays.RealArray2D     cimport CRealArray2D
from pScientific.Geometry3.Coordinates3 cimport Coordinates3       

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "QCDispersionDFTD2.h":

    cdef CReal CQCDispersionDFTD2_Energy "QCDispersionDFTD2_Energy" ( CReal          s6           ,
                                                                      CReal          sR           ,
                                                                      CReal          dR           ,
                                                                      CRealArray1D  *sqrtC6       ,
                                                                      CRealArray1D  *r0           ,
                                                                      CRealArray2D  *coordinates3 ,
                                                                      CRealArray2D  *gradients3   )
