from pCore.CPrimitiveTypes              cimport CBoolean         , \
                                                CFalse           , \
                                                CTrue            , \
                                                CInteger         , \
                                                CReal              
from pCore.Status                       cimport CStatus          , \
                                                CStatus_OK
from pScientific.Arrays.RealArray1D     cimport CRealArray1D     , \
                                                RealArray1D
from pScientific.Arrays.RealArray2D     cimport CRealArray2D     , \
                                                RealArray2D
from pScientific.Arrays.SymmetricMatrix cimport CSymmetricMatrix , \
                                                SymmetricMatrix        

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "OrthogonalizingTransformation.h":

    cdef CInteger COrthogonalizingTransformation "OrthogonalizingTransformation" ( CSymmetricMatrix *S                   ,
                                                                                   CBoolean          doCanonical         ,
                                                                                   CBoolean          preserveInput       ,
                                                                                   CReal            *eigenValueTolerance ,
                                                                                   CRealArray1D     *eigenValues         ,
                                                                                   CRealArray2D     *eigenVectors        ,
                                                                                   CRealArray2D     *X                   ,
                                                                                   CStatus          *status              )
