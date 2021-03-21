from pCore.CPrimitiveTypes                           cimport CBoolean            , \
                                                             CFalse              , \
                                                             CInteger            , \
                                                             CReal               , \
                                                             CTrue                  
from pScientific.Arrays.RealArray1D                  cimport CRealArray1D        , \
                                                             RealArray1D         , \
                                                             RealArray1D_SetItem 
from pScientific.RandomNumbers.RandomNumberGenerator cimport CRandomNumberGenerator , \
                                                             RandomNumberGenerator

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "RandomNumberDistribution.h":

    # . Functions.
    cdef CReal RandomNumberDistribution_GaussianBoxMueller ( CRandomNumberGenerator *rng, CReal mu, CReal sigma )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class NormalDeviateGenerator:

    cdef public RandomNumberGenerator randomNumberGenerator
    cdef public object                mu
    cdef public object                sigma
