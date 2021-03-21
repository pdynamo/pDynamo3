from pCore.CPrimitiveTypes          cimport CBoolean            , \
                                            CCardinal           , \
                                            CFalse              , \
                                            CInteger            , \
                                            CReal               , \
                                            CTrue
from pScientific.Arrays.RealArray1D cimport CRealArray1D        , \
                                            RealArray1D         , \
                                            RealArray1D_SetItem

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "RandomNumberGenerator.h":

    # . The random number generator type type.
    ctypedef struct CRandomNumberGeneratorType "RandomNumberGeneratorType":
        pass

    # . The random number generator type.
    ctypedef struct CRandomNumberGenerator "RandomNumberGenerator":
        pass

    # . Functions.
    cdef CRandomNumberGenerator *RandomNumberGenerator_Allocate     ( CRandomNumberGeneratorType *type )
    cdef CRandomNumberGenerator *RandomNumberGenerator_Clone        ( CRandomNumberGenerator  *self )
    cdef void                    RandomNumberGenerator_Deallocate   ( CRandomNumberGenerator **self )
    cdef CCardinal               RandomNumberGenerator_NextCardinal ( CRandomNumberGenerator  *self )
    cdef CReal                   RandomNumberGenerator_NextReal     ( CRandomNumberGenerator  *self )
    cdef CReal                   RandomNumberGenerator_NextRealOpen ( CRandomNumberGenerator  *self )
    cdef void                    RandomNumberGenerator_SetSeed      ( CRandomNumberGenerator  *self, CCardinal seed )

    # . Type declarations.
    cdef CRandomNumberGeneratorType *RandomNumberGeneratorType_MersenneTwister

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class RandomNumberGenerator:

    cdef CRandomNumberGenerator *cObject
    cdef public object           initialSeed
    cdef public object           isOwner
