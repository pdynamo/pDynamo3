# . Primitive C types used by the library.

from libc.stdint cimport int16_t  , \
                         int32_t  , \
                         int64_t  , \
                         uint16_t , \
                         uint32_t , \
                         uint64_t

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "Boolean.h":

    ctypedef enum CBoolean "Boolean":
        CFalse "False" = 0
        CTrue  "True"  = 1

cdef extern from "Cardinal.h":

    ctypedef unsigned int CCardinal "Cardinal"

cdef extern from "Integer.h":

    ctypedef int CInteger "Integer"

cdef extern from "MachineTypes.h":

    ctypedef uint16_t CCardinal16 "Cardinal16"
    ctypedef uint32_t CCardinal32 "Cardinal32"
    ctypedef uint64_t CCardinal64 "Cardinal64"
    ctypedef char     CCharacter  "Character"
    ctypedef  int16_t CInteger16  "Integer16"
    ctypedef  int32_t CInteger32  "Integer32"
    ctypedef  int64_t CInteger64  "Integer64"
    ctypedef float    CReal32     "Real32"
    ctypedef double   CReal64     "Real64"

cdef extern from "Real.h":

    ctypedef double CReal "Real"

cdef extern from "Size.h":

    ctypedef size_t CSize "Size"
