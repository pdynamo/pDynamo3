# . Primitive C types used by the library.

#from libc.stdint cimport uint32_t

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

cdef extern from "Real.h":

    ctypedef double CReal "Real"

cdef extern from "Size.h":

    ctypedef size_t CSize "Size"
