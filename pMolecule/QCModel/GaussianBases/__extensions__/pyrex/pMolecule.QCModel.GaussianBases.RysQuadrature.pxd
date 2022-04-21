from pCore.CPrimitiveTypes cimport CBoolean   , \
                                   CFalse     , \
                                   CTrue      , \
                                   CInteger   , \
                                   CReal
from pCore.Status          cimport CStatus    , \
                                   CStatus_OK

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "RysQuadrature.h":

    ctypedef struct CRysQuadrature "RysQuadrature":
        CReal *roots
        CReal *weights

    cdef CInteger CRysQuadrature_MaximumRoots "RysQuadrature_MaximumRoots" ( )
    cdef void     CRysQuadrature_Roots        "RysQuadrature_Roots"        ( CRysQuadrature *roots, CInteger nRoots, CReal x )
