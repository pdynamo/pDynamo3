from pCore.CPrimitiveTypes cimport CBoolean, CFalse, CTrue

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "Status.h":

    ctypedef enum CStatus "Status":
        CStatus_OK                    "Status_OK"
        CStatus_AlgorithmError        "Status_AlgorithmError"
        CStatus_IndexOutOfRange       "Status_IndexOutOfRange"
        CStatus_InvalidArgument       "Status_InvalidArgument"
        CStatus_InvalidArrayOperation "Status_InvalidArrayOperation"
        CStatus_MathError             "Status_MathError"
        CStatus_NonConformableArrays  "Status_NonConformableArrays"
        CStatus_OutOfMemory           "Status_OutOfMemory"
