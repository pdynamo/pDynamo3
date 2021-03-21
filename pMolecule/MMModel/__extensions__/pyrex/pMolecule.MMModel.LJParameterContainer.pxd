from pCore.CPrimitiveTypes cimport CBoolean, CFalse, CInteger, CReal, CTrue

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "LJParameterContainer.h":

    ctypedef struct CLJParameterContainer "LJParameterContainer":
        CInteger  ntypes
        CReal    *epsilon
        CReal    *sigma
        CReal    *tableA
        CReal    *tableB
        CInteger *tableindex

    cdef CLJParameterContainer *LJParameterContainer_Allocate              ( CInteger ntypes )
    cdef CLJParameterContainer *LJParameterContainer_Clone                 ( CLJParameterContainer  *self )
    cdef void                   LJParameterContainer_Deallocate            ( CLJParameterContainer **self )
    cdef void                   LJParameterContainer_MakeEpsilonSigmaAMBER ( CLJParameterContainer  *self )
    cdef void                   LJParameterContainer_MakeEpsilonSigmaOPLS  ( CLJParameterContainer  *self )
    cdef void                   LJParameterContainer_MakeTableAMBER        ( CLJParameterContainer  *self )
    cdef void                   LJParameterContainer_MakeTableOPLS         ( CLJParameterContainer  *self )
    cdef CLJParameterContainer *LJParameterContainer_MergeEpsilonSigma     ( CLJParameterContainer  *self, CLJParameterContainer *other )
    cdef void                   LJParameterContainer_Scale                 ( CLJParameterContainer  *self, CReal scale )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class LJParameterContainer:

    cdef CLJParameterContainer *cObject
    cdef public object          analyticForm
    cdef public object          isOwner
    cdef public object          label
    cdef public object          parameterKeys
