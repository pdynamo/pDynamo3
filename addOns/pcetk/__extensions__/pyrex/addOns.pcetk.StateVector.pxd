from pCore.CPrimitiveTypes                           cimport CBoolean                , \
                                                             CFalse                  , \
                                                             CInteger                , \
                                                             CReal                   , \
                                                             CTrue                      
from pCore.Status                                    cimport CStatus                 , \
                                                             CStatus_OK              , \
                                                             CStatus_IndexOutOfRange
from pScientific.RandomNumbers.RandomNumberGenerator cimport CRandomNumberGenerator

# Include StateVector.h in the generated C code
cdef extern from "StateVector.h":
    ctypedef struct CTitrSite "TitrSite":
        CInteger  indexActive
        CInteger  indexFirst
        CInteger  indexLast
        CInteger  indexSite

    ctypedef struct CPairSite "PairSite":
        CTitrSite *a
        CTitrSite *b

    ctypedef struct CStateVector "StateVector":
        CTitrSite   *sites
        CTitrSite  **substateSites
        CInteger      nsites
        CInteger      nssites
        CPairSite   *pairs
        CInteger      npairs


    # Allocation and deallocation
    cdef CStateVector *StateVector_Allocate          (CInteger nsites, CStatus *status)
    cdef void          StateVector_AllocateSubstate  (CStateVector *self, CInteger nssites, CStatus *status)
    cdef void          StateVector_AllocatePairs     (CStateVector *self, CInteger npairs, CStatus *status)
    cdef void          StateVector_Deallocate        (CStateVector **self)

    # Copying and cloning
    cdef CStateVector *StateVector_Clone             (CStateVector *self, CStatus *status)
    cdef void          StateVector_CopyTo            (CStateVector *self, CStateVector *other, CStatus *status)

    # Setting all items at once
    cdef void          StateVector_Reset             (CStateVector *self)
    cdef void          StateVector_ResetSubstate     (CStateVector *self)
    cdef void          StateVector_ResetToMaximum    (CStateVector *self)
    cdef void          StateVector_Randomize         (CStateVector *self, CRandomNumberGenerator *generator)

    # Accessing items
    cdef void          StateVector_SetSite           (CStateVector *self, CInteger indexSite, CInteger indexFirst, CInteger indexLast, CStatus *status)
    cdef void          StateVector_SetPair           (CStateVector *self, CInteger indexPair, CInteger indexFirstSite, CInteger indexSecondSite, CReal Wmax, CStatus *status)
    cdef void          StateVector_GetPair           (CStateVector *self, CInteger indexPair, CInteger *indexFirstSite, CInteger *indexSecondSite, CReal *Wmax, CStatus *status)
    cdef CBoolean       StateVector_IsSubstate        (CStateVector *self, CInteger siteIndex, CStatus *status)
    cdef CInteger       StateVector_GetItem           (CStateVector *self, CInteger siteIndex, CStatus *status)
    cdef void          StateVector_SetItem           (CStateVector *self, CInteger siteIndex, CInteger instanceLocalIndex, CStatus *status)
    cdef CInteger       StateVector_GetActualItem     (CStateVector *self, CInteger siteIndex, CStatus *status)
    cdef void          StateVector_SetActualItem     (CStateVector *self, CInteger siteIndex, CInteger instanceGlobalIndex, CStatus *status)
    cdef CInteger       StateVector_GetSubstateItem   (CStateVector *self, CInteger index, CStatus *status)
    cdef void          StateVector_SetSubstateItem   (CStateVector *self, CInteger selectedSiteIndex, CInteger index, CStatus *status)

    # Incrementation
    cdef CBoolean       StateVector_Increment         (CStateVector *self)
    cdef CBoolean       StateVector_IncrementSubstate (CStateVector *self)

#-------------------------------------------------------------------------------
cdef class StateVector:
    cdef CStateVector   *cObject
    cdef public object   isOwner
    cdef public object   owningModel
