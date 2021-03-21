from addOns.pcetk.StateVector                        cimport CStateVector        , \
                                                             StateVector         , \
                                                             StateVector_SetSite
from pCore.CPrimitiveTypes                           cimport CBoolean            , \
                                                             CFalse              , \
                                                             CInteger            , \
                                                             CReal               , \
                                                             CTrue                  
from pCore.Status                                    cimport CStatus             , \
                                                             CStatus_OK
from pScientific.Arrays.RealArray1D                  cimport CRealArray1D        , \
                                                             RealArray1D

# Include EnergyModel.h in the generated C code
cdef extern from "EnergyModel.h":
    ctypedef struct CEnergyModel "EnergyModel":
        CInteger        nstates
        CInteger        ninstances
        CReal           temperature
        CStateVector  *vector
        CRealArray1D  *probabilities

    # Allocation and deallocation
    cdef CEnergyModel *EnergyModel_Allocate                          (CInteger nsites, CInteger ninstances, CStatus *status)
    cdef void          EnergyModel_Deallocate                        (CEnergyModel **self)
    # Handling of the interactions matrix
    cdef CBoolean       EnergyModel_CheckInteractionsSymmetric        (CEnergyModel *self, CReal tolerance, CReal *maxDeviation)
    cdef void          EnergyModel_SymmetrizeInteractions            (CEnergyModel *self, CStatus *status)
    cdef void          EnergyModel_ResetInteractions                 (CEnergyModel *self)
    cdef void          EnergyModel_ScaleInteractions                 (CEnergyModel *self, CReal scale)
    # Calculation of energies
    cdef CReal          EnergyModel_CalculateMicrostateEnergy         (CEnergyModel *self, CStateVector *vector, CReal pH)
    cdef CReal          EnergyModel_CalculateMicrostateEnergyUnfolded (CEnergyModel *self, CStateVector *vector, CReal pH)
    # Calculation of partition functions
    cdef CReal          EnergyModel_CalculateZunfolded                (CEnergyModel *self, CReal pH, CReal Gzero, CStatus *status)
    cdef CReal          EnergyModel_CalculateZfolded                  (CEnergyModel *self, CReal pH, CReal Gzero, CStatus *status)
    # Functions for getting items
    cdef CReal          EnergyModel_GetGmodel                         (CEnergyModel *self, CInteger instIndexGlobal)
    cdef CReal          EnergyModel_GetGintr                          (CEnergyModel *self, CInteger instIndexGlobal)
    cdef CInteger       EnergyModel_GetProtons                        (CEnergyModel *self, CInteger instIndexGlobal)
    cdef CReal          EnergyModel_GetProbability                    (CEnergyModel *self, CInteger instIndexGlobal)
    cdef CReal          EnergyModel_GetInteraction                    (CEnergyModel *self, CInteger instIndexGlobalA, CInteger instIndexGlobalB)
    cdef CReal          EnergyModel_GetInterSymmetric                 (CEnergyModel *self, CInteger instIndexGlobalA, CInteger instIndexGlobalB)
    cdef CReal          EnergyModel_GetDeviation                      (CEnergyModel *self, CInteger instIndexGlobalA, CInteger instIndexGlobalB)
    # Functions for setting items
    cdef void          EnergyModel_SetGmodel                         (CEnergyModel *self, CInteger instIndexGlobal, CReal value)
    cdef void          EnergyModel_SetGintr                          (CEnergyModel *self, CInteger instIndexGlobal, CReal value)
    cdef void          EnergyModel_SetProtons                        (CEnergyModel *self, CInteger instIndexGlobal, CInteger value)
    cdef void          EnergyModel_SetProbability                    (CEnergyModel *self, CInteger instIndexGlobal, CReal value)
    cdef void          EnergyModel_SetInteraction                    (CEnergyModel *self, CInteger instIndexGlobalA, CInteger instIndexGlobalB, CReal value)

    # Calculation of probabilities
    cdef void EnergyModel_CalculateProbabilitiesAnalytically (CEnergyModel *self, CReal pH, CStatus *status)
    cdef void EnergyModel_CalculateProbabilitiesAnalyticallyUnfolded (CEnergyModel *self, CReal pH, CStatus *status)


#-------------------------------------------------------------------------------
cdef class EnergyModel:
    cdef CEnergyModel  *cObject
    cdef public object  isOwner
    cdef public object  owningModel
