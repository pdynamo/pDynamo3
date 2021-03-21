from addOns.pcetk.EnergyModel                        cimport CEnergyModel          , \
                                                             EnergyModel
from addOns.pcetk.StateVector                        cimport CStateVector          , \
                                                             StateVector           , \
                                                             StateVector_GetPair   , \
                                                             StateVector_Randomize
from pCore.CPrimitiveTypes                           cimport CBoolean              , \
                                                             CFalse                , \
                                                             CInteger              , \
                                                             CReal                 , \
                                                             CTrue                    
from pCore.Status                                    cimport CStatus               , \
                                                             CStatus_OK
from pScientific.Arrays.RealArray1D                  cimport RealArray1D_Scale     , \
                                                             RealArray1D_Set
from pScientific.RandomNumbers.RandomNumberGenerator cimport CRandomNumberGenerator as CGenerator

cdef extern from "MCModelDefault.h":
    ctypedef struct CMCModelDefault "MCModelDefault":
        CReal           limit
        CInteger        nprod
        CInteger        nequil
        CEnergyModel  *energyModel
        CStateVector  *vector
        CGenerator    *generator


    cdef CMCModelDefault *MCModelDefault_Allocate               (CReal limit, CInteger nequil, CInteger nprod, CInteger randomSeed, CStatus *status)
    cdef void             MCModelDefault_Deallocate             (CMCModelDefault **self)
    cdef void             MCModelDefault_LinkToEnergyModel      (CMCModelDefault *self, CEnergyModel *energyModel, CStatus *status)
    cdef CReal             MCModelDefault_MCScan                 (CMCModelDefault *self, CReal pH, CInteger nmoves, CInteger *movesDone, CInteger *movesAccepted, CInteger *flipsDone, CInteger *flipsAccepted)
    cdef CInteger          MCModelDefault_FindPairs              (CMCModelDefault *self, CInteger npairs, CStatus *status)
    cdef void             MCModelDefault_UpdateProbabilities    (CMCModelDefault *self)
    cdef void             MCModelDefault_Equilibration          (CMCModelDefault *self, CReal pH)
    cdef void             MCModelDefault_Production             (CMCModelDefault *self, CReal pH)

#-------------------------------------------------------------------------------
cdef class MCModelDefault:
    cdef CMCModelDefault *cObject
    cdef public object    isOwner
    cdef public object    owningModel
