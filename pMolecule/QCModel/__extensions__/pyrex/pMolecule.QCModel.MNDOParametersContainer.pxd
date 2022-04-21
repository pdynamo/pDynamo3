from pCore.CPrimitiveTypes                                  cimport CBoolean               , \
                                                                    CFalse                 , \
                                                                    CInteger               , \
                                                                    CReal                  , \
                                                                    CTrue
from pCore.Status                                           cimport CStatus
from pMolecule.QCModel.GaussianBases.GaussianBasisContainer cimport GaussianBasisContainer
from pMolecule.QCModel.MNDOParameters                       cimport CMNDOParameters        , \
                                                                    MNDOParameters
from pScientific.Arrays.RealArray1D                         cimport RealArray1D

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "MNDOParametersContainer.h":

    ctypedef struct CMNDOParametersContainer "MNDOParametersContainer":
        CBoolean          isOwner 
        CInteger          capacity
        CMNDOParameters **entries 

    cdef CMNDOParametersContainer *MNDOParametersContainer_Allocate   ( CInteger                   capacity ,
                                                                        CStatus                   *status   )
    cdef CMNDOParametersContainer *MNDOParametersContainer_Clone      ( CMNDOParametersContainer  *self     ,
                                                                        CStatus                   *status   )
    cdef void                      MNDOParametersContainer_Deallocate ( CMNDOParametersContainer **self     )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class MNDOParametersContainer:

    cdef CMNDOParametersContainer *cObject
    cdef public object             isOwner
    cdef public object             label
    cdef public object             uniqueEntries
