from pCore.CPrimitiveTypes             cimport CBoolean        , \
                                               CFalse          , \
                                               CTrue           , \
                                               CInteger        , \
                                               CReal
from pCore.Status                      cimport CStatus         , \
                                               CStatus_OK
from pScientific.Arrays.IntegerArray1D cimport CIntegerArray1D , \
                                               IntegerArray1D        

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "DFTFunctionalModel.h":

    ctypedef struct CDFTFunctionalModel "DFTFunctionalModel":
        CBoolean isSpinRestricted
        CInteger order

    cdef CDFTFunctionalModel *DFTFunctionalModel_Allocate        ( CInteger              numberOfFunctionals ,
                                                                   CStatus              *status              )
    cdef CDFTFunctionalModel *DFTFunctionalModel_Clone           ( CDFTFunctionalModel  *self                ,
                                                                   CStatus              *status              )
    cdef void                 DFTFunctionalModel_Deallocate      ( CDFTFunctionalModel **self                )
    cdef CReal                DFTFunctionalModel_ExchangeScaling ( CDFTFunctionalModel  *self                )
    cdef CDFTFunctionalModel *DFTFunctionalModel_MakeFromIDs     ( CIntegerArray1D      *ids                 ,
                                                                   CBoolean              isSpinRestricted    ,
                                                                   CStatus              *status              )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class DFTFunctionalModel:

    cdef CDFTFunctionalModel *cObject
    cdef IntegerArray1D       ids
    cdef public object        isOwner
