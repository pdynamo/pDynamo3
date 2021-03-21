from pCore.CPrimitiveTypes                     cimport CBoolean                    , \
                                                       CFalse                      , \
                                                       CInteger                    , \
                                                       CReal                       , \
                                                       CTrue
from pCore.Status                              cimport CStatus                     , \
                                                       CStatus_OK
from pMolecule.QCModel.MNDOParametersContainer cimport CMNDOParametersContainer    , \
                                                       MNDOParametersContainer
from pScientific.Arrays.IntegerArray1D         cimport IntegerArray1D              , \
                                                       CIntegerArray1D      
from pScientific.Arrays.RealArray1D            cimport CRealArray1D                , \
                                                       RealArray1D
from pScientific.Arrays.SymmetricMatrix        cimport CSymmetricMatrix            , \
                                                       SymmetricMatrix

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "MNDOMultipoles.h":

    ctypedef enum CMultipoleRepresentation "MultipoleRepresentation":
        MultipoleRepresentation_Buckingham = 1 ,
        MultipoleRepresentation_Cartesian  = 2 ,
        MultipoleRepresentation_Spherical  = 3

    cdef void MNDO_AtomicMultipoles     ( CMNDOParametersContainer *parameters              ,
                                          CIntegerArray1D          *basisIndices            ,
                                          CSymmetricMatrix         *density                 ,
                                          CMultipoleRepresentation  multipoleRepresentation ,
                                          CInteger                  multipoleOrder          ,
                                          CRealArray1D             *multipoles              )
    cdef void MNDO_AtomicMultipolesFock ( CMNDOParametersContainer *parameters              ,
                                          CIntegerArray1D          *basisIndices            ,
                                          CRealArray1D             *potentials              ,
                                          CInteger                  multipoleOrder          ,
                                          CSymmetricMatrix         *fock                    )
    cdef void MNDO_BondOrders           ( CIntegerArray1D          *basisIndices            ,
                                          CSymmetricMatrix         *density                 ,
                                          CSymmetricMatrix         *bondOrders              )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class MNDOMultipoleEvaluator:
    pass
