from pCore.CPrimitiveTypes              cimport CBoolean         , \
                                                CFalse           , \
                                                CInteger         , \
                                                CReal            , \
                                                CTrue
from pCore.Status                       cimport CStatus          , \
                                                CStatus_OK
from pScientific.Arrays.IntegerArray1D  cimport IntegerArray1D   , \
                                                CIntegerArray1D  
from pScientific.Arrays.RealArray1D     cimport CRealArray1D     , \
                                                RealArray1D
from pScientific.Arrays.SymmetricMatrix cimport CSymmetricMatrix , \
                                                SymmetricMatrix

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "Mulliken.h":

    cdef void Mulliken_AtomicCharges            ( CIntegerArray1D  *basisIndices ,
                                                  CSymmetricMatrix *density      ,
                                                  CSymmetricMatrix *overlap      ,
                                                  CRealArray1D     *charges      )
    cdef void Mulliken_BondOrders               ( CIntegerArray1D  *basisIndices ,
                                                  CSymmetricMatrix *density      ,
                                                  CSymmetricMatrix *overlap      ,
                                                  CSymmetricMatrix *bondOrders   ,
                                                  CStatus          *status       )
    cdef void Mulliken_ChargeDensityDerivatives ( CIntegerArray1D  *basisIndices ,
                                                  CRealArray1D     *potentials   ,
                                                  CSymmetricMatrix *overlap      ,
                                                  CSymmetricMatrix *fock         )
    cdef void Mulliken_WeightedDensity          ( CIntegerArray1D  *basisIndices ,
                                                  CRealArray1D     *potentials   ,
                                                  CSymmetricMatrix *density      ,
                                                  CSymmetricMatrix *wDensity     )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class MullikenMultipoleEvaluator:
    pass
