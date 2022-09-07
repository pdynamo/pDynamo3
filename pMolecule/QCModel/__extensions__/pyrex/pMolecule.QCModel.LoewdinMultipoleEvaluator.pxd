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
from pScientific.Arrays.RealArray2D     cimport CRealArray2D     , \
                                                RealArray2D
from pScientific.Arrays.SymmetricMatrix cimport CSymmetricMatrix , \
                                                SymmetricMatrix

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "Loewdin.h":

    cdef void  Loewdin_AtomicCharges                  ( CIntegerArray1D  *basisIndices        ,
                                                        CSymmetricMatrix *loewdinT            ,
                                                        CSymmetricMatrix *density             ,
                                                        CRealArray1D     *charges             ,
                                                        CStatus          *status              )
    cdef void  Loewdin_BondOrders                     ( CIntegerArray1D  *basisIndices        ,
                                                        CSymmetricMatrix *loewdinT            ,
                                                        CSymmetricMatrix *density             ,
                                                        CSymmetricMatrix *bondOrders          ,
                                                        CStatus          *status              )
    cdef void  Loewdin_ChargeDensityDerivatives       ( CIntegerArray1D  *basisIndices        ,
                                                        CRealArray1D     *potentials          ,
                                                        CSymmetricMatrix *loewdinT            ,
                                                        CSymmetricMatrix *fock                )
    cdef CReal Loewdin_ChargeRestraintMatrix          ( CIntegerArray1D  *basisIndices        ,
                                                        CRealArray1D     *nuclearCharges      ,
                                                        CIntegerArray1D  *crIndices           ,
                                                        CRealArray1D     *crWeights           ,
                                                        CBoolean          isSpin              ,
                                                        CSymmetricMatrix *loewdinT            ,
                                                        CSymmetricMatrix *W                   )
    cdef void  Loewdin_ChargeRestraintWeightedDensity ( CIntegerArray1D  *basisIndices        ,
                                                        CIntegerArray1D  *crIndices           ,
                                                        CRealArray1D     *crWeights           ,
                                                        CBoolean          isSpin              ,
                                                        CReal             dRdL                ,
                                                        CRealArray2D     *eigenVectors        ,
                                                        CRealArray2D     *Z                   ,
                                                        CSymmetricMatrix *A                   )
    cdef void  Loewdin_WeightedDensity                ( CIntegerArray1D  *basisIndices        ,
                                                        CRealArray1D     *potentials          ,
                                                        CRealArray1D     *eigenValues         ,
                                                        CRealArray2D     *eigenVectors        ,
                                                        CSymmetricMatrix *loewdinT            ,
                                                        CSymmetricMatrix *density             ,
                                                        CReal            *eigenValueTolerance ,
                                                        CSymmetricMatrix *wDensity            ,
                                                        CStatus          *status              )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class LoewdinMultipoleEvaluator:
    pass
