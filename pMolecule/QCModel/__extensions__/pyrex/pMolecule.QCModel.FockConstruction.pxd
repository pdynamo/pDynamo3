from pCore.CPrimitiveTypes              cimport CBoolean         , \
                                                CFalse           , \
                                                CTrue            , \
                                                CInteger         , \
                                                CReal
from pCore.Status                       cimport CStatus          , \
                                                CStatus_OK
from pMolecule.QCModel.BlockStorage     cimport BlockStorage     , \
                                                CBlockStorage
from pScientific.Arrays.RealArray1D     cimport CRealArray1D     , \
                                                 RealArray1D        
from pScientific.Arrays.SymmetricMatrix cimport CSymmetricMatrix , \
                                                SymmetricMatrix        

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "FockConstruction.h":

    cdef CReal Fock_MakeFromFitIntegrals ( CBlockStorage    *fitIntegrals         ,
                                           CSymmetricMatrix *fitMatrix            ,
                                           CReal             totalCharge          ,
                                           CRealArray1D     *fitPotential         ,
                                           CSymmetricMatrix *dTotal               ,
                                           CSymmetricMatrix *fTotal               ,
                                           CStatus          *status               )
    cdef CReal Fock_MakeFromTEIs         ( CBlockStorage    *twoElectronIntegrals ,
                                           CSymmetricMatrix *dTotal               ,
                                           CSymmetricMatrix *dSpin                ,
                                           CReal             exchangeScaling      ,
                                           CSymmetricMatrix *fTotal               ,
                                           CSymmetricMatrix *fSpin                )
    cdef CReal Fock_MakeFromTEIsCoulomb  ( CBlockStorage    *twoElectronIntegrals ,
                                           CSymmetricMatrix *dTotal               ,
                                           CSymmetricMatrix *fTotal               )
    cdef CReal Fock_MakeFromTEIsExchange ( CBlockStorage    *twoElectronIntegrals ,
                                           CSymmetricMatrix *dTotal               ,
                                           CSymmetricMatrix *dSpin                ,
                                           CReal             exchangeScaling      ,
                                           CSymmetricMatrix *fTotal               ,
                                           CSymmetricMatrix *fSpin                )
