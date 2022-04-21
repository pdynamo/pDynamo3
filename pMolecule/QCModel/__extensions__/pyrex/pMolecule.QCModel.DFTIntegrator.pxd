from pCore.CPrimitiveTypes                                  cimport CBoolean                , \
                                                                    CFalse                  , \
                                                                    CTrue                   , \
                                                                    CInteger                , \
                                                                    CReal
from pCore.Status                                           cimport CStatus                 , \
                                                                    CStatus_OK
from pMolecule.QCModel.DFTFunctionalModel                   cimport CDFTFunctionalModel     , \
                                                                    DFTFunctionalModel
from pMolecule.QCModel.DFTGrid                              cimport CDFTGrid                , \
                                                                    DFTGrid
from pMolecule.QCModel.GaussianBases.GaussianBasisContainer cimport CGaussianBasisContainer , \
                                                                    GaussianBasisContainer
from pScientific.Arrays.RealArray2D                         cimport CRealArray2D
from pScientific.Arrays.SymmetricMatrix                     cimport CSymmetricMatrix        , \
                                                                    SymmetricMatrix      
from pScientific.Geometry3.Coordinates3                     cimport Coordinates3     

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "DFTIntegrator.h":

    cdef void CDFTIntegrator_Integrate "DFTIntegrator_Integrate" ( CDFTFunctionalModel     *functionalModel    ,
                                                                   CDFTGrid                *grid               ,
                                                                   CGaussianBasisContainer *gaussianBases      ,
                                                                   CRealArray2D            *qcCoordinates      ,
                                                                   CSymmetricMatrix        *densityP           ,
                                                                   CSymmetricMatrix        *densityQ           ,
                                                                   CBoolean                 inCore             ,
                                                                   CBoolean                 isSpinUnrestricted ,
                                                                   CReal                   *eQuad              ,
                                                                   CReal                   *rhoQuad            ,
                                                                   CSymmetricMatrix        *fockA              ,
                                                                   CSymmetricMatrix        *fockB              ,
                                                                   CRealArray2D            *gradients3         ,
                                                                   CStatus                 *status             )

