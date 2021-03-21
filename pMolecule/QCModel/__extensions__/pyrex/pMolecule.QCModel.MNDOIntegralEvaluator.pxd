from pCore.CPrimitiveTypes                     cimport CBoolean                    , \
                                                       CFalse                      , \
                                                       CInteger                    , \
                                                       CReal                       , \
                                                       CTrue
from pCore.Status                              cimport CStatus                     , \
                                                       CStatus_OK
from pMolecule.QCModel.BlockStorage            cimport BlockStorage                , \
                                                       CBlockStorage
from pMolecule.QCModel.GaussianBasisContainer  cimport CGaussianBasisContainer     , \
                                                       GaussianBasisContainer
from pMolecule.QCModel.MNDOParametersContainer cimport CMNDOParametersContainer    , \
                                                       MNDOParametersContainer
from pScientific.Arrays.DoubleSymmetricMatrix  cimport CDoubleSymmetricMatrix      , \
                                                       DoubleSymmetricMatrix
from pScientific.Arrays.IntegerArray1D         cimport IntegerArray1D              , \
                                                       CIntegerArray1D      
from pScientific.Arrays.RealArray1D            cimport CRealArray1D                , \
                                                       RealArray1D
from pScientific.Arrays.RealArray2D            cimport CRealArray2D                , \
                                                       RealArray2D
from pScientific.Arrays.SymmetricMatrix        cimport CSymmetricMatrix            , \
                                                       SymmetricMatrix
from pScientific.Geometry3.Coordinates3        cimport Coordinates3
from pScientific.Geometry3.Vector3             cimport Vector3

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "MNDOCoreCore.h":

    cdef CReal MNDO_CoreCoreEnergy               ( CMNDOParametersContainer *parameters               ,
                                                   CRealArray2D             *coordinates3             ,
                                                   CRealArray2D             *gradients3               )

cdef extern from "MNDODipole.h":

    cdef void MNDO_DipoleIntegrals               ( CMNDOParametersContainer *parameters               ,
                                                   CIntegerArray1D          *basisIndices             ,
                                                   CRealArray2D             *coordinates3             ,
                                                   CRealArray1D             *center                   ,
                                                   CSymmetricMatrix         *dX                       ,
                                                   CSymmetricMatrix         *dY                       ,
                                                   CSymmetricMatrix         *dZ                       ,
                                                   CStatus                  *status                   )

cdef extern from "MNDOElectronNuclearTEIs.h":

    cdef void MNDO_ElectronNuclearTEIGradients   ( CMNDOParametersContainer *parameters               ,
                                                   CIntegerArray1D          *basisIndices             ,
                                                   CRealArray2D             *coordinates3             ,
                                                   CSymmetricMatrix         *dTotal                   ,
                                                   CSymmetricMatrix         *dSpin                    ,
                                                   CRealArray2D             *gradients3               )
    cdef void MNDO_ElectronNuclearTEIGradientsCI ( CInteger                  nActive                  ,
                                                   CInteger                  nCore                    ,
                                                   CInteger                  nOrbitals                ,
                                                   CMNDOParametersContainer *parameters               ,
                                                   CIntegerArray1D          *basisIndices             ,
                                                   CRealArray2D             *coordinates3             ,
                                                   CDoubleSymmetricMatrix   *twoPDM                   ,
                                                   CRealArray2D             *orbitals                 ,
                                                   CSymmetricMatrix         *dCore                    ,
                                                   CSymmetricMatrix         *dHF                      ,
                                                   CSymmetricMatrix         *dTotal                   ,
                                                   CSymmetricMatrix         *onePDM                   ,
                                                   CSymmetricMatrix         *zMatrix                  ,
                                                   CRealArray2D             *gradients3               )
    cdef void MNDO_ElectronNuclearTEIIntegrals   ( CMNDOParametersContainer *parameters               ,
                                                   CIntegerArray1D          *basisIndices             ,
                                                   CRealArray2D             *coordinates3             ,
                                                   CSymmetricMatrix         *oneElectronMatrix        ,
                                                   CBlockStorage           **twoElectronIntegrals     )

cdef extern from "MNDOResonance.h":

    cdef void MNDO_ResonanceGradients            ( CMNDOParametersContainer *parameters               ,
                                                   CGaussianBasisContainer  *bases                    ,
                                                   CIntegerArray1D          *basisIndices             ,
                                                   CRealArray2D             *coordinates3             ,
                                                   CSymmetricMatrix         *dTotal                   ,
                                                   CRealArray2D             *gradients3               )
    cdef void MNDO_ResonanceIntegrals            ( CMNDOParametersContainer *parameters               ,
                                                   CGaussianBasisContainer  *bases                    ,
                                                   CIntegerArray1D          *basisIndices             ,
                                                   CRealArray2D             *coordinates3             ,
                                                   CSymmetricMatrix         *oneElectronMatrix        )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class MNDOIntegralEvaluator:
    pass
