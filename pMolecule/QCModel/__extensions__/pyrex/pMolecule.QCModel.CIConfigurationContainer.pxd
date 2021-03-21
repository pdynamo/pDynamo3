from pCore.CPrimitiveTypes                    cimport CBoolean                           , \
                                                      CFalse                             , \
                                                      CInteger                           , \
                                                      CReal                              , \
                                                      CTrue
from pCore.Status                             cimport CStatus                            , \
                                                      CStatus_OK
from pScientific.Arrays.BooleanArray1D        cimport CBooleanArray1D
from pScientific.Arrays.DoubleSymmetricMatrix cimport CDoubleSymmetricMatrix             , \
                                                      DoubleSymmetricMatrix
from pScientific.Arrays.IntegerArray1D        cimport IntegerArray1D_GetItem             , \
                                                      CIntegerArray1D
from pScientific.Arrays.IntegerArray2D        cimport IntegerArray2D_AllocateWithExtents , \
                                                      IntegerArray2D_Deallocate          , \
                                                      IntegerArray2D_Set                 , \
                                                      IntegerArray2D_SetItem             , \
                                                      CIntegerArray2D
from pScientific.Arrays.RealArray1D           cimport CRealArray1D                       , \
                                                      RealArray1D
from pScientific.Arrays.RealArray2D           cimport CRealArray2D                       , \
                                                      RealArray2D
from pScientific.Arrays.SparseSymmetricMatrix cimport CSparseSymmetricMatrix             , \
                                                      SparseSymmetricMatrix
from pScientific.Arrays.SymmetricMatrix       cimport CSymmetricMatrix                   , \
                                                      SymmetricMatrix

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "CIConfigurationContainer.h":

    ctypedef struct CCIConfiguration "CIConfiguration":
        CInteger          nAlphas
        CInteger          nSPQR
        CReal             spin
        CBooleanArray1D  *parity
        CIntegerArray1D  *alphas
        CIntegerArray1D  *betas
        CIntegerArray1D  *spqr

    ctypedef struct CCIConfigurationContainer "CIConfigurationContainer":
        CInteger          nActive
        CInteger          nConfigurations
        CInteger          nElectrons
        CCIConfiguration *configurations


    cdef CCIConfigurationContainer *CIConfigurationContainer_Allocate                ( CInteger                    nActive               ,
                                                                                       CInteger                    nConfigurations       ,
                                                                                       CStatus                    *status                )
    cdef void                       CIConfigurationContainer_Characters              ( CCIConfigurationContainer  *self                  ,
                                                                                       CBoolean                    includeCoreOrbitals   ,
                                                                                       CInteger                    coreOrbitals          ,
                                                                                       CRealArray2D               *orbitalTransformation ,
                                                                                       CRealArray2D               *stateTransformation   ,
                                                                                       CStatus                    *status                )
    cdef CCIConfigurationContainer *CIConfigurationContainer_Clone                   ( CCIConfigurationContainer  *self                  ,
                                                                                       CStatus                    *status                )
    cdef void                       CIConfigurationContainer_Deallocate              ( CCIConfigurationContainer **self                  )
    cdef void                       CIConfigurationContainer_GetCIMatrixSparsity     ( CCIConfigurationContainer  *self                  ,
                                                                                       CInteger                   *nonZero               ,
                                                                                       CReal                      *sparsity              )
    cdef void                       CIConfigurationContainer_MakeCIMatrix            ( CCIConfigurationContainer  *self                  ,
                                                                                       CSymmetricMatrix           *fCoreMO               ,
                                                                                       CDoubleSymmetricMatrix     *moTEIs                ,
                                                                                       CSymmetricMatrix           *ciMatrixFull          ,
                                                                                       CSparseSymmetricMatrix     *ciMatrixSparse        ,
                                                                                       CStatus                    *status                )
    cdef void                       CIConfigurationContainer_MakeDensities           ( CCIConfigurationContainer  *self                  ,
                                                                                       CRealArray1D               *ciVector              ,
                                                                                       CSymmetricMatrix           *onePDMMOa             ,
                                                                                       CSymmetricMatrix           *onePDMMOb             ,
                                                                                       CDoubleSymmetricMatrix     *twoPDM                ,
                                                                                       CStatus                    *status                )
    cdef CCIConfigurationContainer *CIConfigurationContainer_MakeFull                ( CInteger                    nActive               ,
                                                                                       CInteger                    nUp                   ,
                                                                                       CInteger                    nDown                 ,
                                                                                       CStatus                    *status                )
    cdef CCIConfigurationContainer *CIConfigurationContainer_MakeSinglesDoubles      ( CBoolean                    doSingles             ,
                                                                                       CBoolean                    doDoubles             ,
                                                                                       CInteger                    nActive               ,
                                                                                       CInteger                    nClosed               ,
                                                                                       CInteger                    nOpen                 ,
                                                                                       CStatus                    *status                )
    cdef CCIConfigurationContainer *CIConfigurationContainer_MakeUserSpecified       ( CIntegerArray2D            *microStates           ,
                                                                                       CInteger                    activeOrbitals        ,
                                                                                       CInteger                    activeElectrons       ,
                                                                                       CStatus                    *status                )
    cdef CInteger                   CIConfigurationContainer_NumberOfActiveElectrons ( CCIConfigurationContainer  *self                  )
    cdef CInteger                   CIConfigurationContainer_NumberOfActiveOrbitals  ( CCIConfigurationContainer  *self                  )
    cdef CInteger                   CIConfigurationContainer_NumberOfConfigurations  ( CCIConfigurationContainer  *self                  )
    cdef void                       CIConfigurationContainer_StateSpins              ( CCIConfigurationContainer  *self                  ,
                                                                                       CRealArray2D               *vectors               ,
                                                                                       CRealArray1D               *spins                 ,
                                                                                       CStatus                    *status                )
    cdef void                       CIConfigurationContainer_TransitionDipoles       ( CCIConfigurationContainer  *self                  ,
                                                                                       CSymmetricMatrix           *tdMOs                 ,
                                                                                       CSymmetricMatrix           *tdMatrix              )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class CIConfigurationContainer:

    cdef CCIConfigurationContainer *cObject
    cdef public object              isOwner
    cdef        object              _microStates
