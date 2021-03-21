from pCore.CPrimitiveTypes                    cimport CBoolean                 , \
                                                      CFalse                   , \
                                                      CInteger                 , \
                                                      CReal                    , \
                                                      CTrue
from pCore.Status                             cimport CStatus                  , \
                                                      CStatus_OK
from pMolecule.QCModel.BlockStorage           cimport CBlockStorage            , \
                                                      BlockStorage              
from pScientific.Arrays.DoubleSymmetricMatrix cimport CDoubleSymmetricMatrix   , \
                                                      DoubleSymmetricMatrix
from pScientific.Arrays.IntegerArray2D        cimport CIntegerArray2D          , \
                                                      IntegerArray2D
from pScientific.Arrays.RealArray1D           cimport CRealArray1D             , \
                                                      RealArray1D
from pScientific.Arrays.RealArray2D           cimport CRealArray2D             , \
                                                      RealArray2D
from pScientific.Arrays.RealArrayND           cimport CRealArrayND             , \
                                                      RealArrayND
from pScientific.Arrays.SymmetricMatrix       cimport CSymmetricMatrix         , \
                                                      SymmetricMatrix

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "CICPHF.h":

    cdef void CCICPHF_ApplyCPHFMatrix      "CICPHF_ApplyCPHFMatrix"      ( CInteger                n1                        ,
                                                                           CIntegerArray2D        *in1                       ,
                                                                           CInteger                n2                        ,
                                                                           CIntegerArray2D        *in2                       ,
                                                                           CRealArray1D           *aDiagonal                 ,
                                                                           CRealArray1D           *b                         ,
                                                                           CRealArray2D           *orbitals                  ,
                                                                           CBlockStorage          *twoElectronIntegrals      ,
                                                                           CSymmetricMatrix       *work1                     ,
                                                                           CSymmetricMatrix       *work2                     ,
                                                                           CRealArray1D           *x                         )
    cdef void CCICPHF_CalculateCPHFVectors "CICPHF_CalculateCPHFVectors" ( CInteger                nActive                   ,
                                                                           CInteger                nCore                     ,
                                                                           CInteger                nOrbitals                 ,
                                                                           CBlockStorage          *twoElectronIntegrals      ,
                                                                           CDoubleSymmetricMatrix *twoPDM                    ,
                                                                           CRealArray1D           *energies                  ,
                                                                           CRealArray1D           *occupancies               ,
                                                                           CRealArray2D           *orbitals                  ,
                                                                           CRealArrayND           *moTEI234                  ,
                                                                           CSymmetricMatrix       *fCore                     ,
                                                                           CSymmetricMatrix       *onePDM                    ,
                                                                           CSymmetricMatrix       *onePDMMO                  ,
                                                                           CSymmetricMatrix       *work1                     ,
                                                                           CSymmetricMatrix       *work2                     ,
                                                                           CInteger               *numberDegenerateRedundant ,
                                                                           CInteger               *numberNonRedundant        ,
                                                                           CInteger               *numberRedundant           ,
                                                                           CIntegerArray2D        *indicesNR                 ,
                                                                           CIntegerArray2D        *indicesR                  ,
                                                                           CRealArray1D           *aDiagonal                 ,
                                                                           CRealArray1D           *qNR                       ,
                                                                           CRealArray1D           *qR                        ,
                                                                           CRealArray1D           *preconditioner            ,
                                                                           CStatus                *status                    )
    cdef void CCICPHF_Transform            "CICPHF_Transform"            ( CInteger                n1                        ,
                                                                           CIntegerArray2D        *in1                       ,
                                                                           CRealArray1D           *x1                        ,
                                                                           CInteger                n2                        ,
                                                                           CIntegerArray2D        *in2                       ,
                                                                           CRealArray1D           *x2                        ,
                                                                           CRealArray2D           *orbitals                  ,
                                                                           CBoolean                doScale                   ,
                                                                           CSymmetricMatrix       *work                      ,
                                                                           CSymmetricMatrix       *z                         )
