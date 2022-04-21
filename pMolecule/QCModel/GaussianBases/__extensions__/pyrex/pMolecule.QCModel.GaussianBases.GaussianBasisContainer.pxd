from pCore.CPrimitiveTypes                         cimport CBoolean                         , \
                                                           CFalse                           , \
                                                           CInteger                         , \
                                                           CReal                            , \
                                                           CTrue
from pCore.Status                                  cimport CStatus                          , \
                                                           CStatus_OK
from pMolecule.QCModel.GaussianBases.GaussianBasis cimport CGaussianBasis                   , \
                                                           GaussianBasis                    , \
                                                           GaussianBasis_MakeLRotations     , \
                                                           GaussianBasis_MakeRotationMatrix
from pScientific.Arrays.IntegerArray1D             cimport CIntegerArray1D                  , \
                                                           IntegerArray1D
from pScientific.Arrays.RealArray1D                cimport RealArray1D
from pScientific.Arrays.RealArray2D                cimport CRealArray2D                     , \
                                                           RealArray2D

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "GaussianBasisContainer.h":

    ctypedef struct CGaussianBasisContainer "GaussianBasisContainer":
        CBoolean         isOwner 
        CInteger         capacity
        CIntegerArray1D *centerFunctionPointers
        CIntegerArray1D *functionCenters
        CGaussianBasis **entries 

    cdef CGaussianBasisContainer *GaussianBasisContainer_Allocate              ( CInteger                  capacity               ,
                                                                                 CStatus                  *status                 )
    cdef CGaussianBasisContainer *GaussianBasisContainer_Clone                 ( CGaussianBasisContainer  *self                   ,
                                                                                 CStatus                  *status                 )
    cdef void                     GaussianBasisContainer_Deallocate            ( CGaussianBasisContainer **self                   )
    # . For debugging only.
#   cdef void                     GaussianBasisContainer_MakeC2S               ( CGaussianBasisContainer  *self                   ,
#                                                                                CRealArray2D             *T                      ,
#                                                                                CStatus                  *status                 )
    #
    cdef void                     GaussianBasisContainer_MakeIndexArrays       ( CGaussianBasisContainer  *self                   ,
                                                                                 CIntegerArray1D          *centerFunctionPointers ,
                                                                                 CIntegerArray1D          *functionCenters        ,
                                                                                 CStatus                  *status                 )
    cdef CInteger                 GaussianBasisContainer_NumberOfFunctions     ( CGaussianBasisContainer  *self                   )
    cdef CInteger                 GaussianBasisContainer_NumberOfWorkFunctions ( CGaussianBasisContainer  *self                   )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class GaussianBasisContainer:

    cdef CGaussianBasisContainer *cObject
    cdef public object            isOwner
    cdef public object            label
    cdef public object            uniqueEntries
    cdef public object            _centerFunctionPointers
    cdef public object            _functionCenters
    cdef public object            _functionLabels
    cdef public object            _numberOfFunctions
    cdef public object            _numberOfWorkFunctions
