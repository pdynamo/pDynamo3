from pCore.CPrimitiveTypes             cimport CBoolean                         , \
                                               CFalse                           , \
                                               CInteger                         , \
                                               CReal                            , \
                                               CTrue
from pCore.Status                      cimport CStatus                          , \
                                               CStatus_OK
from pMolecule.QCModel.GaussianBasis   cimport CGaussianBasis                   , \
                                               GaussianBasis                    , \
                                               GaussianBasis_MakeLRotations     , \
                                               GaussianBasis_MakeRotationMatrix
from pScientific.Arrays.IntegerArray1D cimport CIntegerArray1D                  , \
                                               IntegerArray1D
from pScientific.Arrays.RealArray1D    cimport RealArray1D
from pScientific.Arrays.RealArray2D    cimport CRealArray2D                     , \
                                               RealArray2D

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "GaussianBasisContainer.h":

    ctypedef struct CGaussianBasisContainer "GaussianBasisContainer":
        CBoolean         isOwner 
        CInteger         capacity
        CGaussianBasis **entries 

    cdef CGaussianBasisContainer *GaussianBasisContainer_Allocate                    ( CInteger                  capacity ,
                                                                                       CStatus                  *status   )
    cdef CGaussianBasisContainer *GaussianBasisContainer_Clone                       ( CGaussianBasisContainer  *self     ,
                                                                                       CStatus                  *status   )
    cdef void                     GaussianBasisContainer_Deallocate                  ( CGaussianBasisContainer **self     )
    cdef void                     GaussianBasisContainer_MakeBasisAtomIndices        ( CGaussianBasisContainer  *self     ,
                                                                                       CBoolean                  doWork   ,
                                                                                       CIntegerArray1D          *indices  ,
                                                                                       CStatus                  *status   )
    cdef void                     GaussianBasisContainer_MakeBasisIndices            ( CGaussianBasisContainer  *self     ,
                                                                                       CBoolean                  doWork   ,
                                                                                       CIntegerArray1D          *indices  ,
                                                                                       CStatus                  *status   )
    cdef void                     GaussianBasisContainer_MakeFunctionTransformations ( CGaussianBasisContainer  *self     ,
                                                                                       CRealArray2D             *w2a      ,
                                                                                       CRealArray2D             *a2w      ,
                                                                                       CStatus                  *status   )
    cdef CInteger                 GaussianBasisContainer_NumberOfBasisFunctions      ( CGaussianBasisContainer  *self     ,
                                                                                       CBoolean                  doWork   )

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
    cdef public object            _numberOfFunctions
    cdef        object            _representation
    cdef RealArray2D              _a2w
    cdef RealArray2D              _w2a
