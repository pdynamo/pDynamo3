from pCore.CPrimitiveTypes              cimport CBoolean         , \
                                                CInteger         , \
                                                CFalse           , \
                                                CInteger         , \
                                                CReal            , \
                                                CTrue
from pCore.RealBlock                    cimport CRealBlock       , \
                                                RealBlock
from pCore.Status                       cimport CStatus          , \
                                                CStatus_OK
from pScientific.Arrays.BaseIterator    cimport CIterator
from pScientific.Arrays.RealArray1D     cimport CRealArray1D     , \
                                                RealArray1D
from pScientific.Arrays.RealArray2D     cimport CRealArray2D     , \
                                                RealArray2D
from pScientific.Arrays.RealIterator    cimport RealIterator
from pScientific.Arrays.SymmetricMatrix cimport CSymmetricMatrix , \
                                                SymmetricMatrix

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "AntisymmetricMatrix.h":

    ctypedef struct CAntisymmetricMatrix "AntisymmetricMatrix":
        CInteger    extent
        CInteger    size
        CRealBlock *block
        CReal      *data

    cdef CReal                 AntisymmetricMatrix_AbsoluteMaximum        ( CAntisymmetricMatrix  *self          )
    cdef void                  AntisymmetricMatrix_Add                    ( CAntisymmetricMatrix  *self          ,
                                                                            CReal                  value         ,
                                                                            CAntisymmetricMatrix  *other         ,
                                                                            CStatus               *status        )
    cdef CAntisymmetricMatrix *AntisymmetricMatrix_Allocate               ( CStatus               *status        )
    cdef CAntisymmetricMatrix *AntisymmetricMatrix_AllocateWithExtent     ( CInteger               extent        ,
                                                                            CStatus               *status        )
    cdef void                  AntisymmetricMatrix_AnticommutatorAS       ( CAntisymmetricMatrix  *self          ,
                                                                            CSymmetricMatrix      *a             ,
                                                                            CAntisymmetricMatrix  *result        ,
                                                                            CStatus               *status        )
    cdef CAntisymmetricMatrix *AntisymmetricMatrix_CloneDeep              ( CAntisymmetricMatrix  *self          ,
                                                                            CStatus               *status        )
    cdef CAntisymmetricMatrix *AntisymmetricMatrix_CloneShallow           ( CAntisymmetricMatrix  *self          ,
                                                                            CStatus               *status        )
    cdef void                  AntisymmetricMatrix_CommutatorAS           ( CAntisymmetricMatrix  *self          ,
                                                                            CSymmetricMatrix      *a             ,
                                                                            CSymmetricMatrix      *result        ,
                                                                            CStatus               *status        )
    cdef void                  AntisymmetricMatrix_CommutatorSS_Fast      ( CAntisymmetricMatrix  *self                                                                      ,
                                                                            CSymmetricMatrix      *a             ,
                                                                            CSymmetricMatrix      *b             ,
                                                                            CRealArray2D          *mA            ,
                                                                            CRealArray2D          *mB            ,
                                                                            CRealArray2D          *mC            ,
                                                                            CStatus               *status        )
    cdef void                  AntisymmetricMatrix_CommutatorSS_Reference ( CAntisymmetricMatrix  *self          ,
                                                                            CSymmetricMatrix      *a             ,
                                                                            CSymmetricMatrix      *b             ,
                                                                            CStatus               *status        )
    cdef void                  AntisymmetricMatrix_CommutatorSSS          ( CAntisymmetricMatrix  *self          ,
                                                                            CSymmetricMatrix      *a             ,
                                                                            CSymmetricMatrix      *b             ,
                                                                            CSymmetricMatrix      *c             ,
                                                                            CStatus               *status        )
    cdef void                  AntisymmetricMatrix_CommutatorTSSST        ( CAntisymmetricMatrix  *self          ,
                                                                            CSymmetricMatrix      *a             ,
                                                                            CSymmetricMatrix      *b             ,
                                                                            CSymmetricMatrix      *c             ,
                                                                            CRealArray2D          *m             ,
                                                                            CBoolean               mTranspose    ,
                                                                            CRealArray2D          *u             ,
                                                                            CRealArray2D          *v             ,
                                                                            CRealArray2D          *w             ,
                                                                            CStatus               *status        )
    cdef void                  AntisymmetricMatrix_CommutatorXSSY         ( CAntisymmetricMatrix  *self          ,
                                                                            CSymmetricMatrix      *a             ,
                                                                            CSymmetricMatrix      *b             ,
                                                                            CRealArray2D          *x             ,
                                                                            CRealArray2D          *y             ,
                                                                            CBoolean               xTranspose    ,
                                                                            CBoolean               yTranspose    ,
                                                                            CRealArray2D          *u             ,
                                                                            CRealArray2D          *v             ,
                                                                            CRealArray2D          *w             ,
                                                                            CStatus               *status        )
    cdef void                  AntisymmetricMatrix_CopyFromRealArray2D    ( CAntisymmetricMatrix  *self          ,
                                                                            CRealArray2D          *other         ,
                                                                            CBoolean               scale         ,
                                                                            CStatus               *status        )
    cdef void                  AntisymmetricMatrix_CopyTo                 ( CAntisymmetricMatrix  *self          ,
                                                                            CAntisymmetricMatrix  *other         ,
                                                                            CStatus               *status        )
    cdef void                  AntisymmetricMatrix_CopyToRealArray2D      ( CAntisymmetricMatrix  *self          ,
                                                                            CRealArray2D          *other         ,
                                                                            CStatus               *status        )
    cdef void                  AntisymmetricMatrix_Deallocate             ( CAntisymmetricMatrix **self          )
    cdef CAntisymmetricMatrix *AntisymmetricMatrix_FromExtentBlock        ( CInteger               extent        ,
                                                                            CRealBlock            *block         ,
                                                                            CBoolean               withReference ,
                                                                            CStatus               *status        )
    cdef void                  AntisymmetricMatrix_GetColumn              ( CAntisymmetricMatrix  *self          ,
                                                                            CInteger               n             ,
                                                                            CRealArray1D          *column        ,
                                                                            CStatus               *status        )
    cdef CReal                 AntisymmetricMatrix_GetItem                ( CAntisymmetricMatrix  *self          ,
                                                                            CInteger               i             ,
                                                                            CInteger               j             ,
                                                                            CStatus               *status        )
    cdef void                  AntisymmetricMatrix_GetItemIndexAndSign    ( CAntisymmetricMatrix  *self          ,
                                                                            CInteger               i             ,
                                                                            CInteger               j             ,
                                                                            CInteger              *index         ,
                                                                            CReal                 *sign          ,
                                                                            CStatus               *status        )
    cdef CIterator            *AntisymmetricMatrix_MakeIterator           ( CAntisymmetricMatrix  *self          ,
                                                                            CStatus               *status        )
    cdef void                  AntisymmetricMatrix_Print                  ( CAntisymmetricMatrix  *self          )
    cdef void                  AntisymmetricMatrix_Scale                  ( CAntisymmetricMatrix  *self          ,
                                                                            CReal                  value         )
    cdef void                  AntisymmetricMatrix_Set                    ( CAntisymmetricMatrix  *self          ,
                                                                            CReal                  value         )
    cdef void                  AntisymmetricMatrix_SetItem                ( CAntisymmetricMatrix  *self          ,
                                                                            CInteger               i             ,
                                                                            CInteger               j             ,
                                                                            CReal                  value         ,
                                                                            CStatus               *status        )
    cdef void                  AntisymmetricMatrix_SymmetricTransform     ( CAntisymmetricMatrix  *self          ,
                                                                            CSymmetricMatrix      *matrix        ,
                                                                            CAntisymmetricMatrix  *result        ,
                                                                            CStatus               *status        )
    cdef CReal                 AntisymmetricMatrix_TraceOfProduct         ( CAntisymmetricMatrix  *self          ,
                                                                            CAntisymmetricMatrix  *other         ,
                                                                            CStatus               *status        )
    cdef void                  AntisymmetricMatrix_Transform              ( CAntisymmetricMatrix  *self          ,
                                                                            CRealArray2D          *matrix        ,
                                                                            CBoolean               useTranspose  ,
                                                                            CAntisymmetricMatrix  *result        ,
                                                                            CStatus               *status        )
    cdef void                  AntisymmetricMatrix_Transpose              ( CAntisymmetricMatrix  *self          )
    cdef void                  AntisymmetricMatrix_ViewOfRaw              ( CAntisymmetricMatrix  *self          ,
                                                                            CInteger               extent        ,
                                                                            CInteger               stride        ,
                                                                            CReal                 *data          ,
                                                                            CInteger               rawstride     ,
                                                                            CStatus               *status        )
    cdef CInteger              AntisymmetricMatrix_ViewSize               ( CInteger               extent        )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class AntisymmetricMatrix:

    cdef CAntisymmetricMatrix *cObject
    cdef public RealBlock      block
    cdef        RealIterator  _iterator
