from pCore.CPrimitiveTypes           cimport CBoolean      , \
                                             CInteger      , \
                                             CFalse        , \
                                             CInteger      , \
                                             CReal         , \
                                             CTrue
from pCore.IntegerBlock              cimport CIntegerBlock , \
                                             IntegerBlock
from pCore.RealBlock                 cimport CRealBlock    , \
                                             RealBlock
from pCore.Status                    cimport CStatus       , \
                                             CStatus_OK
from pScientific.Arrays.BaseIterator cimport CIterator
from pScientific.Arrays.RealArray1D  cimport CRealArray1D  , \
                                             RealArray1D
from pScientific.Arrays.RealArray2D  cimport CRealArray2D  , \
                                             RealArray2D
from pScientific.Arrays.RealIterator cimport RealIterator

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "SymmetricMatrix.h":

    ctypedef struct CSymmetricMatrix "SymmetricMatrix":
        CInteger    extent
        CInteger    size
        CRealBlock *block
        CReal      *data

    ctypedef enum CSymmetricMatrixUpdating_Option "SymmetricMatrixUpdating_Option":
        SymmetricMatrixUpdating_BFGS   = 0
        SymmetricMatrixUpdating_Bofill = 1
        SymmetricMatrixUpdating_MS     = 2
        SymmetricMatrixUpdating_Powell = 3

    cdef CReal             SymmetricMatrix_AbsoluteMaximum          ( CSymmetricMatrix                *self            )
    cdef void              SymmetricMatrix_Add                      ( CSymmetricMatrix                *self            ,
                                                                      CReal                            alpha           ,
                                                                      CSymmetricMatrix                *other           ,
                                                                      CStatus                         *status          )
    cdef CSymmetricMatrix *SymmetricMatrix_Allocate                 ( CStatus                         *status          )
    cdef CSymmetricMatrix *SymmetricMatrix_AllocateWithExtent       ( CInteger                         extent          ,
                                                                      CStatus                         *status          )
    cdef void              SymmetricMatrix_AnticommutatorSS         ( CSymmetricMatrix                *self            ,
                                                                      CSymmetricMatrix                *a               ,
                                                                      CSymmetricMatrix                *b               ,
                                                                      CStatus                         *status          )
    cdef CSymmetricMatrix *SymmetricMatrix_CloneDeep                ( CSymmetricMatrix                *self            ,
                                                                      CStatus                         *status          )
    cdef CSymmetricMatrix *SymmetricMatrix_CloneShallow             ( CSymmetricMatrix                *self            ,
                                                                      CStatus                         *status          )
    cdef void              SymmetricMatrix_CopyFromRealArray2D      ( CSymmetricMatrix                *self            ,
                                                                      CRealArray2D                    *other           ,
                                                                      CStatus                         *status          )
    cdef void              SymmetricMatrix_CopyTo                   ( CSymmetricMatrix                *self            ,
                                                                      CSymmetricMatrix                *other           ,
                                                                      CStatus                         *status          )
    cdef void              SymmetricMatrix_CopyToRealArray2D        ( CSymmetricMatrix                *self            ,
                                                                      CRealArray2D                    *other           ,
                                                                      CStatus *status                                  )
    cdef void              SymmetricMatrix_Deallocate               ( CSymmetricMatrix               **self            )
    cdef void              SymmetricMatrix_DiagonalOfProduct        ( CSymmetricMatrix                *self            ,
                                                                      CSymmetricMatrix                *other           ,
                                                                      CRealArray1D                    *diagonal        ,
                                                                      CStatus                         *status          ) 
    cdef void              SymmetricMatrix_DiagonalOfTransform      ( CSymmetricMatrix                *self            ,
                                                                      CRealArray2D                    *matrix          ,
                                                                      CBoolean                         useTranspose    ,
                                                                      CRealArray1D                    *diagonal        ,
                                                                      CStatus                         *status          )
    cdef CSymmetricMatrix *SymmetricMatrix_FromExtentBlock          ( CInteger                         extent          ,
                                                                      CRealBlock                      *block           ,
                                                                      CBoolean                         withReference   ,
                                                                      CStatus                         *status          )
    cdef void              SymmetricMatrix_GetColumn                ( CSymmetricMatrix                *self            ,
                                                                      CInteger                         n               ,
                                                                      CRealArray1D                    *column          ,
                                                                      CStatus                         *status          )
    cdef CReal             SymmetricMatrix_GetItem                  ( CSymmetricMatrix                *self            ,
                                                                      CInteger                         i               ,
                                                                      CInteger                         j               ,
                                                                      CStatus                         *status          )
    cdef void              SymmetricMatrix_Increment                ( CSymmetricMatrix                *self            ,
                                                                      CReal                            value           )
    cdef void              SymmetricMatrix_IndexedCopyToRealArray2D ( CSymmetricMatrix                *self            ,
                                                                      CIntegerBlock                   *indices         ,
                                                                      CRealArray2D                    *target          ,
                                                                      CStatus                         *status          )
    cdef CBoolean          SymmetricMatrix_IsDiagonal               ( CSymmetricMatrix                *self            ,
                                                                      CReal                            tolerance       )
    cdef void              SymmetricMatrix_MakeFromEigensystem      ( CSymmetricMatrix                *self            ,
                                                                      CBoolean                         zeroMatrix      ,
                                                                      CInteger                         numberOfVectors ,
                                                                      CRealArray1D                    *eigenValues     ,
                                                                      CRealArray2D                    *eigenVectors    ,
                                                                      CStatus                         *status          )
    cdef CIterator        *SymmetricMatrix_MakeIterator             ( CSymmetricMatrix                *self            ,
                                                                      CStatus                         *status          )
    cdef void              SymmetricMatrix_PostMatrixMultiply       ( CSymmetricMatrix                *self            ,
                                                                      CRealArray2D                    *matrix          ,
                                                                      CBoolean                         useTranspose    ,
                                                                      CRealArray2D                    *result          ,
                                                                      CStatus                         *status          )
    cdef void              SymmetricMatrix_PreMatrixMultiply        ( CSymmetricMatrix                *self            ,
                                                                      CRealArray2D                    *matrix          ,
                                                                      CBoolean                         useTranspose    ,
                                                                      CRealArray2D                    *result          ,
                                                                      CStatus                         *status          )
    cdef void              SymmetricMatrix_Print                    ( CSymmetricMatrix                *self            )
    cdef void              SymmetricMatrix_ProjectionMatrix         ( CSymmetricMatrix                *self            ,
                                                                      CRealArray2D                    *vectors         ,
                                                                      CStatus                         *status          )
    cdef void              SymmetricMatrix_ProjectOut               ( CSymmetricMatrix                *self            ,
                                                                      CRealArray2D                    *vectors         ,
                                                                      CSymmetricMatrix                *result          ,
                                                                      CStatus                         *status          )
    cdef void              SymmetricMatrix_Rank1Update              ( CSymmetricMatrix                *self            ,
                                                                      CReal                            alpha           ,
                                                                      CRealArray1D                    *vector          ,
                                                                      CStatus                         *status          )
    cdef void              SymmetricMatrix_Raise                    ( CSymmetricMatrix                *self            ,
                                                                      CRealArray2D                    *vectors         ,
                                                                      CReal                            value           ,
                                                                      CStatus                         *status          )
    cdef CReal             SymmetricMatrix_RootMeanSquare           ( CSymmetricMatrix                *self            )
    cdef void              SymmetricMatrix_Scale                    ( CSymmetricMatrix                *self            ,
                                                                      CReal                            value           )
    cdef void              SymmetricMatrix_ScaleOffDiagonal         ( CSymmetricMatrix                *self            ,
                                                                      CReal                            value           )
    cdef void              SymmetricMatrix_Set                      ( CSymmetricMatrix                *self            ,
                                                                      CReal                            value           )
    cdef void              SymmetricMatrix_SetColumn                ( CSymmetricMatrix                *self            ,
                                                                      CInteger                         column          ,
                                                                      CRealArray1D                    *vector          ,
                                                                      CStatus                         *status          )
    cdef void              SymmetricMatrix_SetItem                  ( CSymmetricMatrix                *self            ,
                                                                      CInteger                         i               ,
                                                                      CInteger                         j               ,
                                                                      CReal                            value           ,
                                                                      CStatus                         *status          )
    cdef CReal             SymmetricMatrix_Sparsity                 ( CSymmetricMatrix                *self            ,
                                                                      CReal                            tolerance       )
    cdef void              SymmetricMatrix_SumDifference            ( CSymmetricMatrix                *self            ,
                                                                      CSymmetricMatrix                *other           ,
                                                                      CStatus                         *status          )
    cdef void              SymmetricMatrix_SymmetricMatrixMultiply  ( CSymmetricMatrix                *self            ,
                                                                      CSymmetricMatrix                *other           ,
                                                                      CRealArray2D                    *result          ,
                                                                      CStatus                         *status          )
    cdef void              SymmetricMatrix_SymmetricTransform       ( CSymmetricMatrix                *self            ,
                                                                      CSymmetricMatrix                *other           ,
                                                                      CSymmetricMatrix                *result          ,
                                                                      CStatus                         *status          )
    cdef CReal             SymmetricMatrix_Trace                    ( CSymmetricMatrix                *self            )
    cdef CReal             SymmetricMatrix_TraceOfProduct           ( CSymmetricMatrix                *self            ,
                                                                      CSymmetricMatrix                *other           ,
                                                                      CStatus                         *status          )
    cdef void              SymmetricMatrix_Transform                ( CSymmetricMatrix                *self            ,
                                                                      CRealArray2D                    *matrix          ,
                                                                      CBoolean                         useTranspose    ,
                                                                      CSymmetricMatrix                *result          ,
                                                                      CStatus                         *status          )
    cdef void              SymmetricMatrix_Update                   ( CSymmetricMatrix                *self            ,
                                                                      CRealArray1D                    *dx              ,
                                                                      CRealArray1D                    *dg              ,
                                                                      CSymmetricMatrixUpdating_Option  option          ,
                                                                      CReal                           *tolerance       ,
                                                                      CStatus                         *status          )
    cdef void              SymmetricMatrix_VectorMultiply           ( CSymmetricMatrix                *self            ,
                                                                      CRealArray1D                    *other           ,
                                                                      CRealArray1D                    *result          ,
                                                                      CStatus                         *status          )
    cdef CInteger          SymmetricMatrix_ViewSize                 ( CInteger                         extent          )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class SymmetricMatrix:

    cdef CSymmetricMatrix    *cObject
    cdef public RealBlock     block
    cdef        RealIterator _iterator
