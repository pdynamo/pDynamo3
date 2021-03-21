from pCore.CPrimitiveTypes          cimport CBoolean     , \
                                            CFalse       , \
                                            CInteger     , \
                                            CReal        , \
                                            CTrue
from pCore.RealBlock                cimport CRealBlock   , \
                                            RealBlock
from pCore.Status                   cimport CStatus      , \
                                            CStatus_OK
from pScientific.Arrays.BaseArray2D cimport BaseArray2D  , \
                                            CView2D
from pScientific.Arrays.RealArray1D cimport CRealArray1D , \
                                            RealArray1D
from pScientific.Arrays.Slicing     cimport CMultiSlice

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "RealArray2D.h":

    # . The array type.
    ctypedef struct CRealArray2D "RealArray2D":
        pass

    # . Functions - general.
    cdef CRealArray2D *RealArray2D_Allocate                 ( CStatus       *status            )
    cdef CRealArray2D *RealArray2D_AllocateWithExtents      ( CInteger       rows              ,
                                                              CInteger       columns           ,
                                                              CStatus       *status            )
    cdef void          RealArray2D_AssignBlock              ( CRealArray2D  *self              ,
                                                              CRealBlock    *block             ,
                                                              CBoolean       withReference     ,
                                                              CStatus       *status            )
    cdef void          RealArray2D_CopyTo                   ( CRealArray2D  *self              ,
                                                              CRealArray2D  *other             ,
                                                              CStatus       *status            )
    cdef void          RealArray2D_Deallocate               ( CRealArray2D **self              )
    cdef CReal         RealArray2D_GetItem                  ( CRealArray2D  *self              ,
                                                              CInteger       i                 ,
                                                              CInteger       j                 ,
                                                              CStatus       *status            )
    cdef CReal         RealArray2D_GetItemMultiSlice        ( CRealArray2D  *self              ,
                                                              CMultiSlice   *multiSlice        ,
                                                              CStatus       *status            )
    cdef CReal        *RealArray2D_PointerToData            ( CRealArray2D  *self              )
    cdef void          RealArray2D_Set                      ( CRealArray2D  *self              ,
                                                              CReal          value             )
    cdef void          RealArray2D_SetItem                  ( CRealArray2D  *self              ,
                                                              CInteger       i                 ,
                                                              CInteger       j                 ,
                                                              CReal          value             ,
                                                              CStatus       *status            )
    cdef void          RealArray2D_SetItemMultiSlice        ( CRealArray2D  *self              ,
                                                              CMultiSlice   *multiSlice        ,
                                                              CReal          value             ,
                                                              CStatus       *status            )

    # . Functions - specific.
    cdef void          RealArray2D_DiagonalOfProduct        ( CRealArray2D  *self              ,
                                                              CBoolean       sTranspose        ,
                                                              CRealArray2D  *other             ,
                                                              CBoolean       oTranspose        ,
                                                              CRealArray1D  *diagonal          ,
                                                              CStatus       *status            )
    cdef CInteger      RealArray2D_GramSchmidtOrthogonalize ( CRealArray2D  *self              ,
                                                              CInteger      *maximumIterations ,
                                                              CInteger      *numberConstant    ,
                                                              CReal         *tolerance         ,
                                                              CStatus       *status            )
    cdef CBoolean      RealArray2D_IsOrthogonal             ( CRealArray2D  *self              ,
                                                              CReal         *tolerance         ,
                                                              CReal         *deviation         ,
                                                              CStatus       *status            )
    cdef CBoolean      RealArray2D_IsSymmetric              ( CRealArray2D  *self              ,
                                                              CReal         *tolerance         ,
                                                              CReal         *deviation         )
    cdef void          RealArray2D_MatrixMultiply           ( CBoolean       aTranspose        ,
                                                              CBoolean       bTranspose        ,
                                                              CReal          alpha             ,
                                                              CRealArray2D  *a                 ,
                                                              CRealArray2D  *b                 ,
                                                              CReal          beta              ,
                                                              CRealArray2D  *c                 ,
                                                              CStatus       *status            )
    cdef void          RealArray2D_ProjectOutOfArray1D      ( CRealArray2D  *self              ,
                                                              CRealArray1D  *vector            ,
                                                              CStatus       *status            )
    cdef CReal         RealArray2D_Trace                    ( CRealArray2D  *self              ,
                                                              CStatus       *status            )
    cdef CReal         RealArray2D_TraceOfProduct           ( CRealArray2D  *self              ,
                                                              CRealArray2D  *other             ,
                                                              CStatus       *status            )
    cdef void          RealArray2D_TransposeSquare          ( CRealArray2D  *self              ,
                                                              CStatus       *status            )
    cdef void          RealArray2D_VectorMultiply           ( CBoolean       aTranspose        ,
                                                              CReal          alpha             ,
                                                              CRealArray2D  *a                 ,
                                                              CRealArray1D  *x                 ,
                                                              CReal          beta              ,
                                                              CRealArray1D  *y                 ,
                                                              CStatus       *status            )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class RealArray2D ( BaseArray2D ):

    cdef CRealArray2D *cObject
