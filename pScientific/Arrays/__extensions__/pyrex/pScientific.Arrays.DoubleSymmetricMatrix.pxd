from pCore.CPrimitiveTypes           cimport CBoolean     , \
                                             CFalse       , \
                                             CInteger     , \
                                             CReal        , \
                                             CTrue
from pCore.RealBlock                 cimport CRealBlock   , \
                                             RealBlock
from pCore.Status                    cimport CStatus      , \
                                             CStatus_OK
from pScientific.Arrays.BaseIterator cimport CIterator
from pScientific.Arrays.RealIterator cimport RealIterator

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "DoubleSymmetricMatrix.h":

    ctypedef struct CDoubleSymmetricMatrix "DoubleSymmetricMatrix":
        CInteger    extent
        CInteger    extent01
        CInteger    size
        CRealBlock *block
        CReal      *data

    cdef CDoubleSymmetricMatrix *DoubleSymmetricMatrix_Allocate            ( CStatus                 *status        )
    cdef CDoubleSymmetricMatrix *DoubleSymmetricMatrix_AllocateWithExtent  ( CInteger                 extent        ,
                                                                             CStatus                 *status        )
    cdef CDoubleSymmetricMatrix *DoubleSymmetricMatrix_CloneDeep           ( CDoubleSymmetricMatrix  *self          ,
                                                                             CStatus                 *status        )
    cdef CDoubleSymmetricMatrix *DoubleSymmetricMatrix_CloneShallow        ( CDoubleSymmetricMatrix  *self          ,
                                                                             CStatus                 *status        )
    cdef void                    DoubleSymmetricMatrix_CopyTo              ( CDoubleSymmetricMatrix  *self          ,
                                                                             CDoubleSymmetricMatrix  *other         ,
                                                                             CStatus                 *status        )
    cdef void                    DoubleSymmetricMatrix_Deallocate          ( CDoubleSymmetricMatrix **self          )
    cdef CDoubleSymmetricMatrix *DoubleSymmetricMatrix_FromExtentBlock     ( CInteger                 extent        ,
                                                                             CRealBlock              *block         ,
                                                                             CBoolean                 withReference ,
                                                                             CStatus                 *status        )
    cdef CReal                   DoubleSymmetricMatrix_GetItem             ( CDoubleSymmetricMatrix  *self          ,
                                                                             CInteger                 i             ,
                                                                             CInteger                 j             ,
                                                                             CInteger                 k             ,
                                                                             CInteger                 l             ,
                                                                             CStatus                 *status        )
    cdef void                    DoubleSymmetricMatrix_IncrementItem       ( CDoubleSymmetricMatrix  *self          ,
                                                                             CInteger                 i             ,
                                                                             CInteger                 j             ,
                                                                             CInteger                 k             ,
                                                                             CInteger                 l             ,
                                                                             CReal                    value         ,
                                                                             CStatus                 *status        )
    cdef CInteger                DoubleSymmetricMatrix_Index               ( CInteger                 i             ,
                                                                             CInteger                 j             ,
                                                                             CInteger                 k             ,
                                                                             CInteger                 l             )
    cdef CIterator              *DoubleSymmetricMatrix_MakeIterator        ( CDoubleSymmetricMatrix  *self          ,
                                                                             CStatus                 *status        )
    cdef void                    DoubleSymmetricMatrix_Set                 ( CDoubleSymmetricMatrix  *self          ,
                                                                             CReal                    value         )
    cdef void                    DoubleSymmetricMatrix_SetItem             ( CDoubleSymmetricMatrix  *self          ,
                                                                             CInteger                 i             ,
                                                                             CInteger                 j             ,
                                                                             CInteger                 k             ,
                                                                             CInteger                 l             ,
                                                                             CReal                    value         ,
                                                                             CStatus                 *status        )
    cdef CInteger                DoubleSymmetricMatrix_ViewSize            ( CInteger                 extent        )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class DoubleSymmetricMatrix:

    cdef CDoubleSymmetricMatrix *cObject
    cdef public RealBlock        block
    cdef        RealIterator    _iterator
