from pCore.CPrimitiveTypes           cimport CBoolean            , \
                                             CFalse              , \
                                             CInteger            , \
                                             CReal               , \
                                             CTrue
from pCore.IntegerBlock              cimport IntegerBlock
from pCore.Status                    cimport CStatus             , \
                                             CStatus_OK
from pScientific.Arrays.BaseIterator cimport BaseIterator        , \
                                             CIterator           , \
                                             Iterator_DataOffSet , \
                                             Iterator_NextIndex

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "IntegerIterator.h":

    # . Functions.
    # . Unary.
    cdef CInteger IntegerIterator_AbsoluteMaximum   ( CIterator *self      ,
                                                      CInteger  *selfData  )
    cdef CInteger IntegerIterator_CountSmall        ( CIterator *self      ,
                                                      CInteger  *selfData  ,
                                                      CInteger   tolerance )
    cdef CInteger  IntegerIterator_DotSelf          ( CIterator *self      ,
                                                      CInteger  *selfData  )
    cdef void     IntegerIterator_FilterGreaterThan ( CIterator *self      ,
                                                      CInteger  *selfData  ,
                                                      CInteger   tolerance ,
                                                      CInteger   value     ,
                                                      CStatus   *status    )
    cdef void     IntegerIterator_FilterLessThan    ( CIterator *self      ,
                                                      CInteger  *selfData  ,
                                                      CInteger   tolerance ,
                                                      CInteger   value     ,
                                                      CStatus   *status    )
    cdef void     IntegerIterator_FilterSmall       ( CIterator *self      ,
                                                      CInteger  *selfData  ,
                                                      CInteger   tolerance ,
                                                      CInteger   value     ,
                                                      CStatus   *status    )
    cdef void     IntegerIterator_Increment         ( CIterator *self      ,
                                                      CInteger  *selfData  ,
                                                      CInteger   value     ,
                                                      CStatus   *status    )
    cdef CInteger IntegerIterator_Maximum           ( CIterator *self      ,
                                                      CInteger  *selfData  )
    cdef CInteger IntegerIterator_Minimum           ( CIterator *self      ,
                                                      CInteger  *selfData  )
    cdef CInteger IntegerIterator_Product           ( CIterator *self      ,
                                                      CInteger  *selfData  )
    cdef void     IntegerIterator_Scale             ( CIterator *self      ,
                                                      CInteger  *selfData  ,
                                                      CInteger   value     ,
                                                      CStatus   *status    )
    cdef void     IntegerIterator_Set               ( CIterator *self      ,
                                                      CInteger  *selfData  ,
                                                      CInteger   value     ,
                                                      CStatus   *status    )
    cdef CReal    IntegerIterator_Sparsity          ( CIterator *self      ,
                                                      CInteger  *selfData  ,
                                                      CInteger   tolerance )
    cdef CInteger IntegerIterator_Sum               ( CIterator *self      ,
                                                      CInteger  *selfData  )

    # . Binary.
    cdef void     IntegerIterator_Add               ( CIterator *self      ,
                                                      CInteger  *selfData  ,
                                                      CIterator *other     ,
                                                      CInteger  *otherData ,
                                                      CInteger   scale     ,
                                                      CStatus   *status    )
    cdef void     IntegerIterator_CopyTo            ( CIterator *self      ,
                                                      CInteger  *selfData  ,
                                                      CIterator *other     ,
                                                      CInteger  *otherData ,
                                                      CStatus   *status    )
    cdef CInteger IntegerIterator_Dot               ( CIterator *self      ,
                                                      CInteger  *selfData  ,
                                                      CIterator *other     ,
                                                      CInteger  *otherData ,
                                                      CStatus   *status    )
    cdef void     IntegerIterator_Multiply          ( CIterator *self      ,
                                                      CInteger  *selfData  ,
                                                      CIterator *other     ,
                                                      CInteger  *otherData ,
                                                      CStatus   *status    )
    cdef void     IntegerIterator_Swap              ( CIterator *self      ,
                                                      CInteger  *selfData  ,
                                                      CIterator *other     ,
                                                      CInteger  *otherData ,
                                                      CStatus   *status    )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class IntegerIterator ( BaseIterator ):

    cdef CInteger *cData
