from pCore.CPrimitiveTypes           cimport CBoolean            , \
                                             CFalse              , \
                                             CInteger            , \
                                             CReal               , \
                                             CTrue
from pCore.RealBlock                 cimport RealBlock
from pCore.Status                    cimport CStatus             , \
                                             CStatus_OK
from pScientific.Arrays.BaseIterator cimport BaseIterator        , \
                                             CIterator           , \
                                             Iterator_DataOffSet , \
                                             Iterator_NextIndex

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "RealIterator.h":

    # . Functions.
    # . Unary.
    cdef void     RealIterator_Absolute            ( CIterator *self      ,
                                                     CReal     *selfData  )
    cdef CReal    RealIterator_AbsoluteMaximum     ( CIterator *self      ,
                                                     CReal     *selfData  )
    cdef CInteger RealIterator_CountSmall          ( CIterator *self      ,
                                                     CReal     *selfData  ,
                                                     CReal      tolerance )
    cdef CReal    RealIterator_DotSelf             ( CIterator *self      ,
                                                     CReal     *selfData  )
    cdef void     RealIterator_Exponential         ( CIterator *self      ,
                                                     CReal     *selfData  ,
                                                     CStatus   *status    )
    cdef void     RealIterator_FilterGreaterThan   ( CIterator *self      ,
                                                     CReal     *selfData  ,
                                                     CReal      tolerance ,
                                                     CReal      value     ,
                                                     CStatus   *status    )
    cdef void     RealIterator_FilterLessThan      ( CIterator *self      ,
                                                     CReal     *selfData  ,
                                                     CReal      tolerance ,
                                                     CReal      value     ,
                                                     CStatus   *status    )
    cdef void     RealIterator_FilterSmall         ( CIterator *self      ,
                                                     CReal     *selfData  ,
                                                     CReal      tolerance ,
                                                     CReal      value     ,
                                                     CStatus   *status    )
    cdef void     RealIterator_Increment           ( CIterator *self      ,
                                                     CReal     *selfData  ,
                                                     CReal      value     ,
                                                     CStatus   *status    )
    cdef CReal    RealIterator_Maximum             ( CIterator *self      ,
                                                     CReal     *selfData  )
    cdef CReal    RealIterator_Minimum             ( CIterator *self      ,
                                                     CReal     *selfData  )
    cdef void     RealIterator_NaturalLogarithm    ( CIterator *self      ,
                                                     CReal     *selfData  ,
                                                     CStatus   *status    )
    cdef CReal    RealIterator_Norm2               ( CIterator *self      ,
                                                     CReal     *selfData  )
    cdef void     RealIterator_Normalize           ( CIterator *self      ,
                                                     CReal     *selfData  ,
                                                     CReal      tolerance ,
                                                     CStatus   *status    )
    cdef void     RealIterator_Power               ( CIterator *self      ,
                                                     CReal     *selfData  ,
                                                     CReal      power     ,
                                                     CStatus   *status    )
    cdef CReal    RealIterator_Product             ( CIterator *self      ,
                                                     CReal     *selfData  )
    cdef void     RealIterator_Reciprocate         ( CIterator *self      ,
                                                     CReal     *selfData  ,
                                                     CReal      tolerance ,
                                                     CReal      value     ,
                                                     CStatus   *status    )
    cdef void     RealIterator_ReciprocatePower    ( CIterator *self      ,
                                                     CReal     *selfData  ,
                                                     CReal      power     ,
                                                     CReal      tolerance ,
                                                     CReal      value     ,
                                                     CStatus   *status    )
    cdef CReal    RealIterator_RootMeanSquare      ( CIterator *self      ,
                                                     CReal     *selfData  )
    cdef void     RealIterator_Scale               ( CIterator *self      ,
                                                     CReal     *selfData  ,
                                                     CReal      value     ,
                                                     CStatus   *status    )
    cdef void     RealIterator_Set                 ( CIterator *self      ,
                                                     CReal     *selfData  ,
                                                     CReal      value     ,
                                                     CStatus   *status    )
    cdef CReal    RealIterator_Sparsity            ( CIterator *self      ,
                                                     CReal     *selfData  ,
                                                     CReal      tolerance )
    cdef void     RealIterator_Square              ( CIterator *self      ,
                                                     CReal     *selfData  ,
                                                     CStatus   *status    )
    cdef void     RealIterator_SquareRoot          ( CIterator *self      ,
                                                     CReal     *selfData  ,
                                                     CStatus   *status    )
    cdef CReal    RealIterator_Sum                 ( CIterator *self      ,
                                                     CReal     *selfData  )

    # . Binary.
    cdef void     RealIterator_Add                 ( CIterator *self      ,
                                                     CReal     *selfData  ,
                                                     CIterator *other     ,
                                                     CReal     *otherData ,
                                                     CReal      scale     ,
                                                     CStatus   *status    )
    cdef void     RealIterator_CopyTo              ( CIterator *self      ,
                                                     CReal     *selfData  ,
                                                     CIterator *other     ,
                                                     CReal     *otherData ,
                                                     CStatus   *status    )
    cdef void     RealIterator_Divide              ( CIterator *self      ,
                                                     CReal     *selfData  ,
                                                     CIterator *other     ,
                                                     CReal     *otherData ,
                                                     CReal      tolerance ,
                                                     CReal      value     ,
                                                     CStatus   *status    )
    cdef CReal    RealIterator_Dot                 ( CIterator *self      ,
                                                     CReal     *selfData  ,
                                                     CIterator *other     ,
                                                     CReal     *otherData ,
                                                     CStatus   *status    )
    cdef void     RealIterator_Multiply            ( CIterator *self      ,
                                                     CReal     *selfData  ,
                                                     CIterator *other     ,
                                                     CReal     *otherData ,
                                                     CStatus   *status    )
    cdef void     RealIterator_Swap                ( CIterator *self      ,
                                                     CReal     *selfData  ,
                                                     CIterator *other     ,
                                                     CReal     *otherData ,
                                                     CStatus   *status    )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class RealIterator ( BaseIterator ):

    cdef CReal *cData
