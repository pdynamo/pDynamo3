from pCore.BooleanBlock              cimport BooleanBlock
from pCore.CPrimitiveTypes           cimport CBoolean            , \
                                             CFalse              , \
                                             CInteger            , \
                                             CReal               , \
                                             CTrue
from pCore.Status                    cimport CStatus             , \
                                             CStatus_OK
from pScientific.Arrays.BaseIterator cimport BaseIterator        , \
                                             CIterator           , \
                                             Iterator_DataOffSet , \
                                             Iterator_NextIndex

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "BooleanIterator.h":

    # . Functions.
    # . Unary.
    cdef CBoolean BooleanIterator_All        ( CIterator *self      ,
                                               CBoolean  *selfData  )
    cdef CBoolean BooleanIterator_Any        ( CIterator *self      ,
                                               CBoolean  *selfData  )
    cdef CInteger BooleanIterator_CountFalse ( CIterator *self      ,
                                               CBoolean  *selfData  )
    cdef CInteger BooleanIterator_CountTrue  ( CIterator *self      ,
                                               CBoolean  *selfData  )
    cdef void     BooleanIterator_Not        ( CIterator *self      ,
                                               CBoolean  *selfData  ,
                                               CStatus   *status    )
    cdef void     BooleanIterator_Set        ( CIterator *self      ,
                                               CBoolean  *selfData  ,
                                               CBoolean   value     ,
                                               CStatus   *status    )
    # . Binary.
    cdef void     BooleanIterator_And        ( CIterator *self      ,
                                               CBoolean  *selfData  ,
                                               CIterator *other     ,
                                               CBoolean  *otherData ,
                                               CStatus   *status    )
    cdef void     BooleanIterator_CopyTo     ( CIterator *self      ,
                                               CBoolean  *selfData  ,
                                               CIterator *other     ,
                                               CBoolean  *otherData ,
                                               CStatus   *status    )
    cdef void     BooleanIterator_Or         ( CIterator *self      ,
                                               CBoolean  *selfData  ,
                                               CIterator *other     ,
                                               CBoolean  *otherData ,
                                               CStatus   *status    )
    cdef void     BooleanIterator_Xor        ( CIterator *self      ,
                                               CBoolean  *selfData  ,
                                               CIterator *other     ,
                                               CBoolean  *otherData ,
                                               CStatus   *status    )
    cdef void     BooleanIterator_Swap       ( CIterator *self      ,
                                               CBoolean  *selfData  ,
                                               CIterator *other     ,
                                               CBoolean  *otherData ,
                                               CStatus   *status    )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class BooleanIterator ( BaseIterator ):

    cdef CBoolean *cData
