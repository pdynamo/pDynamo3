from pCore.BooleanBlock    cimport BooleanBlock, CBooleanBlock
from pCore.CPrimitiveTypes cimport CBoolean, CFalse, CInteger, CReal, CTrue
from pCore.Selection       cimport CSelection, Selection, Selection_Allocate
from pCore.Status          cimport CStatus, CStatus_OK

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "SelectionContainer.h":

    ctypedef struct CSelectionContainer "SelectionContainer":
        CInteger     capacity
        CSelection **items

    cdef CSelectionContainer *SelectionContainer_Allocate            ( CInteger              capacity ,
                                                                       CStatus              *status   )
    cdef CInteger             SelectionContainer_Capacity            ( CSelectionContainer  *self     )
    cdef CSelectionContainer *SelectionContainer_Clone               ( CSelectionContainer  *self     ,
                                                                       CStatus              *status   )
    cdef void                 SelectionContainer_Deallocate          ( CSelectionContainer **self     )
    cdef CSelectionContainer *SelectionContainer_FromCapacity        ( CInteger              capacity ,
                                                                       CStatus              *status   )
    cdef void                 SelectionContainer_FuseItems           ( CSelectionContainer  *self     ,
                                                                       CBooleanBlock        *toFuse   ,
                                                                       CStatus              *status   )
    cdef CBooleanBlock       *SelectionContainer_MakeMembershipFlags ( CSelectionContainer  *self     ,
                                                                       CSelection           *members  ,
                                                                       CBoolean              andTest  ,
                                                                       CStatus              *status   )
    cdef CSelection          *SelectionContainer_UnionOfItems        ( CSelectionContainer  *self     ,
                                                                       CSelection           *toUnion  ,
                                                                       CStatus              *status   )
    cdef CInteger             SelectionContainer_UpperBound          ( CSelectionContainer  *self     )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class SelectionContainer:

    cdef CSelectionContainer *cObject
    cdef public object        itemName
    cdef public object        isOwner
