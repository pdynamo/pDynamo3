from pCore.CPrimitiveTypes    cimport CBoolean, CFalse, CInteger, CReal, CTrue
from pCore.IntegerUtilities   cimport Integer_Allocate, Integer_Deallocate
from pCore.Selection          cimport CSelection, Selection
from pCore.SelectionContainer cimport CSelectionContainer, SelectionContainer
from pCore.Status             cimport CStatus, CStatus_OK

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "PairList.h":

    ctypedef struct CPairList "PairList":
        CBoolean isSelf
        CBoolean isSorted

    ctypedef struct CPairListIterator "PairListIterator":
        pass

    ctypedef struct CPairRecord "PairRecord":
        CInteger  index
        CInteger  capacity
        CInteger *indices

    cdef CPairList           *PairList_Allocate                    ( CInteger           capacity      ,
                                                                     CStatus           *status        )
    cdef void                 PairList_Append                      ( CPairList         *self          ,
                                                                     CPairRecord       *record        ,
                                                                     CStatus           *status        )
    cdef void                 PairList_Deallocate                  ( CPairList        **self          )
    cdef CInteger             PairList_NumberOfPairs               ( CPairList         *self          )

    cdef CPairList           *CrossPairList_MakeFull               ( CInteger           capacity1     ,
                                                                     CSelection        *andSelection1 ,
                                                                     CInteger           capacity2     ,
                                                                     CSelection        *andSelection2 ,
                                                                     CStatus           *status        )
    cdef CPairList           *CrossPairList_MakeFullExcluded       ( CInteger           capacity1     ,
                                                                     CSelection        *andSelection1 ,
                                                                     CInteger           capacity2     ,
                                                                     CSelection        *andSelection2 ,
                                                                     CStatus           *status        )

    cdef CSelectionContainer *SelfPairList_GetConnectedComponents  ( CPairList         *self          ,
                                                                     CInteger           upperBound    ,
                                                                     CStatus           *status        )
    cdef void                *SelfPairList_Renumber                ( CPairList         *self          ,
                                                                     CSelection        *mapping       ,
                                                                     CStatus           *status        )
    cdef CPairList           *SelfPairList_ToCrossPairList         ( CPairList         *self          ,
                                                                     CSelection        *andSelection1 ,
                                                                     CSelection        *andSelection2 ,
                                                                     CSelection        *orSelection   ,
                                                                     CStatus           *status        )
    cdef CPairList           *SelfPairList_ToCrossPairListExcluded ( CPairList         *self          ,
                                                                     CSelection        *andSelection1 ,
                                                                     CSelection        *andSelection2 ,
                                                                     CSelection        *orSelection   ,
                                                                     CStatus           *status        )
    cdef CPairList           *SelfPairList_ToSelfPairList          ( CPairList         *self          ,
                                                                     CSelection        *andSelection  ,
                                                                     CSelection        *orSelection   ,
                                                                     CStatus           *status        )
    cdef CPairList           *SelfPairList_ToSelfPairListExcluded  ( CPairList         *self          ,
                                                                     CInteger           capacity      ,
                                                                     CSelection        *andSelection  ,
                                                                     CSelection        *orSelection   ,
                                                                     CStatus           *status        )

    cdef void                 PairListIterator_Initialize          ( CPairListIterator *self          ,
                                                                     CPairList         *target        )
    cdef CPairRecord         *PairListIterator_Next                ( CPairListIterator *self          )

    cdef CPairRecord         *PairRecord_FromIndices               ( CInteger           index         ,
                                                                     CInteger           capacity      ,
                                                                     CInteger          *indices       ,
                                                                     CStatus           *status        )

#===================================================================================================================================
# . Class and subclasses.
#===================================================================================================================================
cdef class PairList:

    cdef CPairList     *cObject
    cdef public object  label
    cdef public object  isOwner

cdef class CrossPairList ( PairList ):
    pass

cdef class SelfPairList ( PairList ):
    pass

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class PairListIterator:

    cdef CInteger           current
    cdef CPairListIterator  cIterator
    cdef CPairRecord       *cRecord
