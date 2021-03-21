from pCore.BooleanBlock    cimport CBooleanBlock
from pCore.CPrimitiveTypes cimport CBoolean, CFalse, CInteger, CReal, CTrue
from pCore.Status          cimport CStatus, CStatus_OK

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "Selection.h":

    ctypedef struct CSelection "Selection":
        CInteger  capacity
        CInteger *indices

    cdef CSelection    *Selection_Allocate            ( CInteger     capacity   ,
                                                        CStatus     *status     )
    cdef CInteger       Selection_Capacity            ( CSelection  *self       )
    cdef CSelection    *Selection_Clone               ( CSelection  *self       ,
                                                        CStatus     *status     )
    cdef CSelection    *Selection_Complement          ( CSelection  *self       ,
                                                        CInteger     upperBound ,
                                                        CStatus     *status     )
    cdef void           Selection_Deallocate          ( CSelection **self       )
    cdef CSelection    *Selection_Difference          ( CSelection  *self       ,
                                                        CInteger     number     ,
                                                        CSelection **others     ,
                                                        CStatus     *status     )
    cdef CBoolean       Selection_HasItem             ( CSelection  *self       ,
                                                        CInteger     value      ,
                                                        CStatus     *status     )
    cdef CSelection    *Selection_Intersection        ( CInteger     number     ,
                                                        CSelection **others     ,
                                                        CStatus     *status     )
    cdef CBooleanBlock *Selection_MakeFlags           ( CSelection  *self       ,
                                                        CInteger     upperBound ,
                                                        CStatus     *status     )
    cdef CInteger       Selection_PositionOfItem      ( CSelection  *self       ,
                                                        CInteger     value      ,
                                                        CStatus     *status     )
    cdef CSelection    *Selection_Prune               ( CSelection  *self       ,
                                                        CSelection  *toKeep     ,
                                                        CStatus     *status     )
    cdef CSelection    *Selection_SymmetricDifference ( CInteger     number     ,
                                                        CSelection **others     ,
                                                        CStatus     *status     )
    cdef CSelection    *Selection_Union               ( CInteger     number     ,
                                                        CSelection **others     ,
                                                        CStatus     *status     )
    cdef CInteger       Selection_UpperBound          ( CSelection  *self       )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class Selection:

    cdef CSelection   *cObject
    cdef public object isOwner
    cdef public object label
