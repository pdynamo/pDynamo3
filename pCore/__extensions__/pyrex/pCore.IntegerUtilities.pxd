from pCore.CPrimitiveTypes cimport CInteger
from pCore.Status          cimport CStatus

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "IntegerUtilities.h":

    cdef CInteger *Integer_Allocate   ( CInteger   capacity ,
                                        CStatus   *status   )
    cdef CInteger *Integer_Clone      ( CInteger  *self     ,
                                        CInteger   capacity ,
                                        CStatus   *status   )
    cdef void      Integer_CopyTo     ( CInteger  *self     ,
                                        CInteger   capacity ,
                                        CInteger  *other    ,
                                        CStatus   *status   )
    cdef void      Integer_Deallocate ( CInteger **self     )
    cdef CInteger  Integer_GetItem    ( CInteger  *self     ,
                                        CInteger   capacity ,
                                        CInteger   index    ,
                                        CStatus   *status   )
    cdef void      Integer_Set        ( CInteger  *self     ,
                                        CInteger   capacity ,
                                        CInteger   value    )
    cdef void      Integer_SetItem    ( CInteger  *self     ,
                                        CInteger   capacity ,
                                        CInteger   index    ,
                                        CInteger   value    ,
                                        CStatus   *status   )
    cdef void      Integer_Sort       ( CInteger  *self     ,
                                        CInteger   capacity )
