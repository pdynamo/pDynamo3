#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "RegularGridOccupancy.h":

    ctypedef struct CRegularGridOccupancy "RegularGridOccupancy":
        pass

    cdef void RegularGridOccupancy_Deallocate ( CRegularGridOccupancy **self )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class RegularGridOccupancy:

    cdef CRegularGridOccupancy *cObject
    cdef public object          isOwner
