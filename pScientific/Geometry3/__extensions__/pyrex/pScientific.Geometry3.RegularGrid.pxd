from pCore.CPrimitiveTypes cimport CBoolean, CFalse, CInteger, CReal, CTrue
from pCore.Status          cimport CStatus

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "RegularGrid.h":

    # . The regular grid dimension type.
    ctypedef struct CRegularGridDimension "RegularGridDimension":
        CBoolean isPeriodic
        CInteger bins
        CInteger stride
        CReal    binSize
        CReal    lower
        CReal    midPointLower
        CReal    period
        CReal    upper

    # . The regular grid type.
    ctypedef struct CRegularGrid "RegularGrid":
        CInteger               ndimensions
        CRegularGridDimension *dimensions

    # . Functions.
    # . Regular grid dimension.
    cdef void   RegularGridDimension_CopyTo     ( CRegularGridDimension *self, CRegularGridDimension *other )
    cdef void   RegularGridDimension_Initialize ( CRegularGridDimension *self )

    # . Regular grid.
    cdef CRegularGrid *RegularGrid_Allocate           ( CInteger ndimensions, CStatus *status )
    cdef CRegularGrid *RegularGrid_Clone              ( CRegularGrid  *self, CStatus *status )
    cdef void          RegularGrid_Deallocate         ( CRegularGrid **self )
    cdef CInteger      RegularGrid_NumberOfGridPoints ( CRegularGrid  *self )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class RegularGrid:

    cdef CRegularGrid  *cObject
    cdef public object  isOwner
