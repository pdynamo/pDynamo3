from pCore.CPrimitiveTypes              cimport CBoolean, CFalse, CInteger, CReal, CTrue
from pCore.Selection                    cimport Selection, CSelection
from pMolecule.MMModel.MMTerm           cimport MMTerm
from pScientific.Arrays.RealArray1D     cimport CRealArray1D, RealArray1D_AllocateWithExtent , RealArray1D_GetItem, RealArray1D_SetItem
from pScientific.Arrays.RealArray2D     cimport CRealArray2D, RealArray2D_AllocateWithExtents, RealArray2D_GetItem, RealArray2D_SetItem
from pScientific.Geometry3.Coordinates3 cimport Coordinates3
from pScientific.Splines.BicubicSpline  cimport BicubicSpline_MakeFromRealArray2D, BicubicSplineType_Periodic, CBicubicSpline

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "CMAPDihedralContainer.h":

    ctypedef struct CCMAPDihedral "CMAPDihedral":
        CBoolean isActive
        CInteger atom1
        CInteger atom2
        CInteger atom3
        CInteger atom4
        CInteger atom5
        CInteger atom6
        CInteger atom7
        CInteger atom8
        CInteger type

    ctypedef struct CCMAPDihedralContainer "CMAPDihedralContainer":
        CBoolean         isSorted
        CInteger         nParameters
        CInteger         nTerms
        CCMAPDihedral   *terms
        CBicubicSpline **parameters

    cdef void                    CMAPDihedralContainer_ActivateTerms         ( CCMAPDihedralContainer  *self )
    cdef CCMAPDihedralContainer *CMAPDihedralContainer_Allocate              ( CInteger nTerms, CInteger nParameters )
    cdef CCMAPDihedralContainer *CMAPDihedralContainer_Clone                 ( CCMAPDihedralContainer  *self )
    cdef void                    CMAPDihedralContainer_DeactivateTerms       ( CCMAPDihedralContainer  *self, CSelection *selection )
    cdef void                    CMAPDihedralContainer_Deallocate            ( CCMAPDihedralContainer **self )
    cdef CReal                   CMAPDihedralContainer_Energy                ( CCMAPDihedralContainer  *self, CRealArray2D *coordinates3, CRealArray2D *gradients3 )
    cdef CCMAPDihedralContainer *CMAPDihedralContainer_Merge                 ( CCMAPDihedralContainer  *self, CCMAPDihedralContainer *other, CInteger atomincrement )
    cdef CInteger                CMAPDihedralContainer_NumberOfInactiveTerms ( CCMAPDihedralContainer  *self )
    cdef CCMAPDihedralContainer *CMAPDihedralContainer_Prune                 ( CCMAPDihedralContainer  *self, CSelection *selection )
    cdef void                    CMAPDihedralContainer_Sort                  ( CCMAPDihedralContainer  *self )
    cdef CInteger                CMAPDihedralContainer_UpperBound            ( CCMAPDihedralContainer  *self )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class CMAPDihedralContainer ( MMTerm ):

    cdef CCMAPDihedralContainer *cObject
    cdef public object           isOwner
    cdef public object           parameterKeys
