from pCore.CPrimitiveTypes              cimport CBoolean, CFalse, CInteger, CReal, CTrue
from pCore.Selection                    cimport Selection, CSelection
from pMolecule.MMModel.MMTerm           cimport MMTerm
from pScientific.Arrays.RealArray2D     cimport CRealArray2D
from pScientific.Geometry3.Coordinates3 cimport Coordinates3

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "FourierDihedralContainer.h":

    ctypedef struct CFourierDihedral "FourierDihedral":
        CBoolean isActive
        CInteger atom1
        CInteger atom2
        CInteger atom3
        CInteger atom4
        CInteger type

    ctypedef struct CFourierDihedralParameter "FourierDihedralParameter":
         CInteger period
         CReal    fc
         CReal    phase
         CReal    cosphase
         CReal    sinphase

    ctypedef struct CFourierDihedralContainer "FourierDihedralContainer":
       CInteger                   nParameters
       CInteger                   nTerms
       CFourierDihedral          *terms
       CFourierDihedralParameter *parameters

    cdef void                       FourierDihedralContainer_ActivateTerms         ( CFourierDihedralContainer  *self )
    cdef CFourierDihedralContainer *FourierDihedralContainer_Allocate              ( CInteger nTerms, CInteger nParameters )
    cdef CFourierDihedralContainer *FourierDihedralContainer_Clone                 ( CFourierDihedralContainer  *self )
    cdef void                       FourierDihedralContainer_DeactivateTerms       ( CFourierDihedralContainer  *self, CSelection *selection )
    cdef void                       FourierDihedralContainer_Deallocate            ( CFourierDihedralContainer **self )
    cdef void                       FourierDihedralContainer_FillCosSinPhases      ( CFourierDihedralContainer  *self )
    cdef CReal                      FourierDihedralContainer_Energy                ( CFourierDihedralContainer  *self, CRealArray2D *coordinates3, CRealArray2D *gradients3 )
    cdef CFourierDihedralContainer *FourierDihedralContainer_Merge                 ( CFourierDihedralContainer  *self, CFourierDihedralContainer *other, CInteger atomincrement )
    cdef CInteger                   FourierDihedralContainer_NumberOfInactiveTerms ( CFourierDihedralContainer  *self )
    cdef CFourierDihedralContainer *FourierDihedralContainer_Prune                 ( CFourierDihedralContainer  *self, CSelection *selection )
    cdef void                       FourierDihedralContainer_Sort                  ( CFourierDihedralContainer  *self )
    cdef CInteger                   FourierDihedralContainer_UpperBound            ( CFourierDihedralContainer  *self )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class FourierDihedralContainer ( MMTerm ):

    cdef CFourierDihedralContainer *cObject
    cdef public object              isOwner
    cdef public object              parameterKeys
