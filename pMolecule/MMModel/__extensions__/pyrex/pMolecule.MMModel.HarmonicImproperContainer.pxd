from pCore.CPrimitiveTypes              cimport CBoolean, CFalse, CInteger, CReal, CTrue
from pCore.Selection                    cimport Selection, CSelection
from pMolecule.MMModel.MMTerm           cimport MMTerm
from pScientific.Arrays.RealArray2D     cimport CRealArray2D
from pScientific.Geometry3.Coordinates3 cimport Coordinates3

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "HarmonicImproperContainer.h":

    ctypedef struct CHarmonicImproper "HarmonicImproper":
        CBoolean isActive
        CInteger atom1
        CInteger atom2
        CInteger atom3
        CInteger atom4
        CInteger type

    ctypedef struct CHarmonicImproperParameter "HarmonicImproperParameter":
         CReal eq
         CReal fc
         CReal coseq
         CReal sineq

    ctypedef struct CHarmonicImproperContainer "HarmonicImproperContainer":
       CInteger                    nParameters
       CInteger                    nTerms
       CHarmonicImproper          *terms
       CHarmonicImproperParameter *parameters

    cdef void                        HarmonicImproperContainer_ActivateTerms         ( CHarmonicImproperContainer  *self )
    cdef CHarmonicImproperContainer *HarmonicImproperContainer_Allocate              ( CInteger nTerms, CInteger nParameters )
    cdef CHarmonicImproperContainer *HarmonicImproperContainer_Clone                 ( CHarmonicImproperContainer  *self )
    cdef void                        HarmonicImproperContainer_DeactivateTerms       ( CHarmonicImproperContainer  *self, CSelection *selection )
    cdef void                        HarmonicImproperContainer_Deallocate            ( CHarmonicImproperContainer **self )
    cdef void                        HarmonicImproperContainer_FillCosSinValues      ( CHarmonicImproperContainer  *self )
    cdef CReal                       HarmonicImproperContainer_Energy                ( CHarmonicImproperContainer  *self, CRealArray2D *coordinates3, CRealArray2D *gradients3 )
    cdef CHarmonicImproperContainer *HarmonicImproperContainer_Merge                 ( CHarmonicImproperContainer  *self, CHarmonicImproperContainer *other, CInteger atomincrement )
    cdef CInteger                    HarmonicImproperContainer_NumberOfInactiveTerms ( CHarmonicImproperContainer  *self )
    cdef CHarmonicImproperContainer *HarmonicImproperContainer_Prune                 ( CHarmonicImproperContainer  *self, CSelection *selection )
    cdef void                        HarmonicImproperContainer_Sort                  ( CHarmonicImproperContainer  *self )
    cdef CInteger                    HarmonicImproperContainer_UpperBound            ( CHarmonicImproperContainer  *self )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class HarmonicImproperContainer ( MMTerm ):

    cdef CHarmonicImproperContainer *cObject
    cdef public object               isOwner
    cdef public object               parameterKeys
