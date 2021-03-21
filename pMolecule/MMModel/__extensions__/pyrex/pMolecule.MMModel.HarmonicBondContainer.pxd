from pCore.BooleanBlock                 cimport CBooleanBlock
from pCore.CPrimitiveTypes              cimport CBoolean, CFalse, CInteger, CReal, CTrue
from pCore.Selection                    cimport CSelection, Selection, Selection_MakeFlags
from pMolecule.MMModel.MMTerm           cimport MMTerm
from pScientific.Arrays.RealArray2D     cimport CRealArray2D
from pScientific.Geometry3.Coordinates3 cimport Coordinates3

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "HarmonicBondContainer.h":

    ctypedef struct CHarmonicBond "HarmonicBond":
        CBoolean isActive
        CInteger atom1
        CInteger atom2
        CInteger type

    ctypedef struct CHarmonicBondParameter "HarmonicBondParameter":
        CReal eq
        CReal fc

    ctypedef struct CHarmonicBondContainer "HarmonicBondContainer":
       CInteger                nParameters
       CInteger                nTerms
       CHarmonicBond          *terms
       CHarmonicBondParameter *parameters

    cdef void                    HarmonicBondContainer_ActivateTerms         ( CHarmonicBondContainer  *self )
    cdef CHarmonicBondContainer *HarmonicBondContainer_Allocate              ( CInteger nTerms, CInteger nParameters )
    cdef CHarmonicBondContainer *HarmonicBondContainer_Clone                 ( CHarmonicBondContainer  *self )
    cdef void                    HarmonicBondContainer_DeactivateTerms       ( CHarmonicBondContainer  *self, CSelection *selection )
    cdef void                    HarmonicBondContainer_Deallocate            ( CHarmonicBondContainer **self )
    cdef CReal                   HarmonicBondContainer_Energy                ( CHarmonicBondContainer  *self, CRealArray2D *coordinates3, CRealArray2D *gradients3 )
    cdef CInteger                HarmonicBondContainer_IdentifyBoundaryAtoms ( CHarmonicBondContainer  *self, CSelection *qcAtoms, CInteger **mmboundary, CInteger **qcpartners )
    cdef CHarmonicBondContainer *HarmonicBondContainer_Merge                 ( CHarmonicBondContainer  *self, CHarmonicBondContainer *other, CInteger atomincrement )
    cdef CInteger                HarmonicBondContainer_NumberOfInactiveTerms ( CHarmonicBondContainer  *self )
    cdef CHarmonicBondContainer *HarmonicBondContainer_Prune                 ( CHarmonicBondContainer  *self, CSelection *selection )
    cdef void                    HarmonicBondContainer_Sort                  ( CHarmonicBondContainer  *self )
    cdef CInteger                HarmonicBondContainer_UpperBound            ( CHarmonicBondContainer  *self )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class HarmonicBondContainer ( MMTerm ):

    cdef CHarmonicBondContainer *cObject
    cdef public object           isOwner
    cdef public object           is12Interaction
    cdef public object           parameterKeys
