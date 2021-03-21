from pCore.CPrimitiveTypes              cimport CBoolean, CFalse, CInteger, CReal, CTrue
from pCore.Selection                    cimport Selection, CSelection
from pMolecule.MMModel.MMTerm           cimport MMTerm
from pScientific.Arrays.RealArray2D     cimport CRealArray2D
from pScientific.Geometry3.Coordinates3 cimport Coordinates3

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "HarmonicAngleContainer.h":

    ctypedef struct CHarmonicAngle "HarmonicAngle":
        CBoolean isActive
        CInteger atom1
        CInteger atom2
        CInteger atom3
        CInteger type

    ctypedef struct CHarmonicAngleParameter "HarmonicAngleParameter":
        CReal eq
        CReal fc

    ctypedef struct CHarmonicAngleContainer "HarmonicAngleContainer":
       CInteger                 nParameters
       CInteger                 nTerms
       CHarmonicAngle          *terms
       CHarmonicAngleParameter *parameters

    cdef void                     HarmonicAngleContainer_ActivateTerms         ( CHarmonicAngleContainer  *self )
    cdef CHarmonicAngleContainer *HarmonicAngleContainer_Allocate              ( CInteger nTerms, CInteger nParameters )
    cdef CHarmonicAngleContainer *HarmonicAngleContainer_Clone                 ( CHarmonicAngleContainer  *self )
    cdef void                     HarmonicAngleContainer_DeactivateTerms       ( CHarmonicAngleContainer  *self, CSelection *selection )
    cdef void                     HarmonicAngleContainer_Deallocate            ( CHarmonicAngleContainer **self )
    cdef CReal                    HarmonicAngleContainer_Energy                ( CHarmonicAngleContainer  *self, CRealArray2D *coordinates3, CRealArray2D *gradients3 )
    cdef CHarmonicAngleContainer *HarmonicAngleContainer_Merge                 ( CHarmonicAngleContainer  *self, CHarmonicAngleContainer *other, CInteger atomincrement )
    cdef CInteger                 HarmonicAngleContainer_NumberOfInactiveTerms ( CHarmonicAngleContainer  *self )
    cdef CHarmonicAngleContainer *HarmonicAngleContainer_Prune                 ( CHarmonicAngleContainer  *self, CSelection *selection )
    cdef void                     HarmonicAngleContainer_Sort                  ( CHarmonicAngleContainer  *self )
    cdef CInteger                 HarmonicAngleContainer_UpperBound            ( CHarmonicAngleContainer  *self )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class HarmonicAngleContainer ( MMTerm ):

    cdef CHarmonicAngleContainer *cObject
    cdef public object            isOwner
    cdef public object            parameterKeys
