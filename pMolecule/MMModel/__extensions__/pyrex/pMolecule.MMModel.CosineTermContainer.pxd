from pCore.CPrimitiveTypes              cimport CBoolean, CFalse, CInteger, CReal, CTrue
from pCore.Selection                    cimport Selection, CSelection
from pMolecule.MMModel.MMTerm           cimport MMTerm
from pScientific.Arrays.RealArray2D     cimport CRealArray2D
from pScientific.Geometry3.Coordinates3 cimport Coordinates3

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "CosineParameter.h":

    ctypedef struct CCosineParameter "CosineParameter":
        CInteger  nPowers
        CInteger  nTerms
        CInteger *periods
        CReal    *powerCoefficients
        CReal    *termCoefficients

    cdef void CosineParameter_Allocate   ( CCosineParameter *self   ,
                                           CInteger          nTerms )
    cdef void CosineParameter_Clone      ( CCosineParameter *self   ,
                                           CCosineParameter *other  )
    cdef void CosineParameter_Deallocate ( CCosineParameter *self   )
    cdef void CosineParameter_Initialize ( CCosineParameter *self   )
    cdef void CosineParameter_MakePowers ( CCosineParameter *self   )

cdef extern from "CosineTerm.h":

    ctypedef struct CCosineTerm "CosineTerm":
        CBoolean  isActive
        CInteger  type
        CInteger *indices

    cdef void CosineTerm_Allocate   ( CCosineTerm *self     ,
                                      CInteger     nIndices )
    cdef void CosineTerm_Clone      ( CCosineTerm *self     ,
                                      CCosineTerm *other    ,
                                      CInteger     nIndices )
    cdef void CosineTerm_Deallocate ( CCosineTerm *self     )
    cdef void CosineTerm_Initialize ( CCosineTerm *self     )

cdef extern from "CosineTermContainer.h":

    ctypedef struct CCosineTermContainer "CosineTermContainer":
        CInteger          nIndices
        CInteger          nParameters
        CInteger          nTerms
        CCosineParameter *parameters
        CCosineTerm      *terms

    cdef void                  CosineTermContainer_ActivateTerms         ( CCosineTermContainer  *self          )
    cdef CCosineTermContainer *CosineTermContainer_Allocate              ( CInteger               nIndices      ,
                                                                           CInteger               nTerms        ,
                                                                           CInteger               nParameters   )
    cdef CCosineTermContainer *CosineTermContainer_Clone                 ( CCosineTermContainer  *self          )
    cdef void                  CosineTermContainer_DeactivateTerms       ( CCosineTermContainer  *self          ,
                                                                           CSelection            *selection     )
    cdef void                  CosineTermContainer_Deallocate            ( CCosineTermContainer **self          )
    cdef CInteger              CosineTermContainer_FindMaximumPeriod     ( CCosineTermContainer  *self          )
    cdef void                  CosineTermContainer_MakePowers            ( CCosineTermContainer  *self          )
    cdef CInteger              CosineTermContainer_NumberOfInactiveTerms ( CCosineTermContainer  *self          )
    cdef CCosineTermContainer *CosineTermContainer_Prune                 ( CCosineTermContainer  *self          ,
                                                                           CSelection            *selection     )
    cdef CInteger              CosineTermContainer_UpperBound            ( CCosineTermContainer  *self          )

cdef extern from "CosineTermEnergies.h":

    cdef CReal CosineTermEnergy_Angle      ( CCosineTermContainer *self         ,
                                             CRealArray2D         *coordinates3 ,
                                             CRealArray2D         *gradients3   )
    cdef CReal CosineTermEnergy_Dihedral   ( CCosineTermContainer *self         ,
                                             CRealArray2D         *coordinates3 ,
                                             CRealArray2D         *gradients3   )
    cdef CReal CosineTermEnergy_OutOfPlane ( CCosineTermContainer *self         ,
                                             CRealArray2D         *coordinates3 ,
                                             CRealArray2D         *gradients3   )

#===================================================================================================================================
# . Classes.
#===================================================================================================================================
cdef class CosineTermContainer ( MMTerm ):

    cdef CCosineTermContainer *cObject
    cdef public object         isOwner
    cdef public object         parameterKeys

cdef class CosineAngleContainer ( CosineTermContainer ):
    pass

cdef class CosineDihedralContainer ( CosineTermContainer ):
    pass

cdef class CosineOutOfPlaneContainer ( CosineTermContainer ):
    pass
