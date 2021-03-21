from pCore.CPrimitiveTypes                    cimport CBoolean                , \
                                                      CFalse                  , \
                                                      CInteger                , \
                                                      CReal                   , \
                                                      CTrue
from pCore.Selection                          cimport CSelection              , \
                                                      Selection
from pCore.Status                             cimport CStatus                 , \
                                                      CStatus_OK
from pMolecule.NBModel.ImagePairListContainer cimport CImagePairListContainer , \
                                                      ImagePairListContainer
from pScientific.Arrays.RealArray2D           cimport CRealArray2D
from pScientific.Geometry3.Coordinates3       cimport Coordinates3
from pScientific.Symmetry.SymmetryParameters  cimport CSymmetryParameters     , \
                                                      SymmetryParameters

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "UpdateChecker.h":

    ctypedef struct CUpdateChecker "UpdateChecker":
        CReal buffer

    cdef CUpdateChecker *UpdateChecker_Allocate            ( CStatus                 *status              )
    cdef CBoolean        UpdateChecker_CheckForImageUpdate ( CSymmetryParameters     *set1                ,
                                                             CSymmetryParameters     *set2                ,
                                                             CImagePairListContainer *images              ,
                                                             CReal                    buffer              ,
                                                             CReal                    maximumDisplacement )
    cdef CBoolean        UpdateChecker_CheckForUpdate      ( CRealArray2D            *set1                ,
                                                             CRealArray2D            *set2                ,
                                                             CSelection              *freeAtoms           ,
                                                             CReal                    buffer              ,
                                                             CReal                   *maximumDisplacement )
    cdef void            UpdateChecker_Deallocate          ( CUpdateChecker         **self                )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class UpdateChecker:

    cdef CUpdateChecker *cObject
    cdef public object   isOwner
