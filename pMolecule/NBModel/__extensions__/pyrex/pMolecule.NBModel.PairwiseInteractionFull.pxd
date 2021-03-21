from pCore.CPrimitiveTypes                  cimport CBoolean              , \
                                                    CFalse                , \
                                                    CInteger              , \
                                                    CReal                 , \
                                                    CTrue
from pCore.PairList                         cimport CPairList             , \
                                                    PairList
from pCore.Status                           cimport CStatus               , \
                                                    CStatus_OK
from pMolecule.MMModel.LJParameterContainer cimport CLJParameterContainer , \
                                                    LJParameterContainer
from pMolecule.NBModel.PairwiseInteraction  cimport PairwiseInteraction
from pScientific.Arrays.IntegerArray1D      cimport CIntegerArray1D       , \
                                                    IntegerArray1D
from pScientific.Arrays.RealArray1D         cimport CRealArray1D          , \
                                                    RealArray1D
from pScientific.Arrays.RealArray2D         cimport CRealArray2D
from pScientific.Geometry3.Coordinates3     cimport Coordinates3

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "PairwiseInteractionFull.h":

    ctypedef struct CPairwiseInteractionFull "PairwiseInteractionFull":
        CReal dampingCutOff

    cdef CPairwiseInteractionFull *PairwiseInteractionFull_Allocate            ( CStatus                   *status             )
    cdef CPairwiseInteractionFull *PairwiseInteractionFull_Clone               ( CPairwiseInteractionFull  *self               ,
                                                                                 CStatus                   *status             )
    cdef void                      PairwiseInteractionFull_Deallocate          ( CPairwiseInteractionFull **self               )
    cdef void                      PairwiseInteractionFull_InitializeDependent ( CPairwiseInteractionFull  *self               )
    cdef void                      PairwiseInteractionFull_Interactions        ( CPairwiseInteractionFull  *self               ,
                                                                                 CRealArray1D              *r                  ,
                                                                                 CRealArray1D              *lennardJonesA      ,
                                                                                 CRealArray1D              *lennardJonesB      ,
                                                                                 CRealArray1D              *multipole0         ,
                                                                                 CRealArray1D              *multipole1         ,
                                                                                 CRealArray1D              *multipole2         )
    cdef void                      PairwiseInteractionFull_MMMMEnergy          ( CPairwiseInteractionFull  *self               ,
                                                                                 CRealArray1D              *chargesI           ,
                                                                                 CRealArray1D              *chargesJ           ,
                                                                                 CIntegerArray1D           *ljTypesI           ,
                                                                                 CIntegerArray1D           *ljTypesJ           ,
                                                                                 CLJParameterContainer     *ljParameters       ,
                                                                                 CReal                      electrostaticScale ,
                                                                                 CReal                      lennardJonesScale  ,
                                                                                 CRealArray2D              *coordinates3I      ,
                                                                                 CRealArray2D              *coordinates3J      ,
                                                                                 CPairList                 *pairList           ,
                                                                                 CReal                     *eElectrostatic     ,
                                                                                 CReal                     *eLennardJones      ,
                                                                                 CRealArray2D              *gradients3I        ,
                                                                                 CRealArray2D              *gradients3J        ,
                                                                                 CStatus                   *status             )
    cdef void                      PairwiseInteractionFull_QCMMGradients       ( CPairwiseInteractionFull  *self               ,
                                                                                 CInteger                   multipoleOrder     ,
                                                                                 CRealArray1D              *multipolesQ        ,
                                                                                 CRealArray1D              *chargesM           ,
                                                                                 CReal                      electrostaticScale ,
                                                                                 CRealArray2D              *coordinates3Q      ,
                                                                                 CRealArray2D              *coordinates3M      ,
                                                                                 CPairList                 *pairList           ,
                                                                                 CRealArray2D              *gradients3Q        ,
                                                                                 CRealArray2D              *gradients3M        ,
                                                                                 CStatus                   *status             )
    cdef void                      PairwiseInteractionFull_QCMMPotentials      ( CPairwiseInteractionFull  *self               ,
                                                                                 CInteger                   multipoleOrder     ,
                                                                                 CRealArray1D              *chargesM           ,
                                                                                 CReal                      electrostaticScale ,
                                                                                 CRealArray2D              *coordinates3Q      ,
                                                                                 CRealArray2D              *coordinates3M      ,
                                                                                 CPairList                 *pairList           ,
                                                                                 CRealArray1D              *potentials         ,
                                                                                 CStatus                   *status             )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class PairwiseInteractionFull ( PairwiseInteraction ):

    cdef CPairwiseInteractionFull *cObject
    cdef public object             isOwner
