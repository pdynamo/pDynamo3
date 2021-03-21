from pCore.BooleanBlock                      cimport BooleanBlock          , \
                                                     CBooleanBlock
from pCore.CPrimitiveTypes                   cimport CBoolean              , \
                                                     CFalse                , \
                                                     CInteger              , \
                                                     CReal                 , \
                                                     CTrue
from pCore.PairList                          cimport CPairList             , \
                                                     PairList
from pCore.SelectionContainer                cimport CSelectionContainer   , \
                                                     SelectionContainer
from pCore.Status                            cimport CStatus               , \
                                                     CStatus_OK
from pMolecule.MMModel.LJParameterContainer  cimport CLJParameterContainer , \
                                                     LJParameterContainer
from pMolecule.NBModel.PairwiseInteraction   cimport PairwiseInteraction
from pScientific.Arrays.IntegerArray1D       cimport CIntegerArray1D       , \
                                                     IntegerArray1D
from pScientific.Arrays.RealArray1D          cimport CRealArray1D          , \
                                                     RealArray1D
from pScientific.Arrays.RealArray2D          cimport CRealArray2D
from pScientific.Geometry3.Coordinates3      cimport Coordinates3
from pScientific.Symmetry.SymmetryParameters cimport CSymmetryParameters   , \
                                                     SymmetryParameters

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "PairwiseInteractionMonteCarlo.h":

    ctypedef struct CPairwiseInteractionMonteCarlo "PairwiseInteractionMonteCarlo":
        CInteger  isolateScale
        CReal     buffer
        CReal     chargeScale
        CReal     cutOff
        CReal     epsilonScale
        CReal     sigmaScale
        CReal     underFlowL
        CReal     underFlowQ

    cdef CPairwiseInteractionMonteCarlo *PairwiseInteractionMonteCarlo_Allocate                          ( CStatus                         *status             )
    cdef CPairwiseInteractionMonteCarlo *PairwiseInteractionMonteCarlo_Clone                             ( CPairwiseInteractionMonteCarlo  *self               ,
                                                                                                           CStatus                         *status             )
    cdef void                            PairwiseInteractionMonteCarlo_Deallocate                        ( CPairwiseInteractionMonteCarlo **self               )
    cdef void                            PairwiseInteractionMonteCarlo_Interactions                      ( CPairwiseInteractionMonteCarlo  *self               ,
                                                                                                           CRealArray1D                    *r                  ,
                                                                                                           CRealArray1D                    *electrostatic      ,
                                                                                                           CRealArray1D                    *lennardJonesA      ,
                                                                                                           CRealArray1D                    *lennardJonesB      )
    cdef CReal                           PairwiseInteractionMonteCarlo_MMMMEnergy                        ( CPairwiseInteractionMonteCarlo  *self               ,
                                                                                                           CRealArray1D                    *charges            ,
                                                                                                           CIntegerArray1D                 *ljTypes            ,
                                                                                                           CLJParameterContainer           *ljParameters       ,
                                                                                                           CReal                            electrostaticScale ,
                                                                                                           CReal                            lennardJonesScale  ,
                                                                                                           CSelectionContainer             *isolates           ,
                                                                                                           CBooleanBlock                   *isFree             ,
                                                                                                           CRealArray2D                    *coordinates3       ,
                                                                                                           CSymmetryParameters             *symmetryParameters ,
                                                                                                           CReal                           *eElectrostatic     ,
                                                                                                           CReal                           *eLennardJones      ,
                                                                                                           CStatus                         *status             )
    cdef CReal                           PairwiseInteractionMonteCarlo_MMMMIsolateEnergy                 ( CPairwiseInteractionMonteCarlo  *self               ,
                                                                                                           CInteger                         isolate            ,
                                                                                                           CRealArray1D                    *charges            ,
                                                                                                           CIntegerArray1D                 *ljTypes            ,
                                                                                                           CLJParameterContainer           *ljParameters       ,
                                                                                                           CReal                            electrostaticScale ,
                                                                                                           CReal                            lennardJonesScale  ,
                                                                                                           CSelectionContainer             *isolates           ,
                                                                                                           CBooleanBlock                   *isFree             ,
                                                                                                           CRealArray2D                    *coordinates3       ,
                                                                                                           CSymmetryParameters             *symmetryParameters ,
                                                                                                           CReal                           *eElectrostatic     ,
                                                                                                           CReal                           *eLennardJones      ,
                                                                                                           CStatus                         *status             )
    cdef void                            PairwiseInteractionMonteCarlo_ScaleIsolateInteractionParameters ( CPairwiseInteractionMonteCarlo  *self               ,
                                                                                                           CInteger                         isolate            ,
                                                                                                           CReal                            chargeScale        ,
                                                                                                           CReal                            epsilonScale       ,
                                                                                                           CReal                            sigmaScale         )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class PairwiseInteractionMonteCarlo ( PairwiseInteraction ):

    cdef CPairwiseInteractionMonteCarlo *cObject
    cdef public object                   isOwner
