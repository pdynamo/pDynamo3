from pCore.CPrimitiveTypes                           cimport CBoolean                       , \
                                                             CFalse                         , \
                                                             CInteger                       , \
                                                             CReal                          , \
                                                             CTrue
from pCore.SelectionContainer                        cimport CSelectionContainer            , \
                                                             SelectionContainer
from pCore.Status                                    cimport CStatus                        , \
                                                             CStatus_OK
from pMolecule.MMModel.LJParameterContainer          cimport CLJParameterContainer          , \
                                                             LJParameterContainer
from pMolecule.NBModel.PairwiseInteractionMonteCarlo cimport CPairwiseInteractionMonteCarlo , \
                                                             PairwiseInteractionMonteCarlo  , \
                                                             PairwiseInteractionMonteCarlo_MMMMEnergy
from pScientific.Arrays.IntegerArray1D               cimport CIntegerArray1D                , \
                                                             IntegerArray1D
from pScientific.Arrays.RealArray1D                  cimport CRealArray1D                   , \
                                                             RealArray1D
from pScientific.Arrays.RealArray2D                  cimport CRealArray2D
from pScientific.Geometry3.Coordinates3              cimport Coordinates3
from pScientific.Geometry3.Matrix33                  cimport Matrix33
from pScientific.Symmetry.SymmetryParameters         cimport CSymmetryParameters            , \
                                                             SymmetryParameters             , \
                                                             SymmetryParameters_Volume

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "MonteCarloSystemGeometry.h":

    ctypedef struct CMonteCarloSystemGeometry "MonteCarloSystemGeometry":
        CInteger                        blocks
        CInteger                        moves
        CInteger                        nReject
        CInteger                        nRejectM
        CInteger                        nRejectT
        CInteger                        nRejectV
        CInteger                        nTryM
        CInteger                        nTryV
        CReal                           beta
        CReal                           dielectric
        CReal                           eCurrent
        CReal                           pressure
        CReal                           tFactor
        CReal                           volume
        CReal                           acceptanceRatio
        CReal                           rMax
        CReal                           tMax
        CReal                           vMax
        CReal                           eAv
        CReal                           eAv2
        CReal                           eTot
        CReal                           eTot2
        CReal                           eTotB
        CReal                           eTotB2
        CReal                           hAv
        CReal                           hAv2
        CReal                           hTot
        CReal                           hTot2
        CReal                           hTotB
        CReal                           hTotB2
        CReal                           vAv
        CReal                           vAv2
        CReal                           vTot
        CReal                           vTot2
        CReal                           vTotB
        CReal                           vTotB2
        CReal                          *random
        CRealArray2D                   *oldCoordinates3
        CRealArray2D                   *rotation
        CSymmetryParameters            *oldSymmetryParameters
        CRealArray1D                   *translation
        CIntegerArray1D                *ljTypes
        CRealArray2D                   *coordinates3
        CLJParameterContainer          *ljParameters
        CPairwiseInteractionMonteCarlo *pairwiseInteraction
        CRealArray1D                   *charges
        CSelectionContainer            *isolates
        CSymmetryParameters            *symmetryParameters

    cdef CMonteCarloSystemGeometry *MonteCarloSystemGeometry_Allocate                  ( CInteger  numberOfParticles, CInteger  numberOfRandom )
    cdef void                       MonteCarloSystemGeometry_Deallocate                ( CMonteCarloSystemGeometry **self )
    cdef void                       MonteCarloSystemGeometry_AdjustMoveSizes           ( CMonteCarloSystemGeometry  *self )
    cdef CStatus                    MonteCarloSystemGeometry_MoveIsolate               ( CMonteCarloSystemGeometry  *self )
    cdef CStatus                    MonteCarloSystemGeometry_MoveVolume                ( CMonteCarloSystemGeometry  *self )
    cdef void                       MonteCarloSystemGeometry_StatisticsBlockAccumulate ( CMonteCarloSystemGeometry  *self )
    cdef void                       MonteCarloSystemGeometry_StatisticsBlockStart      ( CMonteCarloSystemGeometry  *self )
    cdef void                       MonteCarloSystemGeometry_StatisticsBlockStop       ( CMonteCarloSystemGeometry  *self )
    cdef void                       MonteCarloSystemGeometry_StatisticsStop            ( CMonteCarloSystemGeometry  *self )
    cdef void                       MonteCarloSystemGeometry_StatisticsStart           ( CMonteCarloSystemGeometry  *self )
