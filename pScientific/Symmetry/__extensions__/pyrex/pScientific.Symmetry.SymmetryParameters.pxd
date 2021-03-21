from pCore.BooleanBlock                 cimport CBooleanBlock             , \
                                                BooleanBlock
from pCore.CPrimitiveTypes              cimport CBoolean                  , \
                                                CFalse                    , \
                                                CInteger                  , \
                                                CReal                     , \
                                                CTrue
from pCore.Selection                    cimport Selection                 , \
                                                CSelection
from pCore.SelectionContainer           cimport SelectionContainer        , \
                                                CSelectionContainer
from pCore.Status                       cimport CStatus                   , \
                                                CStatus_OK
from pScientific.Arrays.RealArray1D     cimport CRealArray1D              , \
                                                RealArray1D_PointerToData
from pScientific.Arrays.RealArray2D     cimport CRealArray2D
from pScientific.Geometry3.Coordinates3 cimport Coordinates3
from pScientific.Geometry3.Matrix33     cimport Matrix33
from pScientific.Geometry3.Vector3      cimport Vector3

#===================================================================================================================================
# . Declarations.
#===================================================================================================================================
cdef extern from "SymmetryParameters.h":

    ctypedef struct CSymmetryParameters "SymmetryParameters":
        CReal         a
        CReal         b
        CReal         c
        CReal         alpha
        CReal         beta
        CReal         gamma
        CRealArray2D *H
        CRealArray2D *inverseH

    cdef CSymmetryParameters *SymmetryParameters_AllocateWithMatrices              ( CRealArray2D         *H            ,
                                                                                     CRealArray2D         *inverseH     ,
                                                                                     CStatus              *status       )
    cdef void                 SymmetryParameters_CenterCoordinates3ByFreeIsolate   ( CSymmetryParameters  *self         ,
                                                                                     CSelectionContainer  *isolates     ,
                                                                                     CBooleanBlock        *freeIsolates ,
                                                                                     CRealArray2D         *coordinates3 ,
                                                                                     CStatus              *status       )
    cdef void                 SymmetryParameters_CenterCoordinates3ByIndex         ( CSymmetryParameters  *self         ,
                                                                                     CSelection           *selection    ,
                                                                                     CRealArray2D         *coordinates3 ,
                                                                                     CStatus              *status       )
    cdef void                 SymmetryParameters_CenterCoordinates3ByIsolate       ( CSymmetryParameters  *self         ,
                                                                                     CSelectionContainer  *isolates     ,
                                                                                     CSelection           *selection    ,
                                                                                     CRealArray2D         *coordinates3 ,
                                                                                     CStatus              *status       )
    cdef void                 SymmetryParameters_CopyTo                            ( CSymmetryParameters  *self         ,
                                                                                     CSymmetryParameters  *other        )
    cdef void                 SymmetryParameters_Deallocate                        ( CSymmetryParameters **self         )
    cdef CBoolean             SymmetryParameters_IsMinimumImageConventionSatisfied ( CSymmetryParameters  *self         ,
                                                                                     CReal                 length       )
    cdef void                 SymmetryParameters_MakeMinimumImageVector            ( CSymmetryParameters  *self         ,
                                                                                     CReal                *r            ,
                                                                                     CReal                *dr           )
    cdef void                 SymmetryParameters_SetCrystalParameters              ( CSymmetryParameters  *self         ,
                                                                                     CReal                 a            ,
                                                                                     CReal                 b            ,
                                                                                     CReal                 c            ,
                                                                                     CReal                 alpha        ,
                                                                                     CReal                 beta         ,
                                                                                     CReal                 gamma        )
    cdef CReal                SymmetryParameters_Volume                            ( CSymmetryParameters  *self         )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class SymmetryParameters:

    cdef CSymmetryParameters *cObject
    cdef public object        isOwner
    cdef public Matrix33      H
    cdef public Matrix33      inverseH
