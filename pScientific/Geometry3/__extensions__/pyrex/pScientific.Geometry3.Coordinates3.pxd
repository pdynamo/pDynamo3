from pCore.CPrimitiveTypes                      cimport CBoolean                   , \
                                                        CFalse                     , \
                                                        CInteger                   , \
                                                        CReal                      , \
                                                        CTrue
from pCore.Selection                            cimport CSelection                 , \
                                                        Selection
from pCore.Status                               cimport CStatus                    , \
                                                        CStatus_OK
from pScientific.Arrays.IntegerArray1D          cimport CIntegerArray1D            , \
                                                        IntegerArray1D
from pScientific.Arrays.RealArray1D             cimport CRealArray1D               , \
                                                        RealArray1D
from pScientific.Arrays.RealArray2D             cimport CRealArray2D               , \
                                                        RealArray2D                , \
                                                        RealArray2D_VectorMultiply
from pScientific.Arrays.Slicing                 cimport MultiSlice_Allocate        , \
                                                        MultiSlice_GetExtent
from pScientific.Arrays.SymmetricMatrix         cimport CSymmetricMatrix           , \
                                                        SymmetricMatrix
from pScientific.Geometry3.Matrix33             cimport Matrix33
from pScientific.Geometry3.RegularGrid          cimport CRegularGrid               , \
                                                        RegularGrid
from pScientific.Geometry3.RegularGridOccupancy cimport CRegularGridOccupancy      , \
                                                        RegularGridOccupancy
from pScientific.Geometry3.Transformation3      cimport CTransformation3           , \
                                                        Transformation3
from pScientific.Geometry3.Vector3              cimport Vector3

#===============================================================================
# . Declarations.
#===============================================================================
cdef extern from "Coordinates3.h":

# . As subclassing need to keep CCoordinates3 as CRealArray2D.
#    ctypedef struct CRealArray2D  "Coordinates3":
#        pass

    cdef void           Coordinates3_Add                                     ( CRealArray2D           *self              ,
                                                                               CReal                   alpha             ,
                                                                               CRealArray2D           *other             ,
                                                                               CStatus                *status            )
    cdef CRealArray2D  *Coordinates3_Allocate                                ( CInteger                extent            ,
                                                                               CStatus                *status            )
    cdef CReal          Coordinates3_Angle                                   ( CRealArray2D           *self              ,
                                                                               CInteger                i                 ,
                                                                               CInteger                j                 ,
                                                                               CInteger                k                 )
    cdef CStatus        Coordinates3_BuildPointFromDistance                  ( CRealArray2D           *self              ,
                                                                               CInteger                i                 ,
                                                                               CInteger                j                 ,
                                                                               CReal                   r                 ,
                                                                               CRealArray1D           *direction         )
    cdef CStatus        Coordinates3_BuildPointFromDistanceAngle             ( CRealArray2D           *self              ,
                                                                               CInteger                i                 ,
                                                                               CInteger                j                 ,
                                                                               CInteger                k                 ,
                                                                               CReal                   r                 ,
                                                                               CReal                   theta             ,
                                                                               CRealArray1D           *direction         )
    cdef CStatus        Coordinates3_BuildPointFromDistanceAngleDihedral     ( CRealArray2D           *self              ,
                                                                               CInteger                i                 ,
                                                                               CInteger                j                 ,
                                                                               CInteger                k                 ,
                                                                               CInteger                l                 ,
                                                                               CReal                   r                 ,
                                                                               CReal                   theta             ,
                                                                               CReal                   phi               )
    cdef CStatus        Coordinates3_BuildPointFromDistancePlaneAngle        ( CRealArray2D           *self              ,
                                                                               CInteger                i                 ,
                                                                               CInteger                j                 ,
                                                                               CInteger                k                 ,
                                                                               CInteger                l                 ,
                                                                               CReal                   r                 ,
                                                                               CReal                   planeangle        )
    cdef CStatus        Coordinates3_BuildPointFromDistanceTetrahedralTripod ( CRealArray2D           *self              ,
                                                                               CInteger                i                 ,
                                                                               CInteger                j                 ,
                                                                               CInteger                k                 ,
                                                                               CInteger                l                 ,
                                                                               CInteger                m                 ,
                                                                               CReal                   r                 )
    cdef CStatus        Coordinates3_Center                                  ( CRealArray2D           *self              ,
                                                                               CSelection             *selection         ,
                                                                               CRealArray1D           *weights           ,
                                                                               CRealArray1D           **center           )
    cdef void           Coordinates3_CenterRaw                               ( CRealArray2D           *self              ,
                                                                               CSelection             *selection         ,
                                                                               CRealArray1D           *weights           ,
                                                                               CReal                  *data              ,
                                                                               CStatus                *status            )
    cdef void           Coordinates3_CopyTo                                  ( CRealArray2D           *self              ,
                                                                               CRealArray2D           *other             ,
                                                                               CStatus                *status            )
    cdef CReal          Coordinates3_Dihedral                                ( CRealArray2D           *self              ,
                                                                               CInteger                i                 ,
                                                                               CInteger                j                 ,
                                                                               CInteger                k                 ,
                                                                               CInteger                l                 )
    cdef CReal          Coordinates3_Distance                                ( CRealArray2D           *self              ,
                                                                               CInteger                i                 ,
                                                                               CInteger                j                 )
    cdef void           Coordinates3_EnclosingOrthorhombicBox                ( CRealArray2D           *self              ,
                                                                               CSelection             *selection         ,
                                                                               CRealArray1D           *radii             ,
                                                                               CRealArray1D           *origin            ,
                                                                               CRealArray1D           *extents           )
    cdef void           Coordinates3_FromRegularGrid                         ( CRealArray2D           *self              ,
                                                                               CRegularGrid           *grid              ,
                                                                               CSelection             *selection         ,
                                                                               CStatus                *status            )
    cdef void           Coordinates3_Gather                                  ( CRealArray2D           *self              ,
                                                                               CRealArray2D           *other             ,
                                                                               CSelection             *selection         )
    cdef void           Coordinates3_GatherAdd                               ( CRealArray2D           *self              ,
                                                                               CReal                  alpha              ,
                                                                               CRealArray2D           *other             ,
                                                                               CSelection             *selection         )
    cdef CStatus        Coordinates3_IdentifyOccupiedGridPoints              ( CRealArray2D           *self              ,
                                                                               CRegularGrid           *grid              ,
                                                                               CRealArray1D           *radii             ,
                                                                               CBoolean                QMIDPOINTOVERLAP  ,
                                                                               CSelection            **occupied          )
    cdef void           Coordinates3_InertiaMatrix                           ( CRealArray2D           *self              ,
                                                                               CSelection             *selection         ,
                                                                               CRealArray1D           *weights           ,
                                                                               CSymmetricMatrix       *inertia           )
    cdef void           Coordinates3_MakeConformingGrid                      ( CRealArray2D           *self              ,
                                                                               CSelection             *andSelection      ,
                                                                               CRegularGrid           *grid              ,
                                                                               CRegularGrid          **conformingGrid    ,
                                                                               CIntegerArray1D       **offSet            ,
                                                                               CStatus                *status            )
    cdef void           Coordinates3_MakeConformingGridAndOccupancy          ( CRealArray2D           *self              ,
                                                                               CSelection             *andSelection      ,
                                                                               CRegularGrid           *grid              ,
                                                                               CRegularGrid          **conformingGrid    ,
                                                                               CRegularGridOccupancy **occupancy         ,
                                                                               CIntegerArray1D       **offSet            ,
                                                                               CStatus                *status            )
    cdef CRegularGrid  *Coordinates3_MakeGrid                                ( CRealArray2D           *self              ,
                                                                               CSelection             *andSelection      ,
                                                                               CReal                   gridSize          ,
                                                                               CStatus                *status            )
    cdef void           Coordinates3_MakeGridAndOccupancy                    ( CRealArray2D           *self              ,
                                                                               CSelection             *andSelection      ,
                                                                               CReal                   gridSize          ,
                                                                               CRegularGrid          **grid              ,
                                                                               CRegularGridOccupancy **occupancy         ,
                                                                               CStatus                *status            )
    cdef void           Coordinates3_MakePeriodicGridAndOccupancy            ( CRealArray2D           *self              ,
                                                                               CRealArray1D           *boxSize           ,
                                                                               CReal                   gridSize          ,
                                                                               CRegularGrid          **grid              ,
                                                                               CRegularGridOccupancy **occupancy         ,
                                                                               CStatus                *status            )
    cdef void           Coordinates3_MomentsOfInertia                        ( CRealArray2D           *self              ,
                                                                               CSelection             *selection         ,
                                                                               CRealArray1D           *weights           ,
                                                                               CRealArray1D           *moments           ,
                                                                               CRealArray2D           *axes              )
    cdef CReal          Coordinates3_RadiusOfGyration                        ( CRealArray2D           *self              ,
                                                                               CSelection             *selection         ,
                                                                               CRealArray1D           *weights           )
    cdef CReal          Coordinates3_RootMeanSquareDeviation                 ( CRealArray2D           *self              ,
                                                                               CRealArray2D           *other             ,
                                                                               CSelection             *selection         ,
                                                                               CRealArray1D           *weights           )
    cdef void           Coordinates3_Rotate                                  ( CRealArray2D           *self              ,
                                                                               CRealArray2D           *rotation          ,
                                                                               CSelection             *selection         )
    cdef CInteger       Coordinates3_RotationTranslationVectors              ( CRealArray2D           *self              ,
                                                                               CRealArray1D           *weights           ,
                                                                               CBoolean                QRx               ,
                                                                               CBoolean                QRy               ,
                                                                               CBoolean                QRz               ,
                                                                               CBoolean                QTx               ,
                                                                               CBoolean                QTy               ,
                                                                               CBoolean                QTz               ,
                                                                               CRealArray2D           *vectors           ,
                                                                               CStatus                *status            )
    cdef void           Coordinates3_ScaleRows                               ( CRealArray2D           *self              ,
                                                                               CRealArray1D           *rowScalingFactors ,
                                                                               CStatus                *status            )
    cdef void           Coordinates3_Scatter                                 ( CRealArray2D           *self              ,
                                                                               CRealArray2D           *other             ,
                                                                               CSelection             *selection         )
    cdef void           Coordinates3_ScatterAdd                              ( CRealArray2D           *self              ,
                                                                               CReal                  alpha              ,
                                                                               CRealArray2D           *other             ,
                                                                               CSelection             *selection         )
    cdef void           Coordinates3_SetByRow                                ( CRealArray2D           *self              ,
                                                                               CSelection             *selection         ,
                                                                               CReal                   value             ,
                                                                               CStatus                *status            )
    cdef void           Coordinates3_Superimpose                             ( CRealArray2D           *self              ,
                                                                               CRealArray2D           *other             ,
                                                                               CSelection             *selection         ,
                                                                               CRealArray1D           *weights           ,
                                                                               CRealArray2D           *rotation          ,
                                                                               CRealArray1D           *translation       )
    cdef void           Coordinates3_ToPrincipalAxes                         ( CRealArray2D           *self              ,
                                                                               CSelection             *selection         ,
                                                                               CRealArray1D           *weights           )
    cdef void           Coordinates3_Transform                               ( CRealArray2D           *self              ,
                                                                               CTransformation3       *transformation    ,
                                                                               CSelection             *selection         )
    cdef void           Coordinates3_Translate                               ( CRealArray2D           *self              ,
                                                                               CRealArray1D           *translation       ,
                                                                               CSelection             *selection         )
    cdef void           Coordinates3_TranslateToCenter                       ( CRealArray2D           *self              ,
                                                                               CSelection             *selection         ,
                                                                               CRealArray1D           *weights           )

#===============================================================================
# . Class.
#===============================================================================
cdef class Coordinates3 ( RealArray2D ):

    cdef public object undefined

