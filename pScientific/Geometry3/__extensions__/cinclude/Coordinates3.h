# ifndef _COORDINATES3
# define _COORDINATES3

# include "Boolean.h"
# include "Integer.h"
# include "Matrix33.h"
# include "Real.h"
# include "RealArray1D.h"
# include "RealArray2D.h"
# include "RegularGrid.h"
# include "RegularGridOccupancy.h"
# include "Selection.h"
# include "Status.h"
# include "SymmetricMatrix.h"
# include "Transformation3.h"
# include "Vector3.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Type. */
# define Coordinates3                     RealArray2D

/* . Procedures. */
# define Coordinates3_AbsoluteMaximum     RealArray2D_AbsoluteMaximum
# define Coordinates3_Add                 RealArray2D_Add
# define Coordinates3_Clone               RealArray2D_Clone
# define Coordinates3_CopyTo              RealArray2D_CopyTo
# define Coordinates3_Deallocate          RealArray2D_Deallocate
# define Coordinates3_GetItem             RealArray2D_GetItem
# define Coordinates3_Print               RealArray2D_Print
# define Coordinates3_Prune               RealArray2D_PruneByRow
# define Coordinates3_Resize              RealArray2D_Resize
# define Coordinates3_RootMeanSquare      RealArray2D_RootMeanSquare
# define Coordinates3_Scale               RealArray2D_Scale
# define Coordinates3_Set                 RealArray2D_Set
# define Coordinates3_SetItem             RealArray2D_SetItem
# define Coordinates3_View1D              RealArray2D_View1D
# define Coordinates3_View2D              RealArray2D_View
# define Coordinates3_ViewCoordinates3    RealArray2D_View
# define Coordinates3_ViewOfRaw           RealArray2D_ViewOfRaw
# define Coordinates3_ViewVector3         RealArray2D_View1D

/* . Macros. */
# define Coordinates3_Columns             View2D_Columns
# define Coordinates3_Data                Array2D_Data
# define Coordinates3_Extent              View2D_Extent
# define Coordinates3_Item                Array2D_Item
# define Coordinates3_ItemIndex           View2D_ItemIndex
# define Coordinates3_ItemPointer         Array2D_ItemPointer
# define Coordinates3_RowPointer          Array2D_RowPointer
# define Coordinates3_Rows                View2D_Rows
# define Coordinates3_Size                View2D_Size

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define Coordinates3_DecrementRow( self, i, xij, yij, zij ) \
        { \
	    Coordinates3_Item ( self, i, 0 ) -= xij ; \
	    Coordinates3_Item ( self, i, 1 ) -= yij ; \
	    Coordinates3_Item ( self, i, 2 ) -= zij ; \
        }

# define Coordinates3_DifferenceRow( self, i, j, xij, yij, zij ) \
        { \
	    xij = Coordinates3_Item ( self, i, 0 ) - Coordinates3_Item ( self, j, 0 ) ; \
	    yij = Coordinates3_Item ( self, i, 1 ) - Coordinates3_Item ( self, j, 1 ) ; \
	    zij = Coordinates3_Item ( self, i, 2 ) - Coordinates3_Item ( self, j, 2 ) ; \
        }

# define Coordinates3_GetRow( self, i, x, y, z ) \
        { \
	   x = Coordinates3_Item ( self, i, 0 ) ; \
	   y = Coordinates3_Item ( self, i, 1 ) ; \
	   z = Coordinates3_Item ( self, i, 2 ) ; \
        }

# define Coordinates3_IncrementRow( self, i, xij, yij, zij ) \
        { \
	    Coordinates3_Item ( self, i, 0 ) += xij ; \
	    Coordinates3_Item ( self, i, 1 ) += yij ; \
	    Coordinates3_Item ( self, i, 2 ) += zij ; \
        }

# define Coordinates3_ScaleRow( self, i, value ) \
        { \
	    Coordinates3_Item ( self, i, 0 ) *= value ; \
	    Coordinates3_Item ( self, i, 1 ) *= value ; \
	    Coordinates3_Item ( self, i, 2 ) *= value ; \
        }

# define Coordinates3_SetRow( self, i, x, y, z ) \
        { \
	    Coordinates3_Item ( self, i, 0 ) = x ; \
	    Coordinates3_Item ( self, i, 1 ) = y ; \
	    Coordinates3_Item ( self, i, 2 ) = z ; \
        }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Coordinates3 *Coordinates3_Allocate                                ( const Integer                extent            ,
                                                                                  Status                *status            ) ;
extern Real          Coordinates3_Angle                                   ( const Coordinates3          *self              ,
                                                                            const Integer                i                 ,
                                                                            const Integer                j                 ,
                                                                            const Integer                k                 ) ;
extern Status        Coordinates3_BuildPointFromDistance                  (       Coordinates3          *self              ,
                                                                            const Integer                i                 ,
                                                                            const Integer                j                 ,
                                                                            const Real                   r                 ,
                                                                            const Vector3               *direction         ) ;
extern Status        Coordinates3_BuildPointFromDistanceAngle             (       Coordinates3          *self              ,
                                                                            const Integer                i                 ,
                                                                            const Integer                j                 ,
                                                                            const Integer                k                 ,
                                                                            const Real                   r                 ,
                                                                            const Real                   theta             ,
                                                                            const Vector3               *direction         ) ;
extern Status        Coordinates3_BuildPointFromDistanceAngleDihedral     (       Coordinates3          *self              ,
                                                                            const Integer                i                 ,
                                                                            const Integer                j                 ,
                                                                            const Integer                k                 ,
                                                                            const Integer                l                 ,
                                                                            const Real                   r                 ,
                                                                            const Real                   theta             ,
                                                                            const Real                   phi               ) ;
extern Status        Coordinates3_BuildPointFromDistancePlaneAngle        (       Coordinates3          *self              ,
                                                                            const Integer                i                 ,
                                                                            const Integer                j                 ,
                                                                            const Integer                k                 ,
                                                                            const Integer                l                 ,
                                                                            const Real                   r                 ,
                                                                            const Real                   planeangle        ) ;
extern Status        Coordinates3_BuildPointFromDistanceTetrahedralTripod (       Coordinates3          *self              ,
                                                                            const Integer                i                 ,
                                                                            const Integer                j                 ,
                                                                            const Integer                k                 ,
                                                                            const Integer                l                 ,
                                                                            const Integer                m                 ,
                                                                            const Real                   r                 ) ;
extern Status        Coordinates3_Center                                  ( const Coordinates3          *self              ,
                                                                            const Selection             *selection         ,
                                                                            const RealArray1D           *weights           ,
                                                                                  Vector3               **center           ) ;
extern void          Coordinates3_CenterRaw                               ( const Coordinates3          *self              ,
                                                                            const Selection             *selection         ,
                                                                            const RealArray1D           *weights           ,
                                                                                  Real                  *data              ,
                                                                                  Status                *status            ) ;
extern Real          Coordinates3_Dihedral                                ( const Coordinates3          *self              ,
                                                                            const Integer                i                 ,
                                                                            const Integer                j                 ,
                                                                            const Integer                k                 ,
                                                                            const Integer                l                 ) ;
extern Real          Coordinates3_Distance                                ( const Coordinates3          *self              ,
                                                                            const Integer                i                 ,
                                                                            const Integer                j                 ) ;
extern void          Coordinates3_EnclosingOrthorhombicBox                ( const Coordinates3          *self              ,
                                                                            const Selection             *selection         ,
                                                                            const RealArray1D           *radii             ,
                                                                                  Vector3               *origin            ,
                                                                                  Vector3               *extents           ) ;
extern void          Coordinates3_FromRegularGrid                         (       Coordinates3          *self              ,
                                                                            const RegularGrid           *grid              ,
                                                                                  Selection             *selection         ,
                                                                                  Status                *status            ) ;
extern void          Coordinates3_Gather                                  (       Coordinates3          *self              ,
                                                                            const Coordinates3          *other             ,
                                                                            const Selection             *selection         ) ;
extern void          Coordinates3_GatherAdd                               (       Coordinates3          *self              ,
                                                                            const Real                  alpha              ,
                                                                            const Coordinates3          *other             ,
                                                                            const Selection             *selection         ) ;
extern Status        Coordinates3_IdentifyOccupiedGridPoints              (       Coordinates3          *self              ,
                                                                            const RegularGrid           *grid              ,
                                                                            const RealArray1D           *radii             ,
                                                                            const Boolean                QMIDPOINTOVERLAP  ,
                                                                                  Selection            **occupied          ) ;
extern void          Coordinates3_InertiaMatrix                           ( const Coordinates3          *self              ,
                                                                            const Selection             *selection         ,
                                                                            const RealArray1D           *weights           ,
                                                                                  SymmetricMatrix       *inertia           ) ;
extern void          Coordinates3_MakeConformingGrid                      ( const Coordinates3          *self              ,
                                                                            const Selection             *andSelection      ,
                                                                            const RegularGrid           *grid              ,
                                                                                  RegularGrid          **conformingGrid    ,
                                                                                  IntegerArray1D       **offSet            ,
                                                                                  Status                *status            ) ;
extern void          Coordinates3_MakeConformingGridAndOccupancy          ( const Coordinates3          *self              ,
                                                                            const Selection             *andSelection      ,
                                                                            const RegularGrid           *grid              ,
                                                                                  RegularGrid          **conformingGrid    ,
                                                                                  RegularGridOccupancy **occupancy         ,
                                                                                  IntegerArray1D       **offSet            ,
                                                                                  Status                *status            ) ;
extern RegularGrid  *Coordinates3_MakeGrid                                ( const Coordinates3          *self              ,
                                                                            const Selection             *andSelection      ,
                                                                            const Real                   gridSize          ,
                                                                                  Status                *status            ) ;
extern void          Coordinates3_MakeGridAndOccupancy                    ( const Coordinates3          *self              ,
                                                                            const Selection             *andSelection      ,
                                                                            const Real                   gridSize          ,
                                                                                  RegularGrid          **grid              ,
                                                                                  RegularGridOccupancy **occupancy         ,
                                                                                  Status                *status            ) ;
extern void          Coordinates3_MakePeriodicGridAndOccupancy            ( const Coordinates3          *self              ,
                                                                            const Vector3               *boxSize           ,
                                                                            const Real                   gridSize          ,
                                                                                  RegularGrid          **grid              ,
                                                                                  RegularGridOccupancy **occupancy         ,
                                                                                  Status                *status            ) ;
extern void          Coordinates3_MomentsOfInertia                        ( const Coordinates3          *self              ,
                                                                            const Selection             *selection         ,
                                                                            const RealArray1D           *weights           ,
                                                                                  Vector3               *moments           ,
                                                                                  Matrix33              *axes              ) ;
extern Real          Coordinates3_RadiusOfGyration                        ( const Coordinates3          *self              ,
                                                                            const Selection             *selection         ,
                                                                            const RealArray1D           *weights           ) ;
extern Real          Coordinates3_RootMeanSquareDeviation                 ( const Coordinates3          *self              ,
                                                                            const Coordinates3          *other             ,
                                                                            const Selection             *selection         ,
                                                                            const RealArray1D           *weights           ) ;
extern void          Coordinates3_Rotate                                  (       Coordinates3          *self              ,
                                                                            const Matrix33              *rotation          ,
                                                                            const Selection             *selection         ) ;
extern Integer       Coordinates3_RotationTranslationVectors              ( const Coordinates3          *self              ,
                                                                            const RealArray1D           *weights           ,
                                                                            const Boolean                QRx               ,
                                                                            const Boolean                QRy               ,
                                                                            const Boolean                QRz               ,
                                                                            const Boolean                QTx               ,
                                                                            const Boolean                QTy               ,
                                                                            const Boolean                QTz               ,
                                                                                  RealArray2D           *vectors           ,
                                                                                  Status                *status            ) ;
extern void          Coordinates3_ScaleRows                               (       Coordinates3          *self              ,
                                                                            const RealArray1D           *rowScalingFactors ,
                                                                                  Status                *status            ) ;
extern void          Coordinates3_Scatter                                 ( const Coordinates3          *self              ,
                                                                                  Coordinates3          *other             ,
                                                                            const Selection             *selection         ) ;
extern void          Coordinates3_ScatterAdd                              ( const Coordinates3          *self              ,
                                                                            const Real                  alpha              ,
                                                                                  Coordinates3          *other             ,
                                                                            const Selection             *selection         ) ;
extern void          Coordinates3_SetByRow                                (       Coordinates3          *self              ,
                                                                            const Selection             *selection         ,
                                                                            const Real                   alpha             ,
                                                                                  Status                *status            ) ;
extern void          Coordinates3_Superimpose                             (       Coordinates3          *self              ,
                                                                            const Coordinates3          *other             ,
                                                                            const Selection             *selection         ,
                                                                            const RealArray1D           *weights           ,
                                                                                  Matrix33              *rotation          ,
                                                                                  Vector3               *translation       ) ;
extern void          Coordinates3_ToPrincipalAxes                         (       Coordinates3          *self              ,
                                                                            const Selection             *selection         ,
                                                                            const RealArray1D           *weights           ) ;
extern void          Coordinates3_Transform                               (       Coordinates3          *self              ,
                                                                            const Transformation3       *transformation    ,
                                                                            const Selection             *selection         ) ;
extern void          Coordinates3_Translate                               (       Coordinates3          *self              ,
                                                                            const Vector3               *translation       ,
                                                                            const Selection             *selection         ) ;
extern void          Coordinates3_TranslateToCenter                       (       Coordinates3          *self              ,
                                                                            const Selection             *selection         ,
                                                                            const RealArray1D           *weights           ) ;

# endif
