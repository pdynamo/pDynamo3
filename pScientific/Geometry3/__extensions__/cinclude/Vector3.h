# ifndef _VECTOR3
# define _VECTOR3

# include "Boolean.h"
# include "RealArray1D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Type. */
# define Vector3 RealArray1D

/* . Procedures. */
# define Vector3_AbsoluteMaximum      RealArray1D_AbsoluteMaximum
# define Vector3_AbsoluteMaximumIndex RealArray1D_AbsoluteMaximumIndex
# define Vector3_Add                  RealArray1D_Add
# define Vector3_CloneDeep            RealArray1D_CloneDeep
# define Vector3_CloneShallow         RealArray1D_CloneShallow
# define Vector3_CopyTo               RealArray1D_CopyTo
# define Vector3_Deallocate           RealArray1D_Deallocate
# define Vector3_Divide               RealArray1D_Divide
# define Vector3_Dot                  RealArray1D_Dot
# define Vector3_GetItem              RealArray1D_GetItem
# define Vector3_Increment            RealArray1D_Increment
# define Vector3_Maximum              RealArray1D_Maximum
# define Vector3_Minimum              RealArray1D_Minimum
# define Vector3_Multiply             RealArray1D_Multiply
# define Vector3_Norm2                RealArray1D_Norm2
# define Vector3_Normalize            RealArray1D_Normalize
# define Vector3_RootMeanSquare       RealArray1D_RootMeanSquare
# define Vector3_Scale                RealArray1D_Scale
# define Vector3_Set                  RealArray1D_Set
# define Vector3_SetItem              RealArray1D_SetItem
# define Vector3_Sum                  RealArray1D_Sum

/* . Macros. */
# define Vector3_Data                 Array1D_Data
# define Vector3_Extent               View1D_Extent
# define Vector3_Item                 Array1D_Item
# define Vector3_ItemPointer          Array1D_ItemPointer

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Vector3 *Vector3_Allocate     (       void           ) ;
extern void     Vector3_CrossProduct (       Vector3 *self  ,
                                       const Vector3 *other ) ;
extern Boolean  Vector3_IsEqual      ( const Vector3 *self  ,
                                       const Vector3 *other ) ;
extern Boolean  Vector3_IsNull       ( const Vector3 *self  ) ;

# endif
