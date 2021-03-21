# ifndef _MATRIX33
# define _MATRIX33

# include "Boolean.h"
# include "Integer.h"
# include "Real.h"
# include "RealArray2D.h"
# include "Status.h"
# include "Vector3.h"

/*------------------------------------------------------------------------------
! . Definitions.
!-----------------------------------------------------------------------------*/
/* . Type. */
# define Matrix33 RealArray2D

/* . Procedures. */
# define Matrix33_Add                      RealArray2D_Add
# define Matrix33_CloneDeep                RealArray2D_CloneDeep
# define Matrix33_CloneShallow             RealArray2D_CloneShallow
# define Matrix33_CopyTo                   RealArray2D_CopyTo
# define Matrix33_Deallocate               RealArray2D_Deallocate
# define Matrix33_GetItem                  RealArray2D_GetItem
# define Matrix33_ProjectOutOfVector       RealArray2D_ProjectOutOfVector
# define Matrix33_RootMeanSquare           RealArray2D_RootMeanSquare
# define Matrix33_Scale                    RealArray2D_Scale
# define Matrix33_Set                      RealArray2D_Set
# define Matrix33_SetItem                  RealArray2D_SetItem
# define Matrix33_Transpose                RealArray2D_TransposeSquare

/* . Macros. */
# define Matrix33_Data                     Array2D_Data        
# define Matrix33_Extent                   View2D_Extent
# define Matrix33_Item                     Array2D_Item        
# define Matrix33_ItemIndex                Array2D_ItemIndex   
# define Matrix33_ItemPointer              Array2D_ItemPointer 
# define Matrix33_RowPointer               Array2D_RowPointer 

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define Matrix33_GetRow( self, i, x, y, z ) \
         { \
	   x = Matrix33_Item ( self, i, 0 ) ; \
	   y = Matrix33_Item ( self, i, 1 ) ; \
	   z = Matrix33_Item ( self, i, 2 ) ; \
         }

# define Matrix33_IncrementRow( self, i, a, b, c ) \
         { \
	    Matrix33_Item ( self, i, 0 ) += a ; \
	    Matrix33_Item ( self, i, 1 ) += b ; \
	    Matrix33_Item ( self, i, 2 ) += c ; \
         }

# define Matrix33_SetRow( self, i, x, y, z ) \
         { \
	    Matrix33_Item ( self, i, 0 ) = x ; \
	    Matrix33_Item ( self, i, 1 ) = y ; \
	    Matrix33_Item ( self, i, 2 ) = z ; \
         }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Matrix33 *Matrix33_Allocate               (       void                 ) ;
extern void      Matrix33_ApplyToVector3         ( const Matrix33  *self      ,
                                                         Vector3   *vector3   ) ;
extern void      Matrix33_AngleAxisFromRotation  ( const Matrix33  *self      ,
                                                   const Real      *tolerance ,
                                                         Real      *angle     ,
                                                         Vector3   *axis      ,
                                                         Status    *status    ) ;
extern Real      Matrix33_Determinant            ( const Matrix33  *self      ) ;
extern void      Matrix33_InverseDerivative      ( const Matrix33  *self      ,
                                                   const Integer    i         ,
                                                   const Integer    j         ,
                                                         Matrix33  *other     ) ;
extern void      Matrix33_Invert                 (       Matrix33  *self      ,
                                                   const Matrix33  *other     ) ;
extern Boolean   Matrix33_IsEqual                ( const Matrix33  *self      ,
                                                   const Matrix33  *other     ) ;
extern Boolean   Matrix33_IsIdentity             ( const Matrix33  *self      ) ;
extern Boolean   Matrix33_IsImproperRotation     ( const Matrix33  *self      ) ;
extern Boolean   Matrix33_IsOrthogonal           ( const Matrix33  *self      ) ;
extern Boolean   Matrix33_IsProperRotation       ( const Matrix33  *self      ) ;
extern void      Matrix33_PostMultiplyBy         (       Matrix33  *self      ,
                                                   const Matrix33  *other     ) ;
extern void      Matrix33_PreMultiplyBy          (       Matrix33  *self      ,
                                                   const Matrix33  *other     ) ;
extern Status    Matrix33_Reflection             (       Matrix33 **self      ,
                                                   const Vector3   *normal    ) ;
extern Status    Matrix33_RotationAboutAxis      (       Matrix33 **self      ,
                                                   const Real       angle     ,
                                                   const Real       x         ,
                                                   const Real       y         ,
                                                   const Real       z         ) ;
extern Status    Matrix33_RotationFromQuaternion (       Matrix33 **self      ,
                                                   const Real       q0        ,
                                                   const Real       q1        ,
                                                   const Real       q2        ,
                                                   const Real       q3        ) ;

# endif
