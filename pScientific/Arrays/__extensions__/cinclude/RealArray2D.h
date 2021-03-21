# ifndef _REALARRAY2D
# define _REALARRAY2D

# include "Boolean.h"
# include "Integer.h"
# include "Real.h"
# include "RealArray1D.h"
# include "RealBlock.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Basic definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _ArrayDataType Real
# include "Array2D_Header.i"
# undef  _ArrayDataType

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void         RealArray2D_DiagonalOfProduct        ( const RealArray2D *self              ,
                                                           const Boolean      sTranspose        ,
                                                           const RealArray2D *other             ,
                                                           const Boolean      oTranspose        ,
                                                                 RealArray1D *diagonal          ,
                                                                 Status      *status            ) ;
extern Integer      RealArray2D_GramSchmidtOrthogonalize ( const RealArray2D *self              ,
                                                           const Integer     *maximumIterations ,
                                                           const Integer     *numberConstant    ,
                                                           const Real        *tolerance         ,
                                                                 Status      *status            ) ;
extern Boolean      RealArray2D_IsOrthogonal             ( const RealArray2D *self              ,
                                                           const Real        *tolerance         ,
                                                                 Real        *deviation         ,
                                                                 Status      *status            ) ;
extern Boolean      RealArray2D_IsSymmetric              ( const RealArray2D *self              ,
                                                           const Real        *tolerance         ,
                                                                 Real        *deviation         ) ;
extern void         RealArray2D_MatrixMultiply           ( const Boolean      aTranspose        ,
                                                           const Boolean      bTranspose        ,
                                                           const Real         alpha             ,
                                                           const RealArray2D *a                 ,
                                                           const RealArray2D *b                 ,
                                                           const Real         beta              ,
                                                                 RealArray2D *c                 ,
                                                                 Status      *status            ) ;
extern void         RealArray2D_ProjectOutOfArray1D      ( const RealArray2D *self              ,
                                                                 RealArray1D *vector            ,
                                                                 Status      *status            ) ;
extern Real         RealArray2D_Trace                    ( const RealArray2D *self              ,
                                                                 Status      *status            ) ;
extern Real         RealArray2D_TraceOfProduct           ( const RealArray2D *self              ,
                                                           const RealArray2D *other             ,
                                                                 Status      *status            ) ;
extern RealArray2D *RealArray2D_TransposeClone           ( const RealArray2D *self              ,
                                                                 Status      *status            ) ;
extern void         RealArray2D_TransposeGeneral         (       RealArray2D *self              ,
                                                                 Status      *status            ) ;
extern void         RealArray2D_TransposeSquare          (       RealArray2D *self              ,
                                                                 Status      *status            ) ;
extern void         RealArray2D_VectorMultiply           ( const Boolean      aTranspose        ,
                                                           const Real         alpha             ,
                                                           const RealArray2D *a                 ,
                                                           const RealArray1D *x                 ,
                                                           const Real         beta              ,
                                                                 RealArray1D *y                 ,
                                                                 Status      *status            ) ;

# endif
