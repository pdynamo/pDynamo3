# ifndef _SYMMETRICMATRIX
# define _SYMMETRICMATRIX

# include "Boolean.h"
# include "Integer.h"
# include "IntegerBlock.h"
# include "Iterator.h"
# include "Real.h"
# include "RealArray1D.h"
# include "RealArray2D.h"
# include "RealBlock.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The symmetric matrix type. */
typedef struct {
    Integer    extent;
    Integer    size  ;
    RealBlock *block ;
    Real      *data  ;
} SymmetricMatrix ;

/* . Updating options. */
typedef enum {
    SymmetricMatrixUpdating_BFGS   = 0,
    SymmetricMatrixUpdating_Bofill = 1,
    SymmetricMatrixUpdating_MS     = 2,
    SymmetricMatrixUpdating_Powell = 3
} SymmetricMatrixUpdating_Option ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . A pointer to the start of the data. */
# define SymmetricMatrix_Data( self ) ( &((self)->data[0]) )

/* . The extent. */
# define SymmetricMatrix_Extent( self ) ( (self) == NULL ? 0 : (self)->extent )

/* . The index to an item (i >= j). */
# define SymmetricMatrix_ItemIndex( i, j ) ( ( (i) * ( i + 1 ) ) / 2 + j )

/* . An item. */
# define SymmetricMatrix_Item( self, i, j ) ( (self)->data[ SymmetricMatrix_ItemIndex ( i, j ) ] )

/* . The size. */
# define SymmetricMatrix_Size( self ) ( (self) == NULL ? 0 : (self)->size )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Real             SymmetricMatrix_AbsoluteMaximum          ( const SymmetricMatrix                *self            ) ;
extern void             SymmetricMatrix_Add                      (       SymmetricMatrix                *self            ,
                                                                   const Real                            alpha           ,
                                                                   const SymmetricMatrix                *other           ,
                                                                         Status                         *status          ) ;
extern SymmetricMatrix *SymmetricMatrix_Allocate                 (       Status                         *status          ) ;
extern SymmetricMatrix *SymmetricMatrix_AllocateWithExtent       ( const Integer                         extent          ,
                                                                         Status                         *status          ) ;
extern void             SymmetricMatrix_AnticommutatorSS         (       SymmetricMatrix                *self            ,
                                                                   const SymmetricMatrix                *a               ,
                                                                   const SymmetricMatrix                *b               ,
                                                                         Status                         *status          ) ;
extern SymmetricMatrix *SymmetricMatrix_CloneDeep                ( const SymmetricMatrix                *self            ,
                                                                         Status                         *status          ) ;
extern SymmetricMatrix *SymmetricMatrix_CloneShallow             ( const SymmetricMatrix                *self            ,
                                                                         Status                         *status          ) ;
extern void             SymmetricMatrix_CopyFromRealArray2D      (       SymmetricMatrix                *self            ,
                                                                   const RealArray2D                    *other           ,
                                                                         Status                         *status          ) ;
extern void             SymmetricMatrix_CopyTo                   ( const SymmetricMatrix                *self            ,
                                                                         SymmetricMatrix                *other           ,
                                                                         Status                         *status          ) ;
extern void             SymmetricMatrix_CopyToRealArray2D        ( const SymmetricMatrix                *self            ,
                                                                         RealArray2D                    *other           ,
                                                                         Status                         *status          ) ;
extern void             SymmetricMatrix_Deallocate               (       SymmetricMatrix               **self            ) ;
extern void             SymmetricMatrix_DiagonalOfProduct        ( const SymmetricMatrix                *self            ,
                                                                   const SymmetricMatrix                *other           ,
                                                                         RealArray1D                    *result          ,
                                                                         Status                         *status          ) ; 
extern void             SymmetricMatrix_DiagonalOfTransform      ( const SymmetricMatrix                *self            ,
                                                                   const RealArray2D                    *matrix          ,
                                                                   const Boolean                         useTranspose    ,
                                                                         RealArray1D                    *result          ,
                                                                         Status                         *status          ) ;
extern void             SymmetricMatrix_GetColumn                ( const SymmetricMatrix                *self            ,
                                                                   const Integer                         n               ,
                                                                         RealArray1D                    *column          ,
                                                                         Status                         *status          ) ;
extern SymmetricMatrix *SymmetricMatrix_FromExtentBlock          ( const Integer                         extent          ,
                                                                         RealBlock                      *block           ,
                                                                   const Boolean                         withReference   ,
                                                                         Status                         *status          ) ;
extern Real             SymmetricMatrix_GetItem                  ( const SymmetricMatrix                *self            ,
                                                                   const Integer                         i               ,
                                                                   const Integer                         j               ,
                                                                         Status                         *status          ) ;
extern void             SymmetricMatrix_Increment                (       SymmetricMatrix                *self            ,
                                                                   const Real                            value           ) ;
extern void             SymmetricMatrix_IndexedCopyToRealArray2D ( const SymmetricMatrix                *self            ,
                                                                   const IntegerBlock                   *indices         ,
                                                                         RealArray2D                    *target          ,
                                                                         Status                         *status          ) ;
extern Boolean          SymmetricMatrix_IsDiagonal               ( const SymmetricMatrix                *self            ,
                                                                   const Real                            tolerance       ) ;
extern void             SymmetricMatrix_MakeFromEigensystem      (       SymmetricMatrix                *self            ,
                                                                   const Boolean                         zeroMatrix      ,
                                                                   const Integer                         numberOfVectors ,
                                                                   const RealArray1D                    *eigenValues     ,
                                                                   const RealArray2D                    *eigenVectors    ,
                                                                         Status                         *status          ) ;
extern Iterator        *SymmetricMatrix_MakeIterator             ( const SymmetricMatrix                *self            ,
                                                                         Status                         *status          ) ;
extern void             SymmetricMatrix_PostMatrixMultiply       ( const SymmetricMatrix                *self            ,
                                                                   const RealArray2D                    *matrix          ,
                                                                   const Boolean                         useTranspose    ,
                                                                         RealArray2D                    *result          ,
                                                                         Status                         *status          ) ;
extern void             SymmetricMatrix_PreMatrixMultiply        ( const SymmetricMatrix                *self            ,
                                                                   const RealArray2D                    *matrix          ,
                                                                   const Boolean                         useTranspose    ,
                                                                         RealArray2D                    *result          ,
                                                                         Status                         *status          ) ;
extern void             SymmetricMatrix_Print                    ( const SymmetricMatrix                *self            ) ;
extern void             SymmetricMatrix_ProjectionMatrix         (       SymmetricMatrix                *self            ,
                                                                   const RealArray2D                    *vectors         ,
                                                                         Status                         *status          ) ;
extern void             SymmetricMatrix_ProjectOut               (       SymmetricMatrix                *self            ,
                                                                   const RealArray2D                    *vectors         ,
                                                                         SymmetricMatrix                *result          ,
                                                                         Status                         *status          ) ;
extern void             SymmetricMatrix_Rank1Update              (       SymmetricMatrix                *self            ,
                                                                   const Real                            alpha           ,
                                                                   const RealArray1D                    *vector          ,
                                                                         Status                         *status          ) ;
extern void             SymmetricMatrix_Raise                    (       SymmetricMatrix                *self            ,
                                                                   const RealArray2D                    *vectors         ,
                                                                   const Real                            value           ,
                                                                         Status                         *status          ) ;
extern Real             SymmetricMatrix_RootMeanSquare           ( const SymmetricMatrix                *self            ) ;
extern void             SymmetricMatrix_Scale                    (       SymmetricMatrix                *self            ,
                                                                   const Real                            value           ) ;
extern void             SymmetricMatrix_ScaleDiagonal            (       SymmetricMatrix                *self            ,
                                                                   const Real                            value           ) ;
extern void             SymmetricMatrix_ScaleOffDiagonal         (       SymmetricMatrix                *self            ,
                                                                   const Real                            value           ) ;
extern void             SymmetricMatrix_Set                      (       SymmetricMatrix                *self            ,
                                                                   const Real                            value           ) ;
extern void             SymmetricMatrix_SetColumn                (       SymmetricMatrix                *self            ,
                                                                   const Integer                         column          ,
                                                                   const RealArray1D                    *vector          ,
                                                                         Status                         *status          ) ;
extern void             SymmetricMatrix_SetItem                  (       SymmetricMatrix                *self            ,
                                                                   const Integer                         i               ,
                                                                   const Integer                         j               ,
                                                                   const Real                            value           ,
                                                                         Status                         *status          ) ;
extern Real             SymmetricMatrix_Sparsity                 ( const SymmetricMatrix                *self            ,
                                                                   const Real                            tolerance       ) ;
extern void             SymmetricMatrix_SumDifference            (       SymmetricMatrix                *self            ,
                                                                         SymmetricMatrix                *other           ,
                                                                         Status                         *status          ) ;
extern void             SymmetricMatrix_SymmetricMatrixMultiply  ( const SymmetricMatrix                *self            ,
                                                                   const SymmetricMatrix                *other           ,
                                                                         RealArray2D                    *result          ,
                                                                         Status                         *status          ) ;
extern void             SymmetricMatrix_SymmetricTransform       ( const SymmetricMatrix                *self            ,
                                                                   const SymmetricMatrix                *other           ,
                                                                         SymmetricMatrix                *result          ,
                                                                         Status                         *status          ) ;
extern Real             SymmetricMatrix_Trace                    ( const SymmetricMatrix                *self            ) ;
extern Real             SymmetricMatrix_TraceOfProduct           ( const SymmetricMatrix                *self            ,
                                                                   const SymmetricMatrix                *other           ,
                                                                         Status                         *status          ) ;
extern void             SymmetricMatrix_Transform                ( const SymmetricMatrix                *self            ,
                                                                   const RealArray2D                    *matrix          ,
                                                                   const Boolean                         useTranspose    ,
                                                                         SymmetricMatrix                *result          ,
                                                                         Status                         *status          ) ;
extern void             SymmetricMatrix_Update                   (       SymmetricMatrix                *self            ,
                                                                   const RealArray1D                    *dx              ,
                                                                   const RealArray1D                    *dg              ,
                                                                   const SymmetricMatrixUpdating_Option  option          ,
                                                                   const Real                           *tolerance       ,
                                                                         Status                         *status          ) ;
extern void             SymmetricMatrix_VectorMultiply           ( const SymmetricMatrix                *self            ,
                                                                   const RealArray1D                    *other           ,
                                                                         RealArray1D                    *result          ,
                                                                         Status                         *status          ) ;
extern Integer          SymmetricMatrix_ViewSize                 ( const Integer                         extent          ) ;

# endif
