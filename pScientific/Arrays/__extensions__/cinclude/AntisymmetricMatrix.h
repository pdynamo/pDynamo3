# ifndef _ANTISYMMETRICMATRIX
# define _ANTISYMMETRICMATRIX

# include "Boolean.h"
# include "Integer.h"
# include "Iterator.h"
# include "Real.h"
# include "RealArray1D.h"
# include "RealArray2D.h"
# include "RealBlock.h"
# include "Status.h"
# include "SymmetricMatrix.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The array type. */
typedef struct {
    Integer    extent;
    Integer    size  ;
    RealBlock *block ;
    Real      *data  ;
} AntisymmetricMatrix ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . A pointer to the start of the data. */
# define AntisymmetricMatrix_Data( self ) ( &((self)->data[0]) )

/* . The extent. */
# define AntisymmetricMatrix_Extent( self ) ( (self) == NULL ? 0 : (self)->extent )

/* . The index to an item (i >= j). */
# define AntisymmetricMatrix_ItemIndex( i, j ) ( ( (i) * ( i - 1 ) ) / 2 + j )

/* . An item. */
# define AntisymmetricMatrix_Item( self, i, j ) ( (self)->data[ AntisymmetricMatrix_ItemIndex ( i, j ) ] )

/* . The size. */
# define AntisymmetricMatrix_Size( self ) ( (self) == NULL ? 0 : (self)->size )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Real                 AntisymmetricMatrix_AbsoluteMaximum        ( const AntisymmetricMatrix  *self          ) ;
extern void                 AntisymmetricMatrix_Add                    (       AntisymmetricMatrix  *self          ,
                                                                         const Real                  value         ,
                                                                         const AntisymmetricMatrix  *other         ,
                                                                               Status               *status        ) ;
extern AntisymmetricMatrix *AntisymmetricMatrix_Allocate               (       Status               *status        ) ;
extern AntisymmetricMatrix *AntisymmetricMatrix_AllocateWithExtent     ( const Integer               extent        ,
                                                                               Status               *status        ) ;
extern void                 AntisymmetricMatrix_AnticommutatorAS       ( const AntisymmetricMatrix  *self          ,
                                                                         const SymmetricMatrix      *a             ,
                                                                               AntisymmetricMatrix  *result        ,
                                                                               Status               *status        ) ;
extern AntisymmetricMatrix *AntisymmetricMatrix_CloneDeep              ( const AntisymmetricMatrix  *self          ,  
                                                                               Status               *status        ) ;
extern AntisymmetricMatrix *AntisymmetricMatrix_CloneShallow           ( const AntisymmetricMatrix  *self          ,  
                                                                               Status               *status        ) ;
extern void                 AntisymmetricMatrix_CommutatorAS           ( const AntisymmetricMatrix  *self          ,
                                                                         const SymmetricMatrix      *a             ,
                                                                               SymmetricMatrix      *result        ,
                                                                               Status               *status        ) ;
extern void                 AntisymmetricMatrix_CommutatorSS_Fast      (       AntisymmetricMatrix  *self          ,
                                                                         const SymmetricMatrix      *a             ,
                                                                         const SymmetricMatrix      *b             ,
                                                                               RealArray2D          *mA            ,
                                                                               RealArray2D          *mB            ,
                                                                               RealArray2D          *mC            ,
                                                                               Status               *status        ) ;
extern void                 AntisymmetricMatrix_CommutatorSS_Reference (       AntisymmetricMatrix  *self          ,
                                                                         const SymmetricMatrix      *a             ,
                                                                         const SymmetricMatrix      *b             ,
                                                                               Status               *status        ) ;
extern void                 AntisymmetricMatrix_CommutatorSSS          (       AntisymmetricMatrix  *self          ,
                                                                         const SymmetricMatrix      *a             ,
                                                                         const SymmetricMatrix      *b             ,
                                                                         const SymmetricMatrix      *c             ,
                                                                               Status               *status        ) ;
extern void                 AntisymmetricMatrix_CommutatorTSSST        (       AntisymmetricMatrix  *self          ,
                                                                         const SymmetricMatrix      *a             ,
                                                                         const SymmetricMatrix      *b             ,
                                                                         const SymmetricMatrix      *c             ,
                                                                         const RealArray2D          *m             ,
                                                                         const Boolean               mTranspose    ,
                                                                               RealArray2D          *u             ,
                                                                               RealArray2D          *v             ,
                                                                               RealArray2D          *w             ,
                                                                               Status               *status        ) ;
extern void                 AntisymmetricMatrix_CommutatorXSSY         (       AntisymmetricMatrix  *self          ,
                                                                         const SymmetricMatrix      *a             ,
                                                                         const SymmetricMatrix      *b             ,
                                                                         const RealArray2D          *x             ,
                                                                         const RealArray2D          *y             ,
                                                                         const Boolean               xTranspose    ,
                                                                         const Boolean               yTranspose    ,
                                                                               RealArray2D          *u             ,
                                                                               RealArray2D          *v             ,
                                                                               RealArray2D          *w             ,
                                                                               Status               *status        ) ;
extern void                 AntisymmetricMatrix_CopyFromRealArray2D    (       AntisymmetricMatrix  *self          ,
                                                                         const RealArray2D          *other         ,
                                                                         const Boolean               scale         ,
                                                                               Status               *status        ) ;
extern void                 AntisymmetricMatrix_CopyTo                 ( const AntisymmetricMatrix  *self          ,
                                                                               AntisymmetricMatrix  *other         ,
                                                                               Status               *status        ) ;
extern void                 AntisymmetricMatrix_CopyToRealArray2D      ( const AntisymmetricMatrix  *self          ,
                                                                               RealArray2D          *other         ,
                                                                               Status               *status        ) ;
extern void                 AntisymmetricMatrix_Deallocate             (       AntisymmetricMatrix **self          ) ;
extern AntisymmetricMatrix *AntisymmetricMatrix_FromExtentBlock        ( const Integer               extent        ,
                                                                               RealBlock            *block         ,
                                                                         const Boolean               withReference ,
                                                                               Status               *status        ) ;
extern void                 AntisymmetricMatrix_GetColumn              ( const AntisymmetricMatrix  *self          ,
                                                                         const Integer               n             ,
                                                                               RealArray1D          *column        ,
                                                                               Status               *status        ) ;
extern Real                 AntisymmetricMatrix_GetItem                ( const AntisymmetricMatrix  *self          ,
                                                                         const Integer               i             ,
                                                                         const Integer               j             ,
                                                                               Status               *status        ) ;
extern Boolean              AntisymmetricMatrix_GetItemIndexAndSign    ( const AntisymmetricMatrix  *self          ,
                                                                         const Integer               i             ,
                                                                         const Integer               j             ,
                                                                               Integer              *index         ,
                                                                               Real                 *sign          ,
                                                                               Status               *status        ) ;
extern Iterator            *AntisymmetricMatrix_MakeIterator           ( const AntisymmetricMatrix  *self          ,
                                                                               Status               *status        ) ;
extern void                 AntisymmetricMatrix_Print                  ( const AntisymmetricMatrix  *self          ) ;
extern void                 AntisymmetricMatrix_Scale                  (       AntisymmetricMatrix  *self          ,
                                                                         const Real                  value         ) ;
extern void                 AntisymmetricMatrix_Set                    (       AntisymmetricMatrix  *self          ,
                                                                         const Real                  value         ) ;
extern void                 AntisymmetricMatrix_SetItem                (       AntisymmetricMatrix  *self          ,
                                                                         const Integer               i             ,
                                                                         const Integer               j             ,
                                                                         const Real                  value         ,
                                                                               Status               *status        ) ;
extern void                 AntisymmetricMatrix_SymmetricTransform     ( const AntisymmetricMatrix  *self          ,
                                                                         const SymmetricMatrix      *matrix        ,
                                                                               AntisymmetricMatrix  *result        ,
                                                                               Status               *status        ) ;
extern Real                 AntisymmetricMatrix_TraceOfProduct         ( const AntisymmetricMatrix  *self          ,
                                                                         const AntisymmetricMatrix  *other         ,
                                                                               Status               *status        ) ;
extern void                 AntisymmetricMatrix_Transform              ( const AntisymmetricMatrix  *self          ,
                                                                         const RealArray2D          *matrix        ,
                                                                         const Boolean               useTranspose  ,
                                                                               AntisymmetricMatrix  *result        ,
                                                                               Status               *status        ) ;
extern void                 AntisymmetricMatrix_Transpose              (       AntisymmetricMatrix  *self          ) ;
extern Integer              AntisymmetricMatrix_ViewSize               ( const Integer               extent        ) ;

# endif
