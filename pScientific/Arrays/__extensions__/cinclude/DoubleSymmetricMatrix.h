# ifndef _DOUBLESYMMETRICMATRIX
# define _DOUBLESYMMETRICMATRIX

# include "Boolean.h"
# include "Integer.h"
# include "Iterator.h"
# include "Real.h"
# include "RealBlock.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The array type. */
typedef struct {
    Integer    extent   ; /* . Extent of each dimension. */
    Integer    extent01 ; /* . Extent of first two and last two dimensions. */
    Integer    size     ; /* . Overall size. */
    RealBlock *block    ;
    Real      *data     ;
} DoubleSymmetricMatrix ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Access an item (i >= j, k >= l, ij >= kl). */
/*
# define IJIndex(i,j) ( ( (i) * ( i + 1 ) ) / 2 + j )
# define DoubleSymmetricMatrix_Item(i,j,k,l) IJIndex ( IJIndex ( i, j ), IJIndex ( k, l ) )
*/

/* . The extent. */
# define DoubleSymmetricMatrix_Extent( self ) ( (self) == NULL ? 0 : (self)->extent )

/* . The size. */
# define DoubleSymmetricMatrix_Size( self ) ( (self) == NULL ? 0 : (self)->size )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern DoubleSymmetricMatrix *DoubleSymmetricMatrix_Allocate           (       Status                 *status        ) ;
extern DoubleSymmetricMatrix *DoubleSymmetricMatrix_AllocateWithExtent ( const Integer                 extent        ,
                                                                               Status                 *status        ) ;
extern DoubleSymmetricMatrix *DoubleSymmetricMatrix_CloneDeep          ( const DoubleSymmetricMatrix  *self          ,
                                                                               Status                 *status        ) ;
extern DoubleSymmetricMatrix *DoubleSymmetricMatrix_CloneShallow       ( const DoubleSymmetricMatrix  *self          ,
                                                                               Status                 *status        ) ;
extern void                   DoubleSymmetricMatrix_CopyTo             ( const DoubleSymmetricMatrix  *self          ,
                                                                               DoubleSymmetricMatrix  *other         ,
                                                                               Status                 *status        ) ;
extern void                   DoubleSymmetricMatrix_Deallocate         (       DoubleSymmetricMatrix **self          ) ;
extern DoubleSymmetricMatrix *DoubleSymmetricMatrix_FromExtentBlock    ( const Integer                 extent        ,
                                                                               RealBlock              *block         ,
                                                                         const Boolean                 withReference ,
                                                                               Status                 *status        ) ;
extern Real                   DoubleSymmetricMatrix_GetItem            ( const DoubleSymmetricMatrix  *self          ,
                                                                         const Integer                 i             ,
                                                                         const Integer                 j             ,
                                                                         const Integer                 k             ,
                                                                         const Integer                 l             ,
                                                                               Status                 *status        ) ;
extern void                   DoubleSymmetricMatrix_IncrementItem      ( const DoubleSymmetricMatrix  *self          ,
                                                                         const Integer                 i             ,
                                                                         const Integer                 j             ,
                                                                         const Integer                 k             ,
                                                                         const Integer                 l             ,
                                                                         const Real                    value         ,
                                                                               Status                 *status        ) ;
extern Integer                DoubleSymmetricMatrix_Index              ( const Integer                 i             ,
                                                                         const Integer                 j             ,
                                                                         const Integer                 k             ,
                                                                         const Integer                 l             ) ;
extern Iterator              *DoubleSymmetricMatrix_MakeIterator       ( const DoubleSymmetricMatrix  *self          ,
                                                                               Status                 *status        ) ;
extern void                   DoubleSymmetricMatrix_Print              ( const DoubleSymmetricMatrix  *self          ) ;
extern void                   DoubleSymmetricMatrix_Set                (       DoubleSymmetricMatrix  *self          ,
                                                                               Real                    value         ) ;
extern void                   DoubleSymmetricMatrix_SetItem            ( const DoubleSymmetricMatrix  *self          ,
                                                                         const Integer                 i             ,
                                                                         const Integer                 j             ,
                                                                         const Integer                 k             ,
                                                                         const Integer                 l             ,
                                                                         const Real                    value         ,
                                                                               Status                 *status        ) ;
extern void                   DoubleSymmetricMatrix_Unweight           (       DoubleSymmetricMatrix  *self          ) ;
extern Integer                DoubleSymmetricMatrix_ViewSize           ( const Integer                 extent        ) ;

# endif
