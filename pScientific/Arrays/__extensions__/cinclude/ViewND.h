# ifndef _VIEWND
# define _VIEWND

# include "Boolean.h"
# include "Integer.h"
# include "Iterator.h"
# include "Slice.h"
# include "Status.h"
# include "View1D.h"
# include "View2D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The view type. */
typedef struct {
    Integer  offset  ;
    Integer  rank    ;
    Integer  size    ;
    Integer *extents ;
    Integer *strides ;
} ViewND ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Indexing - these can also be used on arrays of higher rank as all missing indices will take the value 0. */
# define ViewND_Index1D( self, i )          ( (i)*(self)->strides[0]  )
# define ViewND_Index2D( self, i, j )       ( (i)*(self)->strides[0]+ \
                                              (j)*(self)->strides[1]  )
# define ViewND_Index3D( self, i, j, k )    ( (i)*(self)->strides[0]+ \
                                              (j)*(self)->strides[1]+ \
                                              (k)*(self)->strides[2]  )
# define ViewND_Index4D( self, i, j, k, l ) ( (i)*(self)->strides[0]+ \
                                              (j)*(self)->strides[1]+ \
                                              (k)*(self)->strides[2]+ \
                                              (l)*(self)->strides[3]  )

/* . The view size. */
# define ViewND_Size( self ) ( ( (self) == NULL ) ? 0 : (self)->size )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern ViewND   *ViewND_Allocate           (       Status     *status     ) ;
extern ViewND   *ViewND_AllocateWithRank   ( const Integer     rank       ,
                                                   Status     *status     ) ;
extern ViewND   *ViewND_AllocateWithShape  ( const Integer     rank       ,
                                             const Integer    *shape      ,
                                                   Status     *status     ) ;
extern Boolean   ViewND_AreConformable     ( const ViewND     *self       ,
                                             const ViewND     *other      ,
                                                   Status     *status     ) ;
extern Boolean   ViewND_CheckCapacity      ( const ViewND     *self       ,
                                             const Integer     capacity   ) ;
extern ViewND   *ViewND_Clone              ( const ViewND     *self       ,
                                                   Status     *status     ) ;
extern void      ViewND_Deallocate         (       ViewND    **self       ) ;
extern Integer   ViewND_Flatten            ( const ViewND     *self       ,
                                                   Integer    *extents    ,
                                                   Integer    *strides    ) ;
extern ViewND   *ViewND_FromState          ( const Integer     rank       ,
                                             const Integer     offset     ,
                                             const Integer     size       ,
                                                   Status     *status     ) ;
extern Integer   ViewND_GetExtent          ( const ViewND     *self       ,
                                             const Integer     dimension  ,
                                                   Status     *status     ) ;
extern Integer   ViewND_GetIndex           ( const ViewND     *self       ,
                                             const Integer     rank       ,
                                             const Integer    *indices    ,
                                                   Status     *status     ) ;
extern Integer   ViewND_GetIndexMultiSlice ( const ViewND     *self       ,
                                             const MultiSlice *multiSlice ,
                                                   Status     *status     ) ;
extern Integer   ViewND_GetOffset          ( const ViewND     *self       ) ;
extern Integer   ViewND_GetRank            ( const ViewND     *self       ) ;
extern Integer   ViewND_GetSize            ( const ViewND     *self       ) ;
extern Integer   ViewND_GetStride          ( const ViewND     *self       ,
                                             const Integer     dimension  ,
                                                   Status     *status     ) ;
extern void      ViewND_Initialize         (       ViewND     *self       ,
                                             const Integer    *extents    ,
                                                   Status     *status     ) ;
extern Boolean   ViewND_IsCompact          ( const ViewND     *self       ) ;
extern Boolean   ViewND_IsUniform          ( const ViewND     *self       ) ;
extern Iterator *ViewND_MakeIterator       ( const ViewND     *self       ,
                                                   Status     *status     ) ;
extern void      ViewND_SetExtentStride    (       ViewND     *self       ,
                                             const Integer     dimension  ,
                                             const Integer     extent     ,
                                             const Integer     stride     ,
                                                   Status     *status     ) ;
extern void      ViewND_View1DMultiSlice   ( const ViewND     *self       ,
                                             const MultiSlice *multiSlice ,
                                                   View1D     *view       ,
                                                   Status     *status     ) ;
extern void      ViewND_View2DMultiSlice   ( const ViewND     *self       ,
                                             const MultiSlice *multiSlice ,
                                                   View2D     *view       ,
                                                   Status     *status     ) ;
extern ViewND   *ViewND_ViewMultiSlice     ( const ViewND     *self       ,
                                             const MultiSlice *multiSlice ,
                                                   Status     *status     ) ;
extern void      ViewND_ViewTail2D         ( const ViewND     *self       ,
                                             const Integer    *indices    ,
                                                   View2D     *view       ,
                                                   Status     *status     ) ;
# endif
