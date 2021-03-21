# ifndef _VIEW2D
# define _VIEW2D

# include "Boolean.h"
# include "Integer.h"
# include "Iterator.h"
# include "Selection.h"
# include "Slice.h"
# include "Status.h"
# include "View1D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Content macros.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The view fields. */
# define View2D_Fields \
    Integer extent0 ; \
    Integer extent1 ; \
    Integer offset  ; \
    Integer size    ; \
    Integer stride0 ; \
    Integer stride1

/* . Field initialization. */
# define View2D_InitializeFields( self ) \
    (self)->extent0 = 0 ; \
    (self)->extent1 = 0 ; \
    (self)->offset  = 0 ; \
    (self)->size    = 0 ; \
    (self)->stride0 = 1 ; \
    (self)->stride1 = 1

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The view type. */
typedef struct {
    View2D_Fields ;
} View2D ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Conformability of two views. */
# define View2D_AreConformable( self, other ) ( ( (self)->extent0 == (other)->extent0 ) && ( (self)->extent1 == (other)->extent1 ) )

/* . The number of columns. */
# define View2D_Columns( self ) ( ( (self) == NULL ) ? 0 : (self)->extent1 )

/* . Tests. */
# define View2D_IsCompact0( self ) ( ( (self) != NULL ) && ( (self)->stride0 == (self)->extent1 ) )
# define View2D_IsCompact1( self ) ( ( (self) != NULL ) && ( (self)->stride1 == 1               ) )
# define View2D_IsCompact(  self ) ( ( (self) != NULL ) && ( (self)->stride0 == (self)->extent1 ) && ( (self)->stride1 == 1 ) )
# define View2D_IsSquare(   self ) ( ( (self) != NULL ) && ( (self)->extent0 == (self)->extent1 ) )
# define View2D_IsUniform(  self ) ( ( (self) != NULL ) && ( (self)->stride0 == (self)->extent1 * (self)->stride1 ) )

/* . An item index. */
# define View2D_ItemIndex( self, i, j ) ( (i)*(self)->stride0+(j)*(self)->stride1 )

/* . The number of rows. */
# define View2D_Rows( self ) ( ( (self) == NULL ) ? 0 : (self)->extent0 )

/* . The size of the view. */
# define View2D_Size( self ) ( (self) == NULL ? 0 : (self)->size )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern View2D   *View2D_Allocate            (       Status     *status     ) ;
extern Boolean   View2D_CheckCapacity       ( const View2D     *self       ,
                                              const Integer     capacity   ) ;
extern Boolean   View2D_CheckConformability ( const View2D     *self       ,
                                              const View2D     *other      ,
                                                    Status     *status     ) ;
extern void      View2D_ColumnView          ( const View2D     *self       ,
                                              const Integer     i          ,
                                                    View1D     *view       ,
                                                    Status     *status     ) ;
extern void      View2D_CopyTo              ( const View2D     *self       ,
                                                    View2D     *other      ) ;
extern void      View2D_Deallocate          (       View2D    **self       ) ;
extern Integer   View2D_GetColumns          ( const View2D     *self       ) ;
extern Integer   View2D_GetOffset           ( const View2D     *self       ) ;
extern Integer   View2D_GetRows             ( const View2D     *self       ) ;
extern Integer   View2D_GetSize             ( const View2D     *self       ) ;
extern Integer   View2D_GetStride           ( const View2D     *self       ,
                                              const Integer     dimension  ,
                                                    Status     *status     ) ;
extern void      View2D_Initialize          (       View2D     *self       ,
                                                    Integer     rows       ,
                                                    Integer     columns    ,
                                                    Status     *status     ) ;
extern Iterator *View2D_MakeIterator        ( const View2D     *self       ,
                                                    Status     *status     ) ;
extern Iterator *View2D_MakeRowIterator     ( const View2D     *self       ,
                                              const Selection  *rows       ,
                                                    Status     *status     ) ;
extern void      View2D_RowView             ( const View2D     *self       ,
                                              const Integer     i          ,
                                                    View1D     *view       ,
                                                    Status     *status     ) ;
extern void      View2D_SetState            (       View2D     *self       ,
                                              const Integer     extent0    ,
                                              const Integer     extent1    ,
                                              const Integer     offset     ,
                                              const Integer     size       ,
                                              const Integer     stride0    ,
                                              const Integer     stride1    ) ;
extern void      View2D_View                ( const View2D     *self       ,
                                              const Integer     start0     ,
                                              const Integer     start1     ,
                                              const Integer     extent0    ,
                                              const Integer     extent1    ,
                                              const Integer     stride0    ,
                                              const Integer     stride1    ,
                                                    View2D     *view       ,
                                                    Status     *status     ) ;
extern void      View2D_View1D              ( const View2D     *self       ,
                                              const Integer     dimension  ,
                                              const Integer     start0     ,
                                              const Integer     start1     ,
                                              const Integer     extent     ,
                                              const Integer     stride     ,
                                                    View1D     *view       ,
                                                    Status     *status     ) ;
extern void      View2D_View1DMultiSlice    ( const View2D     *self       ,
                                              const MultiSlice *multiSlice ,
                                                    View1D     *view       ,
                                                    Status     *status     ) ;
extern void      View2D_ViewMultiSlice      ( const View2D     *self       ,
                                              const MultiSlice *multiSlice ,
                                                    View2D     *view       ,
                                                    Status     *status     ) ;
# endif
