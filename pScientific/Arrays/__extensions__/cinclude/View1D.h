# ifndef _VIEW1D
# define _VIEW1D

# include "Boolean.h"
# include "Integer.h"
# include "Iterator.h"
# include "Slice.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Content macros.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The view fields. */
# define View1D_Fields \
    Integer extent ; \
    Integer offset ; \
    Integer size   ; \
    Integer stride ;

/* . Field initialization. */
# define View1D_InitializeFields( self ) \
    (self)->extent = 0 ; \
    (self)->offset = 0 ; \
    (self)->size   = 0 ; \
    (self)->stride = 1

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The view type. */
typedef struct {
    View1D_Fields ;
} View1D ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Conformability of two views. */
# define View1D_AreConformable( self, other ) ( (self)->extent == (other)->extent )

/* . The extent. */
# define View1D_Extent( self ) ( ( (self) == NULL ) ? 0 : (self)->extent )

/* . Tests. */
# define View1D_IsCompact( self ) ( ( (self) != NULL ) && ( (self)->stride == 1 ) )
# define View1D_IsUniform( self ) (   (self) != NULL )

/* . An item index into data (without offset). */
# define View1D_ItemIndex( self, i ) ( (i)*(self)->stride )

/* . The size of the view. */
# define View1D_Size( self ) ( (self) == NULL ? 0 : (self)->size )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern View1D   *View1D_Allocate            (       Status   *status   ) ;
extern Boolean   View1D_CheckCapacity       ( const View1D   *self     ,
                                              const Integer   capacity ) ;
extern Boolean   View1D_CheckConformability ( const View1D   *self     ,
                                              const View1D   *other    ,
                                                    Status   *status   ) ;
extern void      View1D_CopyTo              ( const View1D   *self     ,
                                                    View1D   *other    ) ;
extern void      View1D_Deallocate          (       View1D  **self     ) ;
extern Integer   View1D_GetExtent           ( const View1D   *self     ) ;
extern Integer   View1D_GetOffset           ( const View1D   *self     ) ;
extern Integer   View1D_GetStride           ( const View1D   *self     ) ;
extern void      View1D_Initialize          (       View1D   *self     ,
                                                    Integer   extent   ,
                                                    Status   *status   ) ;
extern Iterator *View1D_MakeIterator        ( const View1D   *self     ,
                                                    Status   *status   ) ;
extern void      View1D_SetState            (       View1D   *self     ,
                                              const Integer   extent   ,
                                              const Integer   offset   ,
                                              const Integer   size     ,
                                              const Integer   stride   ) ;
extern void      View1D_View                ( const View1D   *self     ,
                                              const Integer   start    ,
                                              const Integer   extent   ,
                                              const Integer   stride   ,
                                                    View1D   *view     ,
                                                    Status   *status   ) ;
# endif
