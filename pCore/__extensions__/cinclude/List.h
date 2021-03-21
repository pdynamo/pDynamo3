# ifndef _LIST
# define _LIST

# include "Boolean.h"
# include "Integer.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The list element type. */
struct ListElementStructure {
    struct ListElementStructure *next ;
    void                        *node ;
} ;
typedef struct ListElementStructure ListElement ;

/* . The list type. */
typedef struct {
    Integer      nelements ;
    ListElement *first     ;
    ListElement *iterator  ;
    ListElement *last      ;
    Boolean     ( * Element_Match      ) ( void *reference, void *current );
    void        ( * Element_Deallocate ) ( void *current );
} List ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern List        *List_Allocate                ( void ) ;
extern void         List_Deallocate              (       List **self ) ;
extern void         List_Element_Append          (       List  *self , void *node ) ;
extern void         List_Element_Append_By_Index (       List  *self , void *node, const Integer index ) ;
extern void        *List_Element_Find_By_Index   (       List  *self , const Integer index ) ;
extern void        *List_Element_Find_By_Match   (       List  *self , void *reference     ) ;
extern void        *List_Element_Pop_By_Index    (       List  *self , const Integer index ) ;
extern void        *List_Element_Pop_By_Match    (       List  *self , void *reference     ) ;
extern void         List_Empty                   (       List  *self ) ;
extern void         List_Initialize              (       List  *self ) ;
extern void        *List_Iterate                 (       List  *self ) ;
extern ListElement *List_Iterate_Current         ( const List  *self ) ;
extern void         List_Iterate_Initialize      (       List  *self ) ;
extern void         List_Iterate_Set             (       List  *self , ListElement *iterator ) ;
extern Integer      List_Size                    ( const List  *self ) ;

# endif
