# ifndef _SELECTIONCONTAINER
# define _SELECTIONCONTAINER

# include "Boolean.h"
# include "BooleanBlock.h"
# include "Integer.h"
# include "Selection.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The container type. */
typedef struct {
    Integer     capacity ;
    Selection **items    ;
} SelectionContainer ;


/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern SelectionContainer *SelectionContainer_Allocate            ( const Integer              capacity ,
                                                                          Status              *status   ) ;
extern Integer             SelectionContainer_Capacity            ( const SelectionContainer  *self     ) ;
extern SelectionContainer *SelectionContainer_Clone               ( const SelectionContainer  *self     ,
                                                                          Status              *status   ) ;
extern void                SelectionContainer_Deallocate          (       SelectionContainer **self     ) ;
extern SelectionContainer *SelectionContainer_FromCapacity        ( const Integer              capacity ,
                                                                          Status              *status   ) ;
extern Boolean             SelectionContainer_FuseItems           (       SelectionContainer  *self     ,
                                                                          BooleanBlock        *toFuse   ,
                                                                          Status              *status   ) ;
extern BooleanBlock       *SelectionContainer_MakeMembershipFlags (       SelectionContainer  *self     ,
                                                                          Selection           *members  ,
                                                                    const Boolean              andTest  ,
                                                                          Status              *status   ) ;
extern Boolean             SelectionContainer_RemoveItems         (       SelectionContainer  *self     ,
                                                                          Selection           *toRemove ,
                                                                          Status              *status   ) ;
extern Selection          *SelectionContainer_UnionOfItems        ( const SelectionContainer  *self     ,
                                                                          Selection           *toUnion  ,
                                                                          Status              *status   ) ;
extern Integer             SelectionContainer_UpperBound          ( const SelectionContainer  *self     ) ;

# endif
