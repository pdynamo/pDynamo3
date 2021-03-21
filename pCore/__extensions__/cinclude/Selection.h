# ifndef _SELECTION
# define _SELECTION

# include "Boolean.h"
# include "BooleanBlock.h"
# include "Integer.h"
# include "IntegerBlock.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The selection type. */
typedef struct {
    Integer       capacity  ;
    Integer      *indices   ;
    BooleanBlock *flags     ;
    IntegerBlock *positions ;
} Selection ;


/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros - no checking.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The items. */
# define Selection_Items( self ) ( (self)->indices )

/* . An item. */
# define Selection_Item( self, i ) ( (self)->indices[i] )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Selection    *Selection_Allocate             ( const Integer     capacity   ,
                                                            Status     *status     ) ;
extern Integer       Selection_Capacity             ( const Selection  *self       ) ;
extern void          Selection_ClearFlags           (       Selection  *self       ) ;
extern void          Selection_ClearPositions       (       Selection  *self       ) ;
extern void          Selection_ClearRepresentations (       Selection  *self       ) ;
extern Selection    *Selection_Clone                ( const Selection  *self       ,
                                                            Status     *status     ) ;
extern Selection    *Selection_Complement           (       Selection  *self       ,
                                                      const Integer     upperBound ,
                                                            Status     *status     ) ;
extern void          Selection_Deallocate           (       Selection **self       ) ;
extern Selection    *Selection_Difference           (       Selection  *self       ,
                                                      const Integer     number     ,
                                                            Selection **others     ,
                                                            Status     *status     ) ;
extern Selection    *Selection_FromBooleans         ( const Integer     capacity   ,
                                                      const Boolean    *flags      ,
                                                            Status     *status     ) ;
extern Selection    *Selection_FromIntegers         ( const Integer     capacity   ,
                                                      const Integer    *indices    ,
                                                            Status     *status     ) ;
extern Boolean       Selection_HasItem              (       Selection  *self       ,
                                                      const Integer     value      ,
                                                            Status     *status     ) ;
extern Selection    *Selection_Intersection         ( const Integer     number     ,
                                                            Selection **others     ,
                                                            Status     *status     ) ;
extern Boolean       Selection_IsEmpty              ( const Selection  *self       ,
                                                      const Boolean     nullIsFull ) ;
extern BooleanBlock *Selection_MakeFlags            (       Selection  *self       ,
                                                      const Integer     upperBound ,
                                                            Status     *status     ) ;
extern IntegerBlock *Selection_MakePositions        (       Selection  *self       ,
                                                      const Integer     upperBound ,
                                                            Status     *status     ) ;
extern Integer       Selection_PositionOfItem       (       Selection  *self       ,
                                                      const Integer     value      ,
                                                            Status     *status     ) ;
extern Selection    *Selection_Prune                (       Selection  *self       ,
                                                            Selection  *toKeep     ,
                                                            Status     *status     ) ;
extern Selection    *Selection_SymmetricDifference  ( const Integer     number     ,
                                                            Selection **others     ,
                                                            Status     *status     ) ;
extern Selection    *Selection_Union                ( const Integer     number     ,
                                                            Selection **others     ,
                                                            Status     *status     ) ;
extern Integer       Selection_UpperBound           ( const Selection  *self       ) ;

# endif
