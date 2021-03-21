# ifndef _ROWITERATOR2D
# define _ROWITERATOR2D

# include "Boolean.h"
# include "Integer.h"
# include "Iterator.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . N-D iterator. */
typedef struct {
    Integer  counter0 ;
    Integer  counter1 ;
    Integer  extent0  ;
    Integer  extent1  ;
    Integer  next     ;
    Integer  offset   ;
    Integer  size     ;
    Integer  stride0  ;
    Integer  stride1  ;
    Integer *rows     ;
} RowIterator2D ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern RowIterator2D *RowIterator2D_Allocate      ( const Integer        extent0  ,
                                                          Status        *status   ) ;
extern void          *RowIterator2D_Clone         (       void          *vSelf    ,
                                                          Status        *status   ) ;
extern Integer        RowIterator2D_CurrentIndex  (       void          *vSelf    ) ;
extern Integer        RowIterator2D_DataOffSet    (       void          *vSelf    ) ;
extern void           RowIterator2D_Deallocate    (       void         **vSelf    ) ;
extern Integer       *RowIterator2D_Dump          (       void          *vSelf    ,
                                                    const Integer        n0       ,
                                                          Integer       *n        ,
                                                          Status        *status   ) ;
extern void           RowIterator2D_Finalize      (       RowIterator2D *self     ) ;
extern void           RowIterator2D_Initialize    (       RowIterator2D *self     ,
                                                    const Integer        extent1  ,
                                                    const Integer        offset   ,
                                                    const Integer        stride0  ,
                                                    const Integer        stride1  ,
                                                    const Integer       *rows     ,
                                                          Status        *status   ) ;
extern void          *RowIterator2D_Load          ( const Integer        n0       ,
                                                    const Integer        n        ,
                                                    const Integer       *state    ,
                                                          Status        *status   ) ;
extern void           RowIterator2D_MakeIterator  ( const RowIterator2D *self     ,
                                                          Iterator      *iterator ) ;
extern Integer        RowIterator2D_NextIndex     (       void          *vSelf    ) ;
extern Boolean        RowIterator2D_NextInnerLoop (       void          *vSelf    ,
                                                          Integer       *first    ,
                                                          Integer       *extent   ,
                                                          Integer       *stride   ) ;
extern void           RowIterator2D_Reset         (       void          *vSelf    ) ;

# endif
