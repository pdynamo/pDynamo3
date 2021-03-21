# ifndef _ITERATORND
# define _ITERATORND

# include "Boolean.h"
# include "Integer.h"
# include "Iterator.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . N-D iterator. */
typedef struct {
    Integer  next     ;
    Integer  offset   ;
    Integer  rank     ;
    Integer  size     ;
    Integer *counters ;
    Integer *extents  ;
    Integer *strides  ;
} IteratorND ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern IteratorND *IteratorND_Allocate      ( const Integer     rank     ,
                                                    Status     *status   ) ;
extern void       *IteratorND_Clone         (       void       *vSelf    ,
                                                    Status     *status   ) ;
extern Integer     IteratorND_CurrentIndex  (       void       *vSelf    ) ;
extern Integer     IteratorND_DataOffSet    (       void       *vSelf    ) ;
extern void        IteratorND_Deallocate    (       void      **vSelf    ) ;
extern Integer    *IteratorND_Dump          (       void       *vSelf    ,
                                              const Integer     n0       ,
                                                    Integer    *n        ,
                                                    Status     *status   ) ;
extern void        IteratorND_Finalize      (       IteratorND *self     ) ;
extern void        IteratorND_Initialize    (       IteratorND *self     ,
                                              const Integer     offset   ,
                                              const Integer    *extents  ,
                                              const Integer    *strides  ,
                                                    Status     *status   ) ;
extern void       *IteratorND_Load          ( const Integer     n0       ,
                                              const Integer     n        ,
                                              const Integer    *state    ,
                                                    Status     *status   ) ;
extern void        IteratorND_MakeIterator  ( const IteratorND *self     ,
                                                    Iterator   *iterator ) ;
extern Integer     IteratorND_NextIndex     (       void       *vSelf    ) ;
extern Boolean     IteratorND_NextInnerLoop (       void       *vSelf    ,
                                                    Integer    *first    ,
                                                    Integer    *extent   ,
                                                    Integer    *stride   ) ;
extern void        IteratorND_Reset         (       void       *vSelf    ) ;

# endif
