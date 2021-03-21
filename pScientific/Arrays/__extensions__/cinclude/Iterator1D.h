# ifndef _ITERATOR1D
# define _ITERATOR1D

# include "Boolean.h"
# include "Integer.h"
# include "Iterator.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . 1-D iterator. */
typedef struct {
    Integer counter ;
    Integer extent  ;
    Integer offset  ;
    Integer next    ;
    Integer stride  ;
} Iterator1D ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Iterator1D *Iterator1D_Allocate      (       Status     *status   ) ;
extern void       *Iterator1D_Clone         (       void       *vSelf    ,
                                                    Status     *status   ) ;
extern Integer     Iterator1D_CurrentIndex  (       void       *vSelf    ) ;
extern Integer     Iterator1D_DataOffSet    (       void       *vSelf    ) ;
extern void        Iterator1D_Deallocate    (       void      **vSelf    ) ;
extern Integer    *Iterator1D_Dump          (       void       *vSelf    ,
                                              const Integer     n0       ,
                                                    Integer    *n        ,
                                                    Status     *status   ) ;
extern void        Iterator1D_Initialize    (       Iterator1D *self     ,
                                              const Integer     offset   ,
                                              const Integer     extent   ,
                                              const Integer     stride   ,
                                                    Status     *status   ) ;
extern void       *Iterator1D_Load          ( const Integer     n0       ,
                                              const Integer     n        ,
                                              const Integer    *state    ,
                                                    Status     *status   ) ;
extern void        Iterator1D_MakeIterator  ( const Iterator1D *self     ,
                                                    Iterator   *iterator ) ;
extern Integer     Iterator1D_NextIndex     (       void       *vSelf    ) ;
extern Boolean     Iterator1D_NextInnerLoop (       void       *vSelf    ,
                                                    Integer    *first    ,
                                                    Integer    *extent   ,
                                                    Integer    *stride   ) ;
extern void        Iterator1D_Reset         (       void       *vSelf    ) ;

# endif
