# ifndef _ITERATOR
# define _ITERATOR

# include "Boolean.h"
# include "Integer.h"
# include "Status.h"

/* . Notes:

  Inner loops or loops are sequences of elements separated by a constant stride.

    extent              maximum extent of a loop
    numberOfLoops       number of loops required to iterate over all items
    size                number of items in iterator

  An iterator is regular if the extents of all loops are the same and extent > 1, otherwise it is irregular.

*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The iterator type. */
typedef struct {
    void    *( * Clone         ) ( void  *vSelf , Status *status ) ;
    Integer  ( * CurrentIndex  ) ( void  *vSelf ) ;
    Integer  ( * DataOffSet    ) ( void  *vSelf ) ;
    void     ( * Deallocate    ) ( void **vSelf ) ;
    Integer *( * Dump          ) ( void  *vSelf , const Integer n0, Integer *n, Status *status ) ;
    void    *( * Load          ) ( const Integer n0, const Integer n, const Integer *state, Status *status ) ;
    Integer  ( * NextIndex     ) ( void  *vSelf ) ;
    Boolean  ( * NextInnerLoop ) ( void  *vSelf , Integer *first, Integer *extent, Integer *stride ) ;
    void     ( * Reset         ) ( void  *vSelf ) ;
} IteratorType ;

/* . The iterator. */
typedef struct {
/*          Boolean       readOnly      ; */
          Boolean       isRegular     ;
          Integer       extent        ;
          Integer       numberOfLoops ;
          Integer       size          ;
    const IteratorType *type          ;
          void         *vSelf         ;
} Iterator ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern       Iterator     *Iterator_Allocate        (       Status       *status ) ;
extern       Iterator     *Iterator_Clone           ( const Iterator     *self   ,
                                                            Status       *status ) ;
extern       Integer       Iterator_CurrentIndex    ( const Iterator     *self   ) ;
extern       Integer       Iterator_DataOffSet      ( const Iterator     *self   ) ;
extern       void          Iterator_Deallocate      (       Iterator    **self   ) ;
extern       Integer      *Iterator_Dump            ( const Iterator     *self   ,
                                                            Integer      *n      ,
                                                            Status       *status ) ;
extern       Integer       Iterator_GetSize         ( const Iterator     *self   ) ;
extern       Iterator     *Iterator_Load            ( const Integer       n      ,
                                                      const Integer      *state  ,
                                                            Status       *status ) ;
extern       Integer       Iterator_NextIndex       (       Iterator     *self   ) ;
extern       void          Iterator_Reset           (       Iterator     *self   ) ;
extern const IteratorType *Iterator_TypeFromInteger ( const Integer       n      ,
                                                            Status       *status ) ;
extern       Integer       Iterator_TypeToInteger   ( const IteratorType *type   ,
                                                            Status       *status ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Type declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Normal antisymmetric, double symmetric and symmetric iterators are Regular1D. */
extern const IteratorType IteratorType_Regular1D ;
extern const IteratorType IteratorType_RegularND ;
extern const IteratorType IteratorType_Row2D     ;

# endif
