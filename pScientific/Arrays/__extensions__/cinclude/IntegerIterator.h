# ifndef _INTEGERITERATOR
# define _INTEGERITERATOR

# include "Boolean.h"
# include "Integer.h"
# include "Iterator.h"
# include "Real.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Unary. */
extern Integer IntegerIterator_AbsoluteMaximum   ( Iterator *self, const Integer *selfData                                                               ) ;
extern Integer IntegerIterator_CountSmall        ( Iterator *self, const Integer *selfData, const Integer tolerance                                      ) ;
extern Integer IntegerIterator_DotSelf           ( Iterator *self, const Integer *selfData                                                               ) ;
extern void    IntegerIterator_FilterGreaterThan ( Iterator *self,       Integer *selfData, const Integer tolerance, const Integer value, Status *status ) ;
extern void    IntegerIterator_FilterLessThan    ( Iterator *self,       Integer *selfData, const Integer tolerance, const Integer value, Status *status ) ;
extern void    IntegerIterator_FilterSmall       ( Iterator *self,       Integer *selfData, const Integer tolerance, const Integer value, Status *status ) ;
extern void    IntegerIterator_Increment         ( Iterator *self,       Integer *selfData,                          const Integer value, Status *status ) ;
extern Integer IntegerIterator_Maximum           ( Iterator *self, const Integer *selfData                                                               ) ;
extern Integer IntegerIterator_Minimum           ( Iterator *self, const Integer *selfData                                                               ) ;
extern Integer IntegerIterator_Product           ( Iterator *self, const Integer *selfData                                                               ) ;
extern void    IntegerIterator_Scale             ( Iterator *self,       Integer *selfData,                          const Integer value, Status *status ) ; 
extern void    IntegerIterator_Set               ( Iterator *self,       Integer *selfData,                          const Integer value, Status *status ) ; 
extern Real    IntegerIterator_Sparsity          ( Iterator *self, const Integer *selfData, const Integer tolerance                                      ) ;
extern Integer IntegerIterator_Sum               ( Iterator *self, const Integer *selfData                                                               ) ;

/* . Binary. */
extern void    IntegerIterator_Add               ( Iterator *self,       Integer *selfData, Iterator *other, const Integer *otherData, const Integer scale, Status *status ) ;
extern void    IntegerIterator_CopyTo            ( Iterator *self, const Integer *selfData, Iterator *other,       Integer *otherData,                      Status *status ) ;
extern Integer IntegerIterator_Dot               ( Iterator *self, const Integer *selfData, Iterator *other, const Integer *otherData,                      Status *status ) ;
extern void    IntegerIterator_Multiply          ( Iterator *self,       Integer *selfData, Iterator *other, const Integer *otherData,                      Status *status ) ;
extern void    IntegerIterator_Swap              ( Iterator *self,       Integer *selfData, Iterator *other,       Integer *otherData,                      Status *status ) ;

# endif
