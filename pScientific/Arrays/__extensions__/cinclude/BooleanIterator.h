# ifndef _BOOLEANITERATOR
# define _BOOLEANITERATOR

# include "Boolean.h"
# include "Integer.h"
# include "Iterator.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Unary. */
extern Boolean BooleanIterator_All        ( Iterator *self, const Boolean *selfData                                      ) ;
extern Boolean BooleanIterator_Any        ( Iterator *self, const Boolean *selfData                                      ) ;
extern Integer BooleanIterator_CountFalse ( Iterator *self, const Boolean *selfData                                      ) ;
extern Integer BooleanIterator_CountTrue  ( Iterator *self, const Boolean *selfData                                      ) ;
extern void    BooleanIterator_Not        ( Iterator *self,       Boolean *selfData,                      Status *status ) ;
extern void    BooleanIterator_Set        ( Iterator *self,       Boolean *selfData, const Boolean value, Status *status ) ;

/* . Binary. */
extern void    BooleanIterator_And        ( Iterator *self,       Boolean *selfData, Iterator *other, const Boolean *otherData, Status *status ) ;
extern void    BooleanIterator_CopyTo     ( Iterator *self, const Boolean *selfData, Iterator *other,       Boolean *otherData, Status *status ) ;
extern void    BooleanIterator_Or         ( Iterator *self,       Boolean *selfData, Iterator *other, const Boolean *otherData, Status *status ) ;
extern void    BooleanIterator_Xor        ( Iterator *self,       Boolean *selfData, Iterator *other, const Boolean *otherData, Status *status ) ;
extern void    BooleanIterator_Swap       ( Iterator *self,       Boolean *selfData, Iterator *other,       Boolean *otherData, Status *status ) ;

# endif
