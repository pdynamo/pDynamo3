# ifndef _REALITERATOR
# define _REALITERATOR

# include "Boolean.h"
# include "Integer.h"
# include "Iterator.h"
# include "Real.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Unary. */
extern void    RealIterator_Absolute          ( Iterator *self,       Real *selfData                                                                           ) ;    
extern Real    RealIterator_AbsoluteMaximum   ( Iterator *self, const Real *selfData                                                                           ) ;    
extern Integer RealIterator_CountSmall        ( Iterator *self, const Real *selfData,                   const Real tolerance                                   ) ;    
extern Real    RealIterator_DotSelf           ( Iterator *self, const Real *selfData                                                                           ) ;    
extern void    RealIterator_Exponential       ( Iterator *self,       Real *selfData,                                                           Status *status ) ;    
extern void    RealIterator_FilterGreaterThan ( Iterator *self,       Real *selfData,                   const Real tolerance, const Real value, Status *status ) ;    
extern void    RealIterator_FilterLessThan    ( Iterator *self,       Real *selfData,                   const Real tolerance, const Real value, Status *status ) ;    
extern void    RealIterator_FilterSmall       ( Iterator *self,       Real *selfData,                   const Real tolerance, const Real value, Status *status ) ;    
extern void    RealIterator_Increment         ( Iterator *self,       Real *selfData,                                         const Real value, Status *status ) ;    
extern Real    RealIterator_Maximum           ( Iterator *self, const Real *selfData                                                                           ) ;                                                                   
extern Real    RealIterator_Minimum           ( Iterator *self, const Real *selfData                                                                           ) ; 
extern void    RealIterator_NaturalLogarithm  ( Iterator *self,       Real *selfData,                                                           Status *status ) ; 
extern Real    RealIterator_Norm2             ( Iterator *self, const Real *selfData                                                                           ) ; 
extern Boolean RealIterator_Normalize         ( Iterator *self,       Real *selfData,                   const Real tolerance,                   Status *status ) ;             
extern void    RealIterator_Power             ( Iterator *self,       Real *selfData, const Real power,                                         Status *status ) ;             
extern Real    RealIterator_Product           ( Iterator *self, const Real *selfData                                                                           ) ;             
extern void    RealIterator_Reciprocate       ( Iterator *self,       Real *selfData,                   const Real tolerance, const Real value, Status *status ) ;
extern void    RealIterator_ReciprocatePower  ( Iterator *self,       Real *selfData, const Real power, const Real tolerance, const Real value, Status *status ) ;
extern Real    RealIterator_RootMeanSquare    ( Iterator *self, const Real *selfData                                                                           ) ;
extern void    RealIterator_Scale             ( Iterator *self,       Real *selfData,                                         const Real value, Status *status ) ;
extern void    RealIterator_Set               ( Iterator *self,       Real *selfData,                                         const Real value, Status *status ) ;
extern Real    RealIterator_Sparsity          ( Iterator *self, const Real *selfData,                   const Real tolerance                                   ) ;
extern void    RealIterator_Square            ( Iterator *self,       Real *selfData,                                                           Status *status ) ;
extern void    RealIterator_SquareRoot        ( Iterator *self,       Real *selfData,                                                           Status *status ) ;
extern Real    RealIterator_Sum               ( Iterator *self, const Real *selfData                                                                           ) ;

/* . Binary. */
extern void    RealIterator_Add      ( Iterator *self,       Real *selfData, Iterator *other, const Real *otherData,                       const Real scale, Status *status ) ;
extern void    RealIterator_CopyTo   ( Iterator *self, const Real *selfData, Iterator *other,       Real *otherData,                                         Status *status ) ;
extern void    RealIterator_Divide   ( Iterator *self,       Real *selfData, Iterator *other, const Real *otherData, const Real tolerance, const Real value, Status *status ) ;
extern Real    RealIterator_Dot      ( Iterator *self, const Real *selfData, Iterator *other, const Real *otherData,                                         Status *status ) ;
extern void    RealIterator_Multiply ( Iterator *self,       Real *selfData, Iterator *other, const Real *otherData,                                         Status *status ) ;
extern void    RealIterator_Swap     ( Iterator *self,       Real *selfData, Iterator *other,       Real *otherData,                                         Status *status ) ;

# endif
