# ifndef _COSINEPARAMETER
# define _COSINEPARAMETER

# include "Integer.h"
# include "Real.h"

/*----------------------------------------------------------------------------------------------------------------------------------
!
! . A cosine parameter is of the form:
!
!   E = Sum_n c_n cos ( n alpha )
!
!   or
!
!   E = Sum_p c_p cos^p ( alpha )
!
!   The latter is used.
!
!---------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
typedef struct {
    Integer  nPowers           ;
    Integer  nTerms            ;
    Integer *periods           ;
    Real    *powerCoefficients ;
    Real    *termCoefficients  ;
} CosineParameter ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void CosineParameter_Allocate   (       CosineParameter *self   ,
                                         const Integer          nTerms ) ;
extern void CosineParameter_Clone      (       CosineParameter *self   ,
                                         const CosineParameter *other  ) ;
extern void CosineParameter_Deallocate (       CosineParameter *self ) ;
extern void CosineParameter_Initialize (       CosineParameter *self ) ;
extern void CosineParameter_MakePowers (       CosineParameter *self ) ;

# endif
