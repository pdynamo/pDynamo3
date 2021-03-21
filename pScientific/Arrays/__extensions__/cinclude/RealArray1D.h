# ifndef _REALARRAY1D
# define _REALARRAY1D

# include "Real.h"
# include "RealBlock.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Basic definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _ArrayDataType Real
# include "Array1D_Header.i"
# undef  _ArrayDataType

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Real RealArray1D_Norm2     ( const RealArray1D  *self          ) ;
extern void RealArray1D_Normalize (       RealArray1D  *self          ,
                                    const Real         *nullNormValue ,
                                          Status       *status        ) ;
# endif
