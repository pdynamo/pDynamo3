# ifndef _ROTATEORBITAL
# define _ROTATEORBITAL

# include "IntegerArray1D.h"
# include "Matrix33.h"
# include "RealArray1D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void RotateOrbital                (       IntegerArray1D  *orbitalBasisIndices ,
                                           const Matrix33        *rotation            ,
                                           const IntegerArray1D  *mapping             ,
                                           const RealArray1D     *inOrbital           ,
                                                 RealArray1D     *outOrbital          ) ;
extern void RotateOrbital_MakeLRotations ( const Integer          L                   ,
                                           const Matrix33        *R                   ,
                                                 RealArray2D     *T                   ,
                                                 Status          *status              ) ;
# endif
