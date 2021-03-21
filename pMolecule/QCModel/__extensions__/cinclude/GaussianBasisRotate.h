# ifndef _GAUSSIANBASISROTATE
# define _GAUSSIANBASISROTATE

# include "Boolean.h"
# include "GaussianBasis.h"
# include "Integer.h"
# include "Matrix33.h"
# include "RealArray2D.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void GaussianBasis_MakeLRotations     ( const Integer        L      ,
                                               const Matrix33      *R      ,
                                                     RealArray2D   *Tc     ,
                                                     Status        *status ) ;
extern void GaussianBasis_MakeRotationMatrix (       GaussianBasis *self   ,
                                               const RealArray2D   *Tc     ,
                                               const Boolean        doC2O  ,
                                                     RealArray2D   *T      ,
                                                     Status        *status ) ;
# endif
