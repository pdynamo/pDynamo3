# ifndef _CORRELATION
# define _CORRELATION

# include "Boolean.h"
# include "Integer.h"
# include "Real.h"
# include "RealArray1D.h"
# include "RealArray2D.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void Correlation_MakeDotProduct (       RealArray2D *x            ,
                                               RealArray2D *y            ,
                                         const Boolean      useFFT       ,
                                         const Boolean      normalize    ,
                                         const Boolean      removeMean   ,
                                         const Integer      tCorrelation ,
                                         const Real        *tolerance    ,
                                               RealArray1D *f            ,
                                               Status      *status       ) ;
extern void Correlation_MakeSimple     (       RealArray1D *x            ,
                                               RealArray1D *y            ,
                                         const Boolean      useFFT       ,
                                         const Boolean      normalize    ,
                                         const Boolean      removeMean   ,
                                         const Integer      tCorrelation ,
                                         const Real        *tolerance    ,
                                               RealArray1D *f            ,
                                               Status      *status       ) ;

# endif
