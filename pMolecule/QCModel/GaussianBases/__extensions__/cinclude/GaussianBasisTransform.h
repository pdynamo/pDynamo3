# ifndef _GAUSSIANBASISTRANSFORM
# define _GAUSSIANBASISTRANSFORM

# include "Boolean.h"
# include "Integer.h"
# include "Real.h"
# include "RealArray2D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void GaussianBasisTransform1  ( const RealArray2D *tI         ,
                                             Real       **values     ,
                                             Real       **work       ) ;
extern void GaussianBasisTransform1M (       Integer      dI         ,
                                             Integer      dM         ,
                                       const RealArray2D *tI         ,
                                       const Boolean      transpose  ,
                                             Real       **values     ,
                                             Real       **work       ) ;
extern void GaussianBasisTransform2  (       Integer      dI         ,
                                             Integer      dJ         ,
                                       const RealArray2D *tI         ,
                                       const RealArray2D *tJ         ,
                                             Real       **values     ,
                                             Real       **work       ) ;
extern void GaussianBasisTransform3  (       Integer      dI         ,
                                             Integer      dJ         ,
                                             Integer      dK         ,
                                       const RealArray2D *tI         ,
                                       const RealArray2D *tJ         ,
                                       const RealArray2D *tK         ,
                                             Real       **values     ,
                                             Real       **work       ) ;
extern void GaussianBasisTransform4  (       Integer      dI         ,
                                             Integer      dJ         ,
                                             Integer      dK         ,
                                             Integer      dL         ,
                                       const RealArray2D *tI         ,
                                       const RealArray2D *tJ         ,
                                       const RealArray2D *tK         ,
                                       const RealArray2D *tL         ,
                                             Real       **values     ,
                                             Real       **work       ) ;
# endif
