# ifndef _QCDISPERSIONDFTD2
# define _QCDISPERSIONDFTD2

# include "Coordinates3.h"
# include "Real.h"
# include "RealArray1D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Real QCDispersionDFTD2_Energy ( const Real          s6           ,
                                       const Real          sR           ,
                                       const Real          dR           ,
                                       const RealArray1D  *sqrtC6       ,
                                       const RealArray1D  *r0           ,
                                       const Coordinates3 *coordinates3 ,
                                             Coordinates3 *gradients3   ) ;

# endif
