# ifndef _LEBEDEVLAIKOV
# define _LEBEDEVLAIKOV

# include "Coordinates3.h"
# include "Integer.h"
# include "Real.h"
# include "RealArray1D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Integer LebedevLaikov_Angular_Momentum_Value ( const Integer       numberOfPoints ) ;
extern void    LebedevLaikov_GridPointsWeights      ( const Integer       numberOfPoints ,
                                                            Coordinates3 *gridPoints     ,
                                                            RealArray1D  *weights        ) ;
extern Integer LebedevLaikov_Number_Of_Points       ( const Integer       lValue         ) ;
extern Integer LebedevLaikov_Points                 ( const Integer       N              ,
                                                            Real         *X              ,
                                                            Real         *Y              ,
                                                            Real         *Z              ,
                                                            Real         *W              ) ;

# endif
