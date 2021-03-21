# ifndef _CIFOURINDEXTRANSFORMATION
# define _CIFOURINDEXTRANSFORMATION

# include "BlockStorage.h"
# include "DoubleSymmetricMatrix.h"
# include "RealArray2D.h"
# include "RealArrayND.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void CIFourIndexTransformation ( const RealArray2D           *activeMOs            ,
                                              BlockStorage          *twoElectronIntegrals ,
                                              RealArray2D           *moTEI34              ,
                                              RealArrayND           *moTEI234             ,
                                              DoubleSymmetricMatrix *moTEIs               ) ;
# endif
