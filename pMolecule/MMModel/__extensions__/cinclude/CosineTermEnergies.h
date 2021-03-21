# ifndef _COSINETERMENERGIES
# define _COSINETERMENERGIES

# include "Coordinates3.h"
# include "CosineTermContainer.h"
# include "Real.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Real CosineTermEnergy_Angle      ( const CosineTermContainer *self         ,
                                          const Coordinates3        *coordinates3 ,
                                                Coordinates3        *gradients3   ) ;
extern Real CosineTermEnergy_Dihedral   ( const CosineTermContainer *self         ,
                                          const Coordinates3        *coordinates3 ,
                                                Coordinates3        *gradients3   ) ;
extern Real CosineTermEnergy_OutOfPlane ( const CosineTermContainer *self         ,
                                          const Coordinates3        *coordinates3 ,
                                                Coordinates3        *gradients3   ) ;

# endif
