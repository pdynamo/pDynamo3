# ifndef _MNDOCORECORE
# define _MNDOCORECORE

# include "Coordinates3.h"
# include "MNDOParametersContainer.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern Real MNDO_CoreCoreEnergy ( const MNDOParametersContainer *parameters   ,
                                  const Coordinates3            *coordinates3 ,
                                        Coordinates3            *gradients3   ) ;
# endif
