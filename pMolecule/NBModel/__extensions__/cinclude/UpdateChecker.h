# ifndef _UPDATECHECKER
# define _UPDATECHECKER

# include "Boolean.h"
# include "Coordinates3.h"
# include "ImagePairListContainer.h"
# include "Real.h"
# include "Selection.h"
# include "SymmetryParameters.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The update checker type. */
typedef struct {
    Real buffer ;
} UpdateChecker ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern UpdateChecker *UpdateChecker_Allocate            (       Status                 *status              ) ;
extern Boolean        UpdateChecker_CheckForImageUpdate ( const SymmetryParameters     *set1                ,
                                                                SymmetryParameters     *set2                ,
                                                                ImagePairListContainer *images              ,
                                                          const Real                    buffer              ,
                                                          const Real                    maximumDisplacement ) ;
extern Boolean        UpdateChecker_CheckForUpdate      ( const Coordinates3           *set1                ,   
                                                                Coordinates3           *set2                ,   
                                                                Selection              *freeAtoms           ,   
                                                          const Real                    buffer              ,   
                                                                Real                   *maximumDisplacement ) ; 
extern void           UpdateChecker_Deallocate          (       UpdateChecker         **self                ) ;

# endif
