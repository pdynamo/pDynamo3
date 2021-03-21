# ifndef _MNDOPARAMETERSCONTAINER
# define _MNDOPARAMETERSCONTAINER

# include "Boolean.h"
# include "Integer.h"
# include "MNDOParameters.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The container type. */
typedef struct {
    Boolean          isOwner  ;
    Integer          capacity ;
    MNDOParameters **entries  ;
} MNDOParametersContainer ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern MNDOParametersContainer *MNDOParametersContainer_Allocate     ( const Integer                   capacity ,
                                                                             Status                   *status   ) ;
extern MNDOParametersContainer *MNDOParametersContainer_Clone        ( const MNDOParametersContainer  *self     ,
                                                                             Status                   *status   ) ;
extern void                     MNDOParametersContainer_Deallocate   (       MNDOParametersContainer **self     ) ;
extern Integer                  MNDOParametersContainer_LargestBasis ( const MNDOParametersContainer  *self     ) ;

# endif
