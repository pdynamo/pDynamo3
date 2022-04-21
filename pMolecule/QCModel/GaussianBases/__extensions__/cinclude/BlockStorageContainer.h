# ifndef _BLOCKSTORAGECONTAINER
# define _BLOCKSTORAGECONTAINER

# include "BlockStorage.h"
# include "Boolean.h"
# include "Integer.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The container type. */
typedef struct {
    Boolean        isOwner  ;
    Integer        capacity ;
    BlockStorage **entries  ;
} BlockStorageContainer ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern BlockStorageContainer *BlockStorageContainer_Allocate   ( const Integer                 capacity ,
                                                                       Status                 *status   ) ;
extern void                   BlockStorageContainer_Deallocate (       BlockStorageContainer **self     ) ;

# endif
