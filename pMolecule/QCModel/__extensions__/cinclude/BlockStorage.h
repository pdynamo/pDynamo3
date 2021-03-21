# ifndef _BLOCKSTORAGE
# define _BLOCKSTORAGE

# include "Boolean.h"
# include "Integer.h"
# include "List.h"
# include "MachineTypes.h"
# include "Real.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The block type. */
typedef struct {
   Integer     count     ;
   Cardinal16 *indices16 ;
   Cardinal32 *indices32 ;
   Real       *data      ;
} Block ;

/* . The block storage type. */
typedef struct {
    Boolean  checkUnderFlow ;
    Integer  blockSize      ;
    Integer  count          ;
    Integer  nIndices16     ;
    Integer  nIndices32     ;
    Integer  nReal          ;
    Real     underFlow      ;
    List    *blocks         ;
} BlockStorage ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Block. */
extern Block        *Block_Allocate          ( const Integer        blockSize         ,  
                                               const Integer        numberOfIndices16 ,  
                                               const Integer        numberOfIndices32 ,  
                                               const Integer        numberOfReal      ,  
                                                     Status        *status            ) ;
extern void          Block_Deallocate        (       Block        **self              ) ;

/* . Block storage. */
extern void          BlockStorage_AddData    (       BlockStorage  *self              ,
                                               const Integer        count             ,
                                               const Real          *data              ,
                                               const Cardinal16    *indices16         ,
                                               const Cardinal32    *indices32         ,  
                                                     Status        *status            ) ;
extern BlockStorage *BlockStorage_Allocate   (       Status        *status            ) ;
extern Integer       BlockStorage_Count      (       BlockStorage  *self              ) ;
extern void          BlockStorage_Deallocate (       BlockStorage **self              ) ;
extern void          BlockStorage_Empty      (       BlockStorage  *self              ) ;
extern Real          BlockStorage_ByteSize   (       BlockStorage  *self              ) ;
extern Block        *BlockStorage_Iterate    (       BlockStorage  *self              ) ;

# endif
