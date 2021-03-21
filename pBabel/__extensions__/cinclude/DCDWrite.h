# ifndef _DCDWRITE
# define _DCDWRITE

/* . Heavily modified from the VMD DCD plugin. */

# include "DCDHandle.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void      DCDWrite_Close  ( DCDHandle **self ) ;
extern DCDStatus DCDWrite_Frame  ( DCDHandle  *self ) ;
extern void      DCDWrite_Header ( DCDHandle  *self, const char *title ) ;
extern DCDStatus DCDWrite_Open   ( DCDHandle  *self, const char *path  ) ;

#endif
