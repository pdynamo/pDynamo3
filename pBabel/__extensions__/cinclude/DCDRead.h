# ifndef _DCDREAD
# define _DCDREAD

/* . Heavily modified from the VMD DCD plugin. */

# include "DCDHandle.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* READ macros with and without checks. */
# define CHECKREAD(fd, buf, size) { \
                                    auto fio_size_t n, s ; \
                                    s = ( size ) ; \
                                    n = fio_fread ( ( ( void * ) buf ), s, 1, ( fd ) ) ; \
                                    if ( n != 1 ) { (*status) = DCDStatus_BadRead ; goto FinishUp ; } \
                                  }
# define READ(fd, buf, size)  fio_fread ( ( ( void * ) buf ), ( size ), 1, ( fd ) )

/* . Scale factors. */
#define _RecordMarker32BitScale   1
#define _RecordMarker64BitScale   2
#define _RecordMarkerMaximumScale 2

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern void      DCDRead_Close     ( DCDHandle **self ) ;
extern DCDStatus DCDRead_Frame     ( DCDHandle  *self ) ;
extern DCDStatus DCDRead_GotoFrame ( DCDHandle  *self, const Integer f ) ;
extern DCDStatus DCDRead_Header    ( DCDHandle  *self ) ;
extern DCDStatus DCDRead_Open      ( DCDHandle  *self, const char *path ) ;
extern DCDStatus DCDRead_SkipFrame ( DCDHandle  *self ) ;

#endif
