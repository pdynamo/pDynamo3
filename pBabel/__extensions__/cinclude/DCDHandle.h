# ifndef _DCDHANDLE
# define _DCDHANDLE

/* . Heavily modified from the VMD DCD plugin. */

/* . Types. */
# include <stdint.h>
# define Float32   float
# define Float64   double
# define Integer32 int32_t

# include "Boolean.h"
# include "Coordinates3.h"
# include "Integer.h"
# include "Real.h"
# include "SymmetryParameters.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Constants.*/
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661922
#endif

/* . File feature flags. */
#define DCD_IS_XPLOR        0x00
#define DCD_IS_CHARMM       0x01
#define DCD_HAS_4DIMS       0x02
#define DCD_HAS_EXTRA_BLOCK 0x04
#define DCD_HAS_64BIT_REC   0x08

/* . Status codes. */
typedef enum {
    DCDStatus_Normal                  =  0 ,
    DCDStatus_AtomNumberMismatch      =  1 ,
    DCDStatus_BadFormat               =  2 ,
    DCDStatus_BadRead                 =  3 ,
    DCDStatus_BadSeek                 =  4 ,
    DCDStatus_BadWrite                =  5 ,
    DCDStatus_FileAccessFailure       =  6 ,
    DCDStatus_InvalidDataObject       =  7 ,
    DCDStatus_InvalidFrameIndex       =  8 ,
    DCDStatus_OpenFailed              =  9 ,
    DCDStatus_OutOfMemory             = 10
} DCDStatus ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
typedef struct {
    Boolean             has4Dimensions      ;
    Boolean             hasCharges          ;
    Boolean             hasUnitCell         ;
    Boolean             isXPLOR             ;
    Boolean             reverseEndian       ;
    Boolean             useVelocityHeader   ;
    fio_fd              fileDescriptor      ;
    Integer             currentFrame        ;
    fio_size_t          fileSize            ;
    fio_size_t          firstFramePosition  ;
    fio_size_t          firstFrameSize      ;
    fio_size_t          frameSize           ;
    Integer             numberOfAtomIndices ;
    Integer             numberOfAtoms       ;
    Integer             numberOfFixedAtoms  ;
    Integer             numberOfFrames      ;
    Integer             recordMarkerScale   ;
    Integer             saveFrequency       ;
    Integer             startingFrame       ;
    Real                timeStep            ;
    Float64             unitCell[6]         ;
    Float32            *q                   ;
    Float32            *w                   ; 
    Float32            *x                   ;
    Float32            *y                   ;
    Float32            *z                   ;
    Integer32          *atomIndices         ; /* . Indexing starts at 1. */
    Coordinates3       *data3               ;
    SymmetryParameters *symmetryParameters  ;
/* . Eventually:
    char               *title               ;
    Integer             titleSize           ;
*/
} DCDHandle ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern DCDHandle *DCDHandle_Allocate              ( void ) ;
extern DCDStatus  DCDHandle_AllocateQW            ( DCDHandle  *self ) ;
extern DCDStatus  DCDHandle_CheckNumberOfAtoms    ( DCDHandle  *self, const Integer numberOfAtoms ) ;
extern Integer    DCDHandle_CurrentFrame          ( DCDHandle  *self ) ;
extern void       DCDHandle_Deallocate            ( DCDHandle **self ) ;
extern Integer    DCDHandle_NumberOfFrames        ( DCDHandle  *self ) ;
extern DCDStatus  DCDHandle_SetAtomIndices        ( DCDHandle  *self, Selection *selection ) ;
extern DCDStatus  DCDHandle_SetData3              ( DCDHandle  *self, Coordinates3 *data3 ) ;
extern DCDStatus  DCDHandle_SetSymmetryParameters ( DCDHandle  *self, SymmetryParameters *symmetryParameters ) ;

#endif
