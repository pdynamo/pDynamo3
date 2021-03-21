/*==================================================================================================================================
! . DCD trajectory file reading.
! . Heavily modified from the VMD DCD plugin.
!=================================================================================================================================*/

# include "fastio.h"

# include <stdio.h>
# include <sys/stat.h>
# include <sys/types.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>

# include "endianswap.h"
# include "DCDRead.h"
# include "Memory.h"

/*

The DCD format is a Fortran binary format. Each Fortran record consists of the data that is written sandwiched between
two record markers that contain the number of bytes in the record. The record marker sizes are system dependent.

The header consists of:
HDR,ICNTRL(20)             : CHAR*4, INTEGER 
NTITLE, TITLE(NTITLE)      : INTEGER*4, CHAR*80         
NATOMS                     : INTEGER         
ATOMINDICES(NATOMINDICES)  : INTEGER

The first frame consists of:
UNITCELL(6)                : REAL*8
XYZ/W/Q(NATOMS)            : REAL (W and Q optional)

The subsequent frames consist of:
UNITCELL(6)                : REAL*8
XYZ/W/Q(N)                 : REAL where N is NATOMINDICES if there are ATOMINDICES, otherwise NATOMS.

Normally INTEGER refers to the default integer type and REAL to REAL*4. However, it appears that on many 64 bit
machines, CHARMM uses an I4BINARY flag which forces the use of INTEGER*4. Likewise, there are some options that
permit REAL to be REAL*8 although I am uncertain how widely these are used. Note, however, that the record
markers will always be system dependent.

Thus, here we assume INTEGER*4 and REAL*4 with system-dependent record markers of either 32 or 64 bits.

Endian reversal is done in terms of words (either 4 or 8 byte). Strings are not affected.

*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void      DCDRead_CheckFrameCount           ( DCDHandle *self ) ;
static Boolean   DCDRead_DetermineTrajectoryFormat ( DCDHandle *self, DCDStatus *status ) ;
static Boolean   DCDRead_ProcessAtomIndices        ( DCDHandle *self, DCDStatus *status ) ;
static Boolean   DCDRead_ProcessControlFlags       ( DCDHandle *self, DCDStatus *status ) ;
static Boolean   DCDRead_ProcessNumberOfAtoms      ( DCDHandle *self, DCDStatus *status ) ;
static Boolean   DCDRead_ProcessTitle              ( DCDHandle *self, DCDStatus *status ) ;
static DCDStatus DCDRead_Step                      ( fio_fd      fd                ,
                                                     Boolean     has4Dimensions    ,
                                                     Boolean     hasCharges        ,
                                                     Boolean     hasUnitCell       ,
                                                     Boolean     reverseEndian     ,
                                                     Integer     numberOfAtoms     ,
                                                     Integer     numberOfCharges   ,
                                                     Integer     recordMarkerScale ,
                                                     Float32    *q                 ,
                                                     Float32    *w                 ,
                                                     Float32    *x                 ,
                                                     Float32    *y                 ,
                                                     Float32    *z                 ,
                                                     Float64    *u                 ) ;
static Boolean   Record_Read                       ( fio_fd      fd                ,
                                                     void       *buffer            ,
                                                     fio_size_t  numberOfBytes     ,
                                                     Integer     recordMarkerScale ,
                                                     Boolean     reverseEndian     ,
                                                     DCDStatus  *status            ) ;
static Boolean   Record_Skip                       ( fio_fd      fd                ,
                                                     fio_size_t *numberOfBytes     ,
                                                     Integer     recordMarkerScale ,
                                                     Boolean     reverseEndian     ,
                                                     DCDStatus  *status            ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Close after reading.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DCDRead_Close ( DCDHandle **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        fio_fclose ( (*self)->fileDescriptor ) ;
        DCDHandle_Deallocate ( self ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Read a frame.
!---------------------------------------------------------------------------------------------------------------------------------*/
DCDStatus DCDRead_Frame ( DCDHandle *self )
{
    DCDStatus status = DCDStatus_Normal ;
    if ( self != NULL )
    {
        Boolean reducedData   ;
        Integer numberOfAtoms ;
        /* . Read the data. */
        reducedData   = ( ( self->currentFrame > 0 ) && ( self->atomIndices != NULL ) && ( self->numberOfAtomIndices > 0 ) ) ;
        numberOfAtoms = reducedData ? self->numberOfAtomIndices : self->numberOfAtoms ;
        self->currentFrame++ ;
        status = DCDRead_Step ( self->fileDescriptor    ,
                                self->has4Dimensions    ,
                                self->hasCharges        ,
                                self->hasUnitCell       ,
                                self->reverseEndian     ,
                                numberOfAtoms           ,
                                self->numberOfAtoms     ,
                                self->recordMarkerScale ,
                                self->q                 ,
                                self->w                 ,
                                self->x                 ,
                                self->y                 ,
                                self->z                 ,
                                self->unitCell          ) ;
        /* . Process the data. */
        if ( status == DCDStatus_Normal )
        {
            auto Integer             i, i3 ;
            auto Real               *data = self->data3->data ;
            auto SymmetryParameters *p = self->symmetryParameters ;
            /* . Unit cell. */
            if ( ( self->hasUnitCell ) && ( p != NULL ) )
            {
                auto Float64 *u = self->unitCell ;
                /* . Cosines of the angles. */
                if ( ( u[1] >= -1.0 ) && ( u[1] <= 1.0 ) &&
                     ( u[3] >= -1.0 ) && ( u[3] <= 1.0 ) &&
                     ( u[4] >= -1.0 ) && ( u[4] <= 1.0 ) )
                {
                    u[4] = 90.0 - asin ( u[4] ) * 90.0 / M_PI_2 ;
                    u[3] = 90.0 - asin ( u[3] ) * 90.0 / M_PI_2 ;
                    u[1] = 90.0 - asin ( u[1] ) * 90.0 / M_PI_2 ;
                }
                SymmetryParameters_SetCrystalParameters ( p, u[0], u[2], u[5], u[4], u[3], u[1] ) ;
            }
            /* . Reduced data. */
            if ( reducedData )
            {
                for ( i = 0 ; i < numberOfAtoms ; i++ )
                {
                    i3 = 3 * ( self->atomIndices[i] - 1 ) ;
                   data[i3  ] =  self->x[i] ;
                   data[i3+1] =  self->y[i] ;
                   data[i3+2] =  self->z[i] ;
                }
            }
            /* . All data. */
            else
            {
                for ( i = i3 = 0 ; i < numberOfAtoms ; i++, i3 += 3 )
                {
                    data[i3  ] = self->x[i] ;
                    data[i3+1] = self->y[i] ;
                    data[i3+2] = self->z[i] ;
                }
            }
        }
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Go to the beginning of a frame.
!---------------------------------------------------------------------------------------------------------------------------------*/
DCDStatus DCDRead_GotoFrame ( DCDHandle *self, const Integer f )
{
    DCDStatus status = DCDStatus_Normal ;
    if ( self != NULL )
    {
        if ( ( f < 0 ) || ( f >= self->numberOfFrames ) ) { status = DCDStatus_InvalidFrameIndex ; }
        else
        {
            fio_size_t p = self->firstFramePosition ;
            if ( f > 0 ) p += self->firstFrameSize + ( f - 1 ) * self->frameSize ;
            if ( fio_fseek ( self->fileDescriptor, p, FIO_SEEK_SET ) != 0 ) status = DCDStatus_BadSeek ;
            self->currentFrame = f ;
        }
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Read header.
!---------------------------------------------------------------------------------------------------------------------------------*/
DCDStatus DCDRead_Header ( DCDHandle *self )
{
    DCDStatus status = DCDStatus_Normal ;
    if ( self != NULL )
    {
        /* . Process the various records in the header. */
        if ( ! DCDRead_DetermineTrajectoryFormat ( self, &status ) ) goto FinishUp ;
        if ( ! DCDRead_ProcessControlFlags       ( self, &status ) ) goto FinishUp ;
        if ( ! DCDRead_ProcessTitle              ( self, &status ) ) goto FinishUp ;
        if ( ! DCDRead_ProcessNumberOfAtoms      ( self, &status ) ) goto FinishUp ;
        if ( ! DCDRead_ProcessAtomIndices        ( self, &status ) ) goto FinishUp ;

        /* . Find the current position. */
        self->firstFramePosition = fio_ftell ( self->fileDescriptor ) ;
        if ( self->firstFramePosition < 0 ) status = DCDStatus_BadSeek ;

        /* . Check the frame count. */
        DCDRead_CheckFrameCount ( self ) ;
    }
FinishUp:
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Open for reading.
!---------------------------------------------------------------------------------------------------------------------------------*/
DCDStatus DCDRead_Open ( DCDHandle *self, const char *path )
{
    DCDStatus status = DCDStatus_Normal ;
    if ( ( self != NULL ) && ( path != NULL ) )
    {
        fio_fd fileDescriptor ;
        struct stat stbuf ;
        /* . See if the file exists and get its size. */
        memset ( &stbuf, 0, sizeof ( struct stat ) ) ;
        if ( stat ( path, &stbuf ) ) { status = DCDStatus_FileAccessFailure ; }
        /* . OK. */
        else
        {
            /* . Open the file and get its size. */
            if ( fio_open ( path, FIO_READ, &fileDescriptor ) < 0 ) { status = DCDStatus_OpenFailed ; }
            else
            {
                self->fileDescriptor = fileDescriptor ;
#if defined(_MSC_VER) && defined(FASTIO_NATIVEWIN32)
                /*
                ! . The stat() call is not 64-bit savvy on Windows so we have to use the fastio fseek/ftell 
                ! . routines for this until we add a portable filesize routine for this purpose.
                */
                {
                    auto fio_size_t currentPosition = fio_ftell ( fileDescriptor ) ;
                    fio_fseek ( fileDescriptor , 0, FIO_SEEK_END ) ; /* . Seek to the end of the file. */
                    self->fileSize = fio_ftell ( fileDescriptor ) ;
                    fio_fseek ( fileDescriptor, currentPosition, FIO_SEEK_SET ) ;  /* . Return to the current position. */
                }
#else
                self->fileSize = stbuf.st_size ; /* . This works OK on Unix machines. */
#endif
            }
        }
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Skip a frame.
!---------------------------------------------------------------------------------------------------------------------------------*/
DCDStatus DCDRead_SkipFrame ( DCDHandle *self )
{
    DCDStatus status = DCDStatus_Normal ;
    if ( self != NULL )
    {
        fio_size_t offset ;
        if ( self->currentFrame == 0 ) offset = self->firstFrameSize ;
        else                           offset = self->frameSize      ;
        if ( fio_fseek ( self->fileDescriptor, offset, FIO_SEEK_CUR ) != 0 ) status = DCDStatus_BadSeek ;
        self->currentFrame += 1 ;
    }
    return status ;
}

/*==================================================================================================================================
! . Local DCD-specific procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Check the frame count.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void DCDRead_CheckFrameCount ( DCDHandle *self )
{
    fio_size_t d = 3, fullRecord, markers, reducedRecord, trajectorySize, unitCellRecord ;
    /* . Basic record sizes. */
    markers = 2 * self->recordMarkerScale * sizeof ( Integer32 ) ;
                                         fullRecord     = self->numberOfAtoms       * sizeof ( Float32 ) + markers ;
    if ( self->numberOfAtomIndices > 0 ) reducedRecord  = self->numberOfAtomIndices * sizeof ( Float32 ) + markers ;
    else                                 reducedRecord  = fullRecord ;
                                         unitCellRecord = self->hasUnitCell ? 6 * sizeof ( Float64 ) + markers : 0 ;
    /* . Frame sizes. */
    if ( self->has4Dimensions ) d += 1 ;
    self->firstFrameSize =    fullRecord * d + unitCellRecord ;
    self->     frameSize = reducedRecord * d + unitCellRecord ;
    if ( self->hasCharges )
    {
        self->firstFrameSize += fullRecord ;
        self->     frameSize += fullRecord ;
    }
    /* . Determine the number of frames on the trajectory. */
    trajectorySize = self->fileSize - self->firstFramePosition - self->firstFrameSize ;
    if ( trajectorySize < 0 ) { self->numberOfFrames = 0 ; }
    else                      { self->numberOfFrames = ( trajectorySize / self->frameSize ) + 1 ; }
    self->currentFrame = 0 ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Determine the format of the trajectory from the header.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Boolean DCDRead_DetermineTrajectoryFormat ( DCDHandle *self, DCDStatus *status )
{
    auto Integer32 recordMarker[2] ;

    /* . The length of the first record should be 84 (CHAR*4 + 20 * INTEGER*4). */
    CHECKREAD ( self->fileDescriptor, recordMarker, 2 * sizeof ( Integer32 ) ) ;

    /* . Check the magic number in the file header and determine the byte order. */
    if ( ( recordMarker[0] + recordMarker[1] ) == 84 )
    {
        self->reverseEndian     = False             ;
        self->recordMarkerScale = _RecordMarker64BitScale ;
    }
    else if ( recordMarker[0] == 84 )
    {
        self->reverseEndian     = False             ;
        self->recordMarkerScale = _RecordMarker32BitScale ;
    }
    else
    {
        swap4_aligned ( recordMarker, 2 ) ;
        if ( ( recordMarker[0] + recordMarker[1] ) == 84 )
        {
           self->reverseEndian     = True              ;
           self->recordMarkerScale = _RecordMarker64BitScale ;
        }
        else if ( recordMarker[0] == 84 )
        {
            self->reverseEndian     = True              ;
            self->recordMarkerScale = _RecordMarker32BitScale ;
        }
        else { (*status) = DCDStatus_BadFormat ; goto FinishUp ; }
    }

    /* . Return to the beginning of the file. */
    if ( fio_fseek ( self->fileDescriptor, 0, FIO_SEEK_SET ) != 0 ) (*status) = DCDStatus_BadSeek ;
FinishUp:
    return ( (*status) == DCDStatus_Normal ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Process the atom indices in the header.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Boolean DCDRead_ProcessAtomIndices ( DCDHandle *self, DCDStatus *status )
{
    self->numberOfAtomIndices = ( self->numberOfFixedAtoms > 0 ) ? ( self->numberOfAtoms - self->numberOfFixedAtoms ) : 0 ;
    if ( self->numberOfAtomIndices > 0 )
    {
        /* . Allocate space. */
        self->atomIndices = Memory_AllocateArrayOfTypes ( self->numberOfAtomIndices, Integer32 ) ;
        if ( self->atomIndices == NULL ) { (*status) = DCDStatus_OutOfMemory ; }
        else
        {
            if ( Record_Read ( self->fileDescriptor      ,
                               self->atomIndices         ,
                               self->numberOfAtomIndices * sizeof ( Integer32 ) ,
                               self->recordMarkerScale   ,
                               self->reverseEndian       ,
                               status                    ) )
            {
                if ( self->reverseEndian ) swap4_aligned ( self->atomIndices, self->numberOfAtomIndices ) ;
            }
        }
    }
    return ( (*status) == DCDStatus_Normal ) ; 
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Process the control flags in the header.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . 84 is 4 (the header) and 20 * 4 (control flags).
! . Note that the control flags start at flags[1] due to the presence of the header.
*/
static Boolean DCDRead_ProcessControlFlags ( DCDHandle *self, DCDStatus *status )
{
    Integer32 flags[21] ;
    char *headerBuffer = ( char * ) flags ;
    if ( Record_Read ( self->fileDescriptor, headerBuffer, 84, self->recordMarkerScale, self->reverseEndian, status ) )
    {
        auto Integer32 dcdCORDMagic, dcdVELDMagic ;
        auto char *cord = ( char * ) &dcdCORDMagic, *veld = ( char * ) &dcdVELDMagic ; 

        /* . DCD file magic strings "CORD" and "VELD". */
        cord[0] = 'C' ; cord[1] = 'O' ; cord[2] = 'R' ; cord[3] = 'D' ;
        veld[0] = 'V' ; veld[1] = 'E' ; veld[2] = 'L' ; veld[3] = 'D' ;

        /* . Check the magic string. */
             if ( flags[0] == dcdVELDMagic ) { self->useVelocityHeader = True ; }
        else if ( flags[0] != dcdCORDMagic ) { (*status) = DCDStatus_BadFormat ; goto FinishUp ; }

        /* . Check for a CHARMM or XPLOR file. */
        /* . XPLOR file has no version number - assumed to be zero so endianness irrelevant. */
        self->isXPLOR = ( flags[20] == 0 ) ;
        if ( self->isXPLOR )
        {
            self->timeStep = * ( ( Float64 * ) ( headerBuffer + 40 ) ) ;
            if ( self->reverseEndian ) swap8_unaligned ( &(self->timeStep), 1 ) ;
        }
        else
        {
            auto Float32 t ;
            t = * ( ( Float32 * ) ( headerBuffer + 40 ) ) ;
            if ( self->reverseEndian ) swap4_aligned ( &t, 1 ) ;
            self->timeStep = ( Float64 ) t ;
            self->hasUnitCell    = flags[11] ;
            self->has4Dimensions = flags[12] ;
            self->hasCharges     = flags[13] ;
        }

        /* . Other counters. */
        self->numberOfFrames     = flags[1] ;
        self->startingFrame      = flags[2] ;
        self->saveFrequency      = flags[3] ;
        self->numberOfFixedAtoms = flags[9] ;
        if ( self->reverseEndian )
        {
            swap4_unaligned ( &(self->numberOfFrames    ) , 1 ) ;
            swap4_unaligned ( &(self->startingFrame     ) , 1 ) ;
            swap4_unaligned ( &(self->saveFrequency     ) , 1 ) ;
            swap4_unaligned ( &(self->numberOfFixedAtoms) , 1 ) ;
        }
    }
FinishUp:
    return ( (*status) == DCDStatus_Normal ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Process the number of atoms in the header.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Boolean DCDRead_ProcessNumberOfAtoms ( DCDHandle *self, DCDStatus *status )
{
    Integer32 n ;
    if ( Record_Read ( self->fileDescriptor, &n, sizeof ( Integer32 ), self->recordMarkerScale, self->reverseEndian, status ) )
    {
        if ( self->reverseEndian ) swap4_aligned ( &n, 1 ) ;
        self->numberOfAtoms = ( Integer ) n ;
    }
    return ( (*status) == DCDStatus_Normal ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Process the title in the header (not saved).
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . NTITLE (INTEGER*4) + TITLE*80(1:NTITLE).
! . The record is currently skipped.
! . It would be straightforward to read the raw record data and save the data.
! . A useful check on the number of bytes read would be:
! . ( ( ( ( recordMarker[0] + recordMarker[1] ) - sizeof ( Integer32 ) ) % 80 ) == 0 )
*/
static Boolean DCDRead_ProcessTitle ( DCDHandle *self, DCDStatus *status )
{
    return Record_Skip ( self->fileDescriptor    ,
                         NULL                    ,
                         self->recordMarkerScale ,
                         self->reverseEndian     ,
                         status                  ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Read a step.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _MaximumItems                 13
# define _MaximumNumberOfRecordMarkers 12
static DCDStatus DCDRead_Step ( fio_fd   fd                ,
                                Boolean  has4Dimensions    ,
                                Boolean  hasCharges        ,
                                Boolean  hasUnitCell       ,
                                Boolean  reverseEndian     ,
                                Integer  numberOfAtoms     ,
                                Integer  numberOfCharges   ,
                                Integer  recordMarkerScale ,
                                Float32 *q                 ,
                                Float32 *w                 ,
                                Float32 *x                 ,
                                Float32 *y                 ,
                                Float32 *z                 ,
                                Float64 *u                 )
{
    DCDStatus  status = DCDStatus_Normal ;
    fio_iovec  iov[_MaximumItems] ;
    fio_size_t numberOfBytes = 0, numberOfBytesRead ;
    Integer32  b = 0, buffer[_MaximumNumberOfRecordMarkers*_RecordMarkerMaximumScale], i, n = 0 ;
    /* . Set up the I/O data. */
    /* . First record marker. */
    iov[n].iov_base = ( fio_caddr_t ) &buffer[b] ;
    iov[n].iov_len  = recordMarkerScale * sizeof ( Integer32 ) ;
    b += recordMarkerScale ; n += 1 ;
    /* . Unit cell and two record markers (stop/start). */
    if ( hasUnitCell )
    {
        iov[n].iov_base = ( fio_caddr_t ) u ;
        iov[n].iov_len  = 6 * sizeof ( Float64 ) ;
        n += 1 ;
        iov[n].iov_base = ( fio_caddr_t ) &buffer[b];
        iov[n].iov_len  = 2 * recordMarkerScale * sizeof ( Integer32 ) ;
        b += ( 2 * recordMarkerScale ) ; n += 1 ;
    }
    /* . X, Y, Z and optionally W and Q. */
    iov[n].iov_base = ( fio_caddr_t ) x ;
    iov[n].iov_len  = numberOfAtoms * sizeof ( Float32 ) ;
    n += 1 ;
    iov[n].iov_base = ( fio_caddr_t ) &buffer[b] ;
    iov[n].iov_len  = 2 * recordMarkerScale * sizeof ( Integer32 ) ;
    b += ( 2 * recordMarkerScale ) ; n += 1 ;
    iov[n].iov_base = ( fio_caddr_t ) y ;
    iov[n].iov_len  = numberOfAtoms * sizeof ( Float32 ) ;
    n += 1 ;
    iov[n].iov_base = ( fio_caddr_t ) &buffer[b] ;
    iov[n].iov_len  = 2 * recordMarkerScale * sizeof ( Integer32 ) ;
    b += ( 2 * recordMarkerScale ) ; n += 1 ;
    iov[n].iov_base = ( fio_caddr_t ) z ;
    iov[n].iov_len  = numberOfAtoms * sizeof ( Float32 ) ;
    n += 1 ;
    if ( has4Dimensions )
    {
        iov[n].iov_base = ( fio_caddr_t ) &buffer[b] ;
        iov[n].iov_len  = 2 * recordMarkerScale * sizeof ( Integer32 ) ;
        b += ( 2 * recordMarkerScale ) ; n += 1 ;
        iov[n].iov_base = ( fio_caddr_t ) w ;
        iov[n].iov_len  = numberOfAtoms * sizeof ( Float32 ) ;
        n += 1 ;
    }
    if ( hasCharges )
    {
        iov[n].iov_base = ( fio_caddr_t ) &buffer[b] ;
        iov[n].iov_len  = 2 * recordMarkerScale * sizeof ( Integer32 ) ;
        b += ( 2 * recordMarkerScale ) ; n += 1 ;
        iov[n].iov_base = ( fio_caddr_t ) q ;
        iov[n].iov_len  = numberOfCharges * sizeof ( Float32 ) ;
        n += 1 ;
    }
    iov[n].iov_base = ( fio_caddr_t ) &buffer[b] ;
    iov[n].iov_len  = recordMarkerScale * sizeof ( Integer32 ) ;
    b += recordMarkerScale ; n += 1 ;
    for ( i = 0 ; i < n ; i++ ) numberOfBytes += iov[i].iov_len ;
    /* . Read the data. */
    numberOfBytesRead = fio_readv ( fd, &iov[0], n ) ;
    /* . Everything OK so far. */
    if ( numberOfBytesRead == numberOfBytes )
    {
        auto Boolean isOK  = True ;
        auto Integer e, n, nQ, s = 0 ;
        /* . Endianism conversion. */
        if ( reverseEndian )
        {
            swap4_aligned ( &buffer[0] , b ) ;
            swap4_aligned ( x, numberOfAtoms ) ;
            swap4_aligned ( y, numberOfAtoms ) ;
            swap4_aligned ( z, numberOfAtoms ) ;
            if ( has4Dimensions ) swap4_aligned ( w, numberOfAtoms   ) ;
            if ( hasCharges     ) swap4_aligned ( q, numberOfCharges ) ;
            if ( hasUnitCell    ) swap8_aligned ( u, 6 ) ;
        }
        /* . Check the Fortran format size values for safety. */
        e  = hasCharges ? b - 2 : b ;
        n  = numberOfAtoms   * sizeof ( Float32 ) ;
        nQ = numberOfCharges * sizeof ( Float32 ) ;
        if ( recordMarkerScale == 1 )
        {
            if ( hasUnitCell ) { isOK = isOK && ( buffer[0] == 48 ) && ( buffer[1] == 48 ) ;  s+= 2 ; }
            for ( i = s ; i < e ; i++ ) isOK = isOK && ( buffer[i] == n  ) ;
            for ( i = e ; i < b ; i++ ) isOK = isOK && ( buffer[i] == nQ ) ;
        }
        else
        {
            if ( hasUnitCell ) { isOK = isOK && ( ( ( buffer[0] + buffer[1] ) == 48 ) && ( ( buffer[2] + buffer[3] ) == 48 ) ) ;  s+= 2 ; }
            for ( i = s ; i < e ; i++ ) isOK = isOK && ( ( buffer[2*i] + buffer[2*i+1] ) == n  ) ;
            for ( i = e ; i < b ; i++ ) isOK = isOK && ( ( buffer[2*i] + buffer[2*i+1] ) == nQ ) ;
        }
        if ( ! isOK ) status = DCDStatus_BadFormat ;
    }
    /* . Problem reading. */
    else { status = DCDStatus_BadRead ; }
    return status ;
}

/*==================================================================================================================================
! . Local record procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Read a record of known size. No changes are made to the buffer data.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Boolean Record_Read ( fio_fd      fd                ,
                             void       *buffer            ,
                             fio_size_t  numberOfBytes     ,
                             Integer     recordMarkerScale ,
                             Boolean     reverseEndian     ,
                             DCDStatus  *status            )
{
    Integer32 recordMarker[2] ;
    recordMarker[1] = 0 ;
    CHECKREAD ( fd, recordMarker, recordMarkerScale * sizeof ( Integer32 ) ) ;
    if ( reverseEndian ) swap4_aligned ( recordMarker, recordMarkerScale ) ;
    if ( ( recordMarker[0] + recordMarker[1] ) != numberOfBytes ) { (*status) = DCDStatus_BadFormat ; goto FinishUp ; }
    CHECKREAD ( fd, buffer, ( numberOfBytes ) ) ;
    recordMarker[1] = 0 ;
    CHECKREAD ( fd, recordMarker, recordMarkerScale * sizeof ( Integer32 ) ) ;
    if ( reverseEndian ) swap4_aligned ( recordMarker, recordMarkerScale ) ;
    if ( ( recordMarker[0] + recordMarker[1] ) != numberOfBytes ) { (*status) = DCDStatus_BadFormat ; goto FinishUp ; }
FinishUp:
    return ( (*status) == DCDStatus_Normal ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Skip a record. No knowledge of the record size is required.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Boolean Record_Skip ( fio_fd      fd                ,
                             fio_size_t *numberOfBytes     ,
                             Integer     recordMarkerScale ,
                             Boolean     reverseEndian     ,
                             DCDStatus  *status            )
{
    fio_size_t size = 0 ;
    Integer32  recordMarker[2] ;
    recordMarker[1] = 0 ;
    CHECKREAD ( fd, recordMarker, recordMarkerScale * sizeof ( Integer32 ) ) ;
    if ( reverseEndian ) swap4_aligned ( recordMarker, recordMarkerScale ) ;
    size = ( recordMarker[0] + recordMarker[1] ) ;
    if ( fio_fseek ( fd, size, FIO_SEEK_CUR ) != 0 ) { (*status) = DCDStatus_BadSeek ; goto FinishUp ; }
    recordMarker[1] = 0 ;
    CHECKREAD ( fd, recordMarker, recordMarkerScale * sizeof ( Integer32 ) ) ;
    if ( reverseEndian ) swap4_aligned ( recordMarker, recordMarkerScale ) ;
    if ( ( recordMarker[0] + recordMarker[1] ) != size ) { (*status) = DCDStatus_BadFormat ; goto FinishUp ; }
FinishUp:
    if ( numberOfBytes != NULL ) (*numberOfBytes) = size ;
    return ( (*status) == DCDStatus_Normal ) ;
}
