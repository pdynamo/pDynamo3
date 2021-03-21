/*==================================================================================================================================
! . DCD trajectory file writing.
! . Heavily modified from the VMD DCD plugin.
!=================================================================================================================================*/

# include "fastio.h"

# include <stdio.h>
# include <sys/stat.h>
# include <sys/types.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>
# include <time.h>

# include "DCDWrite.h"
# include "Memory.h"

/* . All I/O, including record markers, is done with 32 bit words except for unit cell which requires 64 bit words. */

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _CharmmVersion 34 /* . Pretend to be CHARMM version 34. */

/* . Defines used by DCDWrite_Step. */
# define NFILE_POS  8L
# define NSTEP_POS 20L

/* . WRITE macro. */
# define WRITE(fd, buf, size) fio_fwrite ( ( ( void * ) buf ), ( size ), 1, ( fd ) )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static DCDStatus DCDWrite_Step (       fio_fd   fd            ,
                                       Boolean  hasUnitCell   ,
                                       Integer  currentFrame  ,
                                       Integer  currentStep   ,
                                       Integer  numberOfAtoms ,
                                 const Float32 *x             ,
                                 const Float32 *y             ,
                                 const Float32 *z             ,
                                 const Float64 *unitCell      ) ;

static int FastIO_WriteInteger32 ( fio_fd fd, Integer32 i ) { return ( fio_fwrite ( &i, 4, 1, fd ) != 1 ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Close after writing.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DCDWrite_Close ( DCDHandle **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        fio_fclose ( (*self)->fileDescriptor ) ;
        DCDHandle_Deallocate ( self ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Write a frame.
!---------------------------------------------------------------------------------------------------------------------------------*/
DCDStatus DCDWrite_Frame ( DCDHandle *self )
{
    DCDStatus status = DCDStatus_Normal ;
    if ( self != NULL )
    {
        auto Integer currentStep, i, i3, numberOfAtoms ;
        auto Real   *data = self->data3->data ;
        /* . Counters. */
        self->numberOfFrames++ ;
        currentStep = self->startingFrame + self->numberOfFrames * self->saveFrequency ;
        /* . All data. */
        if ( ( self->numberOfFrames == 1 ) || ( self->atomIndices == NULL ) || ( self->numberOfAtomIndices <= 0 ) )
        {
            numberOfAtoms = self->numberOfAtoms ;
            for ( i = i3 = 0 ; i < numberOfAtoms ; i++, i3 += 3 )
            {
                self->x[i] = data[i3  ] ;
                self->y[i] = data[i3+1] ;
                self->z[i] = data[i3+2] ;
            }
        }
        /* . Reduced data. */
        else
        {
            numberOfAtoms = self->numberOfAtomIndices ;
            for ( i = 0 ; i < numberOfAtoms ; i++ )
            {
                i3 = 3 * ( self->atomIndices[i] - 1 ) ;
                self->x[i] = data[i3  ] ;
                self->y[i] = data[i3+1] ;
                self->z[i] = data[i3+2] ;
            }
        }
        /* . Unit cell. */
        if ( self->hasUnitCell )
        {
            auto Float64            *u = self->unitCell           ;
            auto SymmetryParameters *p = self->symmetryParameters ;
            u[0] = p->a;
            u[2] = p->b;
            u[5] = p->c;
            u[1] = sin ( ( M_PI_2 / 90.0 ) * ( 90.0 - p->gamma ) ) ;
            u[3] = sin ( ( M_PI_2 / 90.0 ) * ( 90.0 - p->beta  ) ) ;
            u[4] = sin ( ( M_PI_2 / 90.0 ) * ( 90.0 - p->alpha ) ) ;
        }
        /* . Writing. */
        status = DCDWrite_Step ( self->fileDescriptor ,
                                 self->hasUnitCell    ,
                                 self->numberOfFrames ,
                                 currentStep          ,
                                 numberOfAtoms        ,
                                 self->x              ,
                                 self->y              ,
                                 self->z              ,
                                 self->unitCell       ) ;
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Write header.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . CHARMM control flags:
!
!  0 NFILE           : number of frames in the file.
!  1 NPRIV           : starting timestep (ISTART).
!  2 NSAVC           : timesteps between the frames written to the file.
!  3 NSTEP           : number of timesteps in the simulation.
!  4 NSAVC           : for velocity trajectory.
!  7 NDEGF           : number of degrees of freedom.
!  8 NATOM - LNFREAT : number of fixed atoms.
!  9 DELTA4          : the timestep (single precision).
! 10 QCRYS           : is there a unitcell?
! 11 DIM4            : is there a fourth dimension?
! 12 QCG             : are there charges?
! 19 VERNUM          : the version number.
! 
! . XPLOR modifications:
!
!  9/10 DELTA8       : the timestep (double precision).
! 19                 : no version number.
!
! . It should be possible to handle the title better here (i.e. divide into several 80 character chunks).
*/
void DCDWrite_Header ( DCDHandle *self, const char *title )
{
    if ( self != NULL )
    {
        auto fio_fd     fd = self->fileDescriptor ;
        auto Float32    outFloat32 ;
        auto Float64    outFloat64 ;
        auto Integer    i ;
        auto Integer32  controlFlags[23], numberOfBytes ;
        auto char       timeString[81], titleString[200] ;
        auto time_t     currentTime  ;
        auto struct tm *timeBuffer   ;
        /* . Header. */
        numberOfBytes = 84 ;
        FastIO_WriteInteger32 ( fd, numberOfBytes ) ;
        if ( self->useVelocityHeader ) strcpy ( titleString, "VELD" ) ;
        else                           strcpy ( titleString, "CORD" ) ;
        numberOfBytes =  4 ;
        WRITE ( fd, titleString, numberOfBytes ) ;
        /* . Control flags with trailing record counters. */
        for ( i = 0 ; i < 20 ; i++ ) controlFlags[i] = 0 ;
        controlFlags[1] = self->startingFrame ;
        controlFlags[2] = self->saveFrequency ;
        controlFlags[8] = ( self->atomIndices == NULL ) ? 0 : self->numberOfAtoms - self->numberOfAtomIndices ;
        if ( self->isXPLOR )
        {
            outFloat64 = self->timeStep ; Memory_CopyTo ( &outFloat64, &controlFlags[9], sizeof ( Float64 ) ) ;
        }
        else
        {
            outFloat32 = self->timeStep ; Memory_CopyTo ( &outFloat32, &controlFlags[9], sizeof ( Float32 ) ) ;
            controlFlags[10] = self->hasUnitCell ;
            controlFlags[19] = _CharmmVersion    ;
        }
        controlFlags[20] =  84 ;
        controlFlags[21] = 164 ;
        controlFlags[22] =   2 ;
        numberOfBytes    =  92 ;
        WRITE ( fd, controlFlags, numberOfBytes ) ;
        /* . Title and time string. */
        strncpy ( titleString, title, 80 ) ;
        titleString[79] = '\0';
        numberOfBytes   =  80 ;
        WRITE ( fd, titleString, numberOfBytes ) ;
        currentTime = time      ( NULL         ) ;
        timeBuffer  = localtime ( &currentTime ) ;
        strftime ( timeString, 80, "REMARKS Created %d %B, %Y at %R", timeBuffer ) ;
        WRITE    ( fd, timeString, numberOfBytes ) ;
        FastIO_WriteInteger32 ( fd, 164 ) ;
        /* . Number of atoms. */
        FastIO_WriteInteger32 ( fd, 4 ) ;
        FastIO_WriteInteger32 ( fd, self->numberOfAtoms ) ;
        FastIO_WriteInteger32 ( fd, 4 ) ;
        /* . Atom indices. */
        if ( self->numberOfAtomIndices > 0 )
        {
           numberOfBytes = self->numberOfAtomIndices * sizeof ( Integer32 ) ;
           FastIO_WriteInteger32 ( fd, numberOfBytes ) ;
           WRITE                 ( fd, self->atomIndices, numberOfBytes ) ;
           FastIO_WriteInteger32 ( fd, numberOfBytes ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Open for writing.
!---------------------------------------------------------------------------------------------------------------------------------*/
DCDStatus DCDWrite_Open ( DCDHandle *self, const char *path )
{
    DCDStatus status = DCDStatus_Normal ;
    if ( ( self != NULL ) && ( path != NULL ) )
    {
        fio_fd fileDescriptor ;
        if ( fio_open ( path, FIO_WRITE, &fileDescriptor ) < 0 ) { status = DCDStatus_OpenFailed ; }
        else self->fileDescriptor = fileDescriptor ;
    }
    return status ;
}

/*==================================================================================================================================
! . Local procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Write step.
!---------------------------------------------------------------------------------------------------------------------------------*/
static DCDStatus DCDWrite_Step (       fio_fd   fd            ,
                                       Boolean  hasUnitCell   ,
                                       Integer  currentFrame  ,
                                       Integer  currentStep   ,
                                       Integer  numberOfAtoms ,
                                 const Float32 *x             ,
                                 const Float32 *y             ,
                                 const Float32 *z             ,
                                 const Float64 *unitCell      )
{
    Integer32 numberOfBytes ;
    /* . Unit cell. */
    if ( hasUnitCell )
    {
        numberOfBytes = 48 ;
        FastIO_WriteInteger32 ( fd, numberOfBytes ) ;
        WRITE ( fd, unitCell, numberOfBytes ) ;
        FastIO_WriteInteger32 ( fd, numberOfBytes ) ;
    }
    /* . Coordinates. */
    numberOfBytes = 4 * numberOfAtoms ;
    FastIO_WriteInteger32 ( fd, numberOfBytes ) ;
    if ( fio_fwrite ( ( void * ) x, numberOfBytes, 1, fd ) != 1 ) return DCDStatus_BadWrite ;
    FastIO_WriteInteger32 ( fd, numberOfBytes ) ;
    FastIO_WriteInteger32 ( fd, numberOfBytes ) ;
    if ( fio_fwrite ( ( void * ) y, numberOfBytes, 1, fd ) != 1 ) return DCDStatus_BadWrite ;
    FastIO_WriteInteger32 ( fd, numberOfBytes ) ;
    FastIO_WriteInteger32 ( fd, numberOfBytes ) ;
    if ( fio_fwrite ( ( void * ) z, numberOfBytes, 1, fd ) != 1 ) return DCDStatus_BadWrite ;
    FastIO_WriteInteger32 ( fd, numberOfBytes ) ;
    /* . Update the header information. */
    fio_fseek             ( fd, NFILE_POS, FIO_SEEK_SET ) ;
    FastIO_WriteInteger32 ( fd, currentFrame            ) ;
    fio_fseek             ( fd, NSTEP_POS, FIO_SEEK_SET ) ;
    FastIO_WriteInteger32 ( fd, currentStep             ) ;
    fio_fseek             ( fd, 0, FIO_SEEK_END         ) ;
    return DCDStatus_Normal ;
}
