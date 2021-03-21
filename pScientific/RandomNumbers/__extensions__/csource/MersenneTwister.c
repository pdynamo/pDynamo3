/*==================================================================================================================================
! . This module implements a Mersenne-Twister random number generator.
!=================================================================================================================================*/

# include <stdio.h>
# include <stdlib.h>

# include "Cardinal.h"
# include "Integer.h"
# include "Memory.h"
# include "RandomNumberGenerator.h"
# include "Real.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Period parameters. */
# define _MersenneTwister_M 397
# define _MersenneTwister_N 624

/* . Least significant r bits. */
static const Cardinal _LowerMask = 0x7fffffffUL ;

/* . Most significant w-r bits. */
static const Cardinal _UpperMask = 0x80000000UL ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Structures.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The Mersenne-Twister structure. */
typedef struct {
    Cardinal mt[_MersenneTwister_N] ;
    Integer  mti ;
} MersenneTwister ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void *MersenneTwister_Allocate ( void )
{
    MersenneTwister *self = Memory_AllocateType ( MersenneTwister ) ;
    return ( void * ) self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void MersenneTwister_Deallocate ( void **self )
{
    if ( (*self) != NULL ) Memory_Deallocate ( (*self) ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Next cardinal.
!---------------------------------------------------------------------------------------------------------------------------------*/
#define MAGIC(y) ( ( ( y ) & 0x1 ) ? 0x9908b0dfUL : 0 )
static Cardinal MersenneTwister_NextCardinal ( void *vSelf )
{
    Cardinal k = 0 ;
    if ( vSelf != NULL )
    {
        MersenneTwister *self = ( MersenneTwister * ) vSelf ;
        Cardinal        *mt   = self->mt ;
        if ( self->mti >= _MersenneTwister_N )
        {
            Integer kk;
            for ( kk = 0 ; kk < _MersenneTwister_N - _MersenneTwister_M ; kk++ )
            {
                Cardinal y = ( mt[kk] & _UpperMask ) | ( mt[kk + 1] & _LowerMask ) ;
                mt[kk]     = mt[kk + _MersenneTwister_M] ^ ( y >> 1 ) ^ MAGIC ( y ) ;
            }
            for ( ; kk < _MersenneTwister_N - 1; kk++ )
            {
                Cardinal y = ( mt[kk] & _UpperMask ) | ( mt[kk + 1] & _LowerMask ) ;
                mt[kk]     = mt[kk + ( _MersenneTwister_M - _MersenneTwister_N ) ] ^ ( y >> 1 ) ^ MAGIC ( y ) ;
            }
            {
                Cardinal y = ( mt[_MersenneTwister_N - 1] & _UpperMask ) | ( mt[0] & _LowerMask ) ;
                mt[_MersenneTwister_N - 1] = mt[_MersenneTwister_M - 1] ^ ( y >> 1 ) ^ MAGIC ( y ) ;
            }
            self->mti = 0 ;
        }
        k = mt[self->mti] ;
        k ^= ( k >> 11 ) ;
        k ^= ( k <<  7 ) & 0x9d2c5680UL ;
        k ^= ( k << 15 ) & 0xefc60000UL ;
        k ^= ( k >> 18 ) ;
        self->mti++ ;
    }
    return k;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Next real on the half-open interval [0,1).
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real MersenneTwister_NextReal ( void *vSelf ) { return MersenneTwister_NextCardinal ( vSelf ) / 4294967296.0 ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set the seed.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MersenneTwister_SetSeed ( void *vSelf, Cardinal seed )
{
    if ( vSelf != NULL )
    {
        Integer i ;
        MersenneTwister *self = ( MersenneTwister * ) vSelf ;
        if ( seed == 0 ) seed = 4357 ;   /* . This is the default. */
        self->mt[0] = seed & 0xffffffffUL;
        for ( i = 1 ; i < _MersenneTwister_N ; i++ )
        {
            self->mt[i]  = ( 1812433253UL * ( self->mt[i-1] ^ ( self->mt[i-1] >> 30 ) ) + i ) ;
            self->mt[i] &= 0xffffffffUL ;
        }
        self->mti = i;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Type definition.
!---------------------------------------------------------------------------------------------------------------------------------*/
static const RandomNumberGeneratorType RNGType_MersenneTwister = { "Mersenne-Twister"            ,
                                                                   0xffffffffUL                  , /* . 2^32 - 1. */
                                                                   0                             ,
                                                                   sizeof ( MersenneTwister )    ,
                                                                   &MersenneTwister_Allocate     ,
                                                                   &MersenneTwister_Deallocate   ,
                                                                   &MersenneTwister_NextCardinal ,
                                                                   &MersenneTwister_NextReal     ,
                                                                   &MersenneTwister_SetSeed      } ;
const RandomNumberGeneratorType *RandomNumberGeneratorType_MersenneTwister = &RNGType_MersenneTwister ;
