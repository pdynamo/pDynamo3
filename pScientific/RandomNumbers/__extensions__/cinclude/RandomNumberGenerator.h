# ifndef _RANDOMNUMBERGENERATOR
# define _RANDOMNUMBERGENERATOR

# include "Boolean.h"
# include "Cardinal.h"
# include "MachineTypes.h"
# include "Real.h"
# include "Size.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The random number generator type. */
typedef struct {
    const Character *name    ;
          Cardinal   maximum ;
          Cardinal   minimum ;
          Size       size    ;
          void    *( * Allocate     ) ( void ) ;
          void     ( * Deallocate   ) ( void **vState ) ;
          Cardinal ( * NextCardinal ) ( void  *vState ) ;
          Real     ( * NextReal     ) ( void  *vState ) ;
          void     ( * SetSeed      ) ( void  *vState, Cardinal seed ) ;
} RandomNumberGeneratorType ;

/* . The random number generator. */
typedef struct {
          Boolean                    hasGaussian ;
          Real                       gaussian    ;
    const RandomNumberGeneratorType *type        ;
          void                      *vState      ;
} RandomNumberGenerator ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedure declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern RandomNumberGenerator *RandomNumberGenerator_Allocate     ( const RandomNumberGeneratorType *type ) ;
extern RandomNumberGenerator *RandomNumberGenerator_Clone        ( const RandomNumberGenerator     *self ) ;
extern void                   RandomNumberGenerator_Deallocate   (       RandomNumberGenerator    **self ) ;
extern Cardinal               RandomNumberGenerator_NextCardinal ( const RandomNumberGenerator     *self ) ;
extern Real                   RandomNumberGenerator_NextReal     ( const RandomNumberGenerator     *self ) ;
extern Real                   RandomNumberGenerator_NextRealOpen ( const RandomNumberGenerator     *self ) ;
extern void                   RandomNumberGenerator_SetSeed      ( const RandomNumberGenerator     *self, Cardinal seed ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Type declarations.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern const RandomNumberGeneratorType *RandomNumberGeneratorType_MersenneTwister ;

# endif
