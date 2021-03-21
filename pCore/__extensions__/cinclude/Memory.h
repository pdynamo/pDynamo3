/*==================================================================================================================================
! . Memory handling.
!=================================================================================================================================*/
# ifndef _MEMORY
# define _MEMORY

# include <stdlib.h>
# include <string.h>

# include "Integer.h"
# include "Size.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define Memory_Allocate(                                  size           ) ( type * ) malloc (                   ( Size ) size )
# define Memory_AllocateArray(                     extent, size           )            calloc ( ( Size ) extent, ( Size ) size )
# define Memory_AllocateArrayOfReferences(         extent, type           ) ( type ** ) calloc ( ( Size ) extent, sizeof ( type * ) )
# define Memory_AllocateArrayOfTypes(              extent, type           ) ( type * ) calloc ( ( Size ) extent, sizeof ( type ) )
# define Memory_AllocateFlexibleArray(             extent, type, itemType ) ( type * ) malloc ( sizeof ( type ) + ( ( Size ) extent ) * sizeof ( itemType ) )
# define Memory_AllocateType(                              type           ) ( type * ) malloc (                   sizeof ( type ) )
# define Memory_Deallocate(                  self                         ) { free ( ( void * ) self ) ; self = NULL ; }                                   
# define Memory_ReallocateArray(             self, extent, size           ) ( type * ) realloc ( ( void * ) self, ( ( Size ) extent   * size          ) ) 
# define Memory_ReallocateArrayOfReferences( self, extent, type           ) ( type ** ) realloc ( ( void * ) self, ( ( Size ) extent ) * sizeof ( type * ) )
# define Memory_ReallocateArrayOfTypes(      self, extent, type           ) ( type * ) realloc ( ( void * ) self, ( ( Size ) extent ) * sizeof ( type ) )
# define Memory_ReallocateFlexibleArray(     self, extent, type, itemType ) ( type * ) realloc ( ( void * ) self, sizeof ( type ) + ( ( Size ) extent ) * sizeof ( itemType ) )

/* . Non-overlapping memory blocks - restrict is used in definition. */
# define Memory_CopyTo( self, other, size ) memcpy ( other, self, ( Size ) size )                                         

/* . Unsafe - ensure that iLocal is not required by one of the arguments. */
# define Memory_Set( self, extent, initializer ) if ( self != NULL ) { auto Integer  iLocal ; for ( iLocal = 0 ; iLocal < extent ; iLocal++ ) self[iLocal] = initializer ; }

# endif
