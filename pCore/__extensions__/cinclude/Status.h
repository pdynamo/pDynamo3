/*==================================================================================================================================
! . The status type is used to indicate program state.
! . The messages differ in severity - some may be recoverable and others not.
!=================================================================================================================================*/
# ifndef _STATUS
# define _STATUS

/*

It would be possible to replace the current Status enum with a structure.

__FILE__ and __LINE__ macros useful?

Also assert?

*/

# include "Boolean.h"
# include "Size.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
typedef enum {
    Status_OK                    =  0 ,
    Status_AlgorithmError        =  1 ,
    Status_IndexOutOfRange       =  2 ,
    Status_InvalidArgument       =  3 ,
    Status_InvalidArrayOperation =  4 ,
    Status_MathError             =  5 ,
    Status_NonConformableArrays  =  6 ,
    Status_OutOfMemory           =  7
} Status ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Macros.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Check whether a pointer to a status flag is OK. */
# define Status_IsOK( self ) ( ( (self) == NULL ) || ( (*(self)) == Status_OK ) )

/* . Check whether a status flag is OK. */
# define Status_IsValueOK( self ) ( (self) == Status_OK )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Functions.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*extern Boolean Status_IsOK ( Status *self ) ;*/
extern void Status_Set  ( Status *self   ,
                          Status  status ) ;
extern Size Status_Size ( void           ) ;

# endif
