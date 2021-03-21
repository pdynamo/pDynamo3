/*==================================================================================================================================
! . The block storage module handles cardinal and real data stored in blocks.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "BlockStorage.h"
# include "Memory.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Definitions.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define _BlockStorage_DefaultSize 1024

/*==================================================================================================================================
! . Blocks.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
Block *Block_Allocate ( const Integer  blockSize         ,
                        const Integer  numberOfIndices16 ,
                        const Integer  numberOfIndices32 ,
                        const Integer  numberOfReal      ,
                              Status  *status            )
{
    Block *self = NULL ;
    if (   ( blockSize         > 0 )   &&
         ( ( numberOfIndices16 > 0 ) ||
           ( numberOfIndices32 > 0 ) ||
           ( numberOfReal      > 0 ) ) &&
           Status_IsOK ( status ) )
    {
        auto Boolean isOK = False ;
        self = Memory_AllocateType ( Block ) ;
        if ( self != NULL )
        {
            isOK            = True ;
            self->count     = 0    ;
            self->data      = NULL ;
            self->indices16 = NULL ;
            self->indices32 = NULL ;
            if ( numberOfIndices16 > 0 )
            {
                self->indices16 = Memory_AllocateArrayOfTypes ( blockSize * numberOfIndices16, Cardinal16 ) ;
                isOK            = isOK && ( self->indices16 != NULL ) ;
            }
            if ( numberOfIndices32 > 0 )
            {
                self->indices32 = Memory_AllocateArrayOfTypes ( blockSize * numberOfIndices32, Cardinal32 ) ;
                isOK            = isOK && ( self->indices32 != NULL ) ;
            }
            if ( numberOfReal      > 0 )
            {
                self->data      = Memory_AllocateArrayOfTypes ( blockSize * numberOfReal     , Real       ) ;
                isOK            = isOK && ( self->data != NULL ) ;
            }
        }
        if ( ! isOK )
        {
            Block_Deallocate ( &self ) ;
            Status_Set ( status, Status_OutOfMemory ) ;
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Block_Deallocate ( Block **self )
{
    if ( (*self) != NULL )
    {
        Memory_Deallocate ( (*self)->data      ) ;
        Memory_Deallocate ( (*self)->indices16 ) ;
        Memory_Deallocate ( (*self)->indices32 ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

static void Block_DeallocateVoid ( void *vSelf )
{
    if ( vSelf != NULL ) { Block *self ; self = ( Block * ) vSelf ; Block_Deallocate ( &self ) ; }
}

/*==================================================================================================================================
! . Block storage.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
BlockStorage *BlockStorage_Allocate ( Status *status )
{
    BlockStorage *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        auto Boolean isOK = False ;
        self = Memory_AllocateType ( BlockStorage ) ;
        if ( self != NULL )
        {
            isOK = True ;
            self->blockSize      = _BlockStorage_DefaultSize ;
            self->count          = 0 ;
            self->nIndices16     = 0 ;
            self->nIndices32     = 0 ;
            self->nReal          = 0 ;
            self->checkUnderFlow = False ;
            self->underFlow      = 0.0e+00 ;
            self->blocks         = List_Allocate ( ) ;
            if ( self->blocks == NULL ) isOK = False ;
            else self->blocks->Element_Deallocate = Block_DeallocateVoid ;
        }
        if ( ! isOK )
        {
            BlockStorage_Deallocate ( &self ) ;
            Status_Set ( status, Status_OutOfMemory ) ;
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Add data to the block storage.
!---------------------------------------------------------------------------------------------------------------------------------*/
void BlockStorage_AddData (       BlockStorage *self      ,
                            const Integer       count     ,
                            const Real         *data      ,
                            const Cardinal16   *indices16 ,
                            const Cardinal32   *indices32 ,
                                  Status       *status    )
{
    if (   ( self      != NULL ) &&
           ( count     >  0    ) &&
         ( ( data      != NULL ) ||
           ( indices16 != NULL ) ||
           ( indices32 != NULL ) ) &&
           Status_IsOK ( status ) )
    {
        auto Block   *block ;
        auto Integer  i, j ;
        /* . Get the current block. */
        if ( self->blocks->last == NULL )
        {
            block = Block_Allocate ( self->blockSize, self->nIndices16, self->nIndices32, self->nReal, status ) ;
            if ( block == NULL ) goto FinishUp ;
            List_Element_Append ( self->blocks, ( void * ) block ) ;
        }
        else block = ( Block * ) self->blocks->last->node ;
        /* . Copy data of a certain size only - data are excluded only if all the data of a certain count underflow. */
        if ( self->checkUnderFlow )
        {
            auto Boolean isOK ;
            for ( i = 0 ; i < count ; i++ )
            {
                if ( block->count >= self->blockSize )
                {
                    block = Block_Allocate ( self->blockSize, self->nIndices16, self->nIndices32, self->nReal, status ) ;
                    if ( block == NULL ) goto FinishUp ;
                    List_Element_Append ( self->blocks, ( void * ) block ) ;
                }
                isOK = False ;
                for ( j = 0 ; j < self->nReal ; j++ )
                {
                    if ( fabs ( data[self->nReal*i+j] ) > self->underFlow ) { isOK = True ; break ; }
                }
                if ( isOK )
                {
                    for ( j = 0 ; j < self->nIndices16 ; j++ ) block->indices16[self->nIndices16*block->count+j] = indices16[self->nIndices16*i+j] ;
                    for ( j = 0 ; j < self->nIndices32 ; j++ ) block->indices32[self->nIndices32*block->count+j] = indices32[self->nIndices32*i+j] ;
                    for ( j = 0 ; j < self->nReal      ; j++ ) block->data     [self->nReal     *block->count+j] = data     [self->nReal     *i+j] ;
                    block->count++ ;
                    self->count++ ;
                }
            }
        }
        /* . Copy all data. */
        else
        {
            for ( i = 0 ; i < count ; i++ )
            {
                if ( block->count >= self->blockSize )
                {
                    block = Block_Allocate ( self->blockSize, self->nIndices16, self->nIndices32, self->nReal, status ) ;
                    if ( block == NULL ) goto FinishUp ;
                    List_Element_Append ( self->blocks, ( void * ) block ) ;
                }
                for ( j = 0 ; j < self->nIndices16 ; j++ ) block->indices16[self->nIndices16*block->count+j] = indices16[self->nIndices16*i+j] ;
                for ( j = 0 ; j < self->nIndices32 ; j++ ) block->indices32[self->nIndices32*block->count+j] = indices32[self->nIndices32*i+j] ;
                for ( j = 0 ; j < self->nReal      ; j++ ) block->data     [self->nReal     *block->count+j] = data     [self->nReal     *i+j] ;
                block->count++ ;
           }
           self->count += count ;
        }
FinishUp:
        ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Size of the block storage in bytes.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real BlockStorage_ByteSize ( BlockStorage *self )
{
    Real size = 0.0e+00 ;
    if  ( self != NULL )
    {
        auto Block *block = NULL ;
        /* . Initialization. */
        size = sizeof ( BlockStorage ) ;
        /* . Blocks. */
        if ( self->blocks != NULL )
        {
            size += sizeof ( List ) ;
            List_Iterate_Initialize ( self->blocks ) ;
            while ( ( block = BlockStorage_Iterate ( self ) ) != NULL )
            {
                size += block->count * ( self->nIndices16 * sizeof ( Cardinal16 ) + self->nIndices32 * sizeof ( Cardinal32 ) + self->nReal * sizeof ( Real ) ) ;
            }
        }
    }
    return size ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Count.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer BlockStorage_Count ( BlockStorage *self ) { return ( ( self == NULL ) ? 0 : self->count ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void BlockStorage_Deallocate ( BlockStorage **self )
{
    if ( (*self) != NULL )
    {
       List_Deallocate ( &((*self)->blocks) ) ;
       Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Emptying.
!---------------------------------------------------------------------------------------------------------------------------------*/
void BlockStorage_Empty ( BlockStorage *self )
{
    if ( self != NULL )
    {
       List_Empty ( self->blocks ) ;
       self->count = 0 ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Iterate over the interactions in a list.
!---------------------------------------------------------------------------------------------------------------------------------*/
Block *BlockStorage_Iterate ( BlockStorage *self )
{
    if ( self == NULL ) return NULL ;
    else                return ( Block * ) List_Iterate ( self->blocks ) ;
}
