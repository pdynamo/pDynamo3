/*==================================================================================================================================
! . Double symmetric arrays.
!
! . These are symmetric four-dimensional arrays. The following items are equivalent: ijkl, ijlk, jikl, jilk, klij, klji, lkij, lkji.
!
! . Packing is similar to SymmetricMatrix: i >= j, k >= l and (ij) >= (kl).
!
!=================================================================================================================================*/

# include <stdio.h>
# include <stdlib.h>

# include "Array_Macros.h"
# include "DoubleSymmetricMatrix.h"
# include "Iterator1D.h"
# include "Memory.h"
# include "NumericalMacros.h"

# define _UseCBLAS
# include "Utilities1D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Raw allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
DoubleSymmetricMatrix *DoubleSymmetricMatrix_Allocate ( Status *status )
{
    DoubleSymmetricMatrix *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( DoubleSymmetricMatrix ) ;
        if ( self != NULL )
        {
            self->block    = NULL ;
            self->data     = NULL ;
            self->extent   = 0 ;
            self->extent01 = 0 ;
            self->size     = 0 ;
        }
        else Status_Set ( status, Status_OutOfMemory ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation with extent.
!---------------------------------------------------------------------------------------------------------------------------------*/
DoubleSymmetricMatrix *DoubleSymmetricMatrix_AllocateWithExtent ( const Integer extent, Status *status )
{
    DoubleSymmetricMatrix *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        auto RealBlock *block = NULL ;
        block = RealBlock_Allocate ( DoubleSymmetricMatrix_ViewSize ( extent ), status ) ;
        self  = DoubleSymmetricMatrix_FromExtentBlock ( extent, block, True   , status ) ;
        if ( self == NULL ) RealBlock_Deallocate ( &block ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deep cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
DoubleSymmetricMatrix *DoubleSymmetricMatrix_CloneDeep ( const DoubleSymmetricMatrix *self, Status *status )
{
    DoubleSymmetricMatrix *clone = NULL ;
    if ( self != NULL )
    {
        clone = DoubleSymmetricMatrix_AllocateWithExtent ( self->extent, status ) ;
        DoubleSymmetricMatrix_CopyTo ( self, clone, status ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Shallow cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
DoubleSymmetricMatrix *DoubleSymmetricMatrix_CloneShallow ( const DoubleSymmetricMatrix *self, Status *status )
{
    DoubleSymmetricMatrix *clone = NULL ;
    if ( self != NULL ) clone = DoubleSymmetricMatrix_FromExtentBlock ( self->extent, self->block, True, status ) ;
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Copying.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DoubleSymmetricMatrix_CopyTo ( const DoubleSymmetricMatrix *self, DoubleSymmetricMatrix *other, Status *status )
{
    if ( ( self != NULL ) && ( other != NULL ) && Status_IsOK ( status ) )
    {
        Utilities1D_CopyTo ( self->size, self->data, 1, other->size, other->data, 1, status ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DoubleSymmetricMatrix_Deallocate ( DoubleSymmetricMatrix **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        RealBlock_Deallocate ( &((*self)->block) ) ;
        Memory_Deallocate    (   (*self)         ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Constructor from an extent and, optionally, a block.
!---------------------------------------------------------------------------------------------------------------------------------*/
DoubleSymmetricMatrix *DoubleSymmetricMatrix_FromExtentBlock ( const Integer    extent        ,
                                                                     RealBlock *block         ,
                                                               const Boolean    withReference ,
                                                                     Status    *status        )
{
    DoubleSymmetricMatrix *self = DoubleSymmetricMatrix_Allocate ( status ) ;
    if ( ( self != NULL ) && ( extent > 0 ) )
    {
        self->extent   = extent ;
        self->extent01 = ( extent * ( extent + 1 ) ) / 2 ;
        self->size     = self->extent01 * self->extent01 ;
        if ( block != NULL )
        {
            if ( withReference ) { Array_AssignBlockWithReference    ( self, block, 0 ) ; }
            else                 { Array_AssignBlockWithoutReference ( self, block, 0 ) ; }
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real DoubleSymmetricMatrix_GetItem ( const DoubleSymmetricMatrix *self   ,
                                     const Integer                i      ,
                                     const Integer                j      ,
                                     const Integer                k      ,
                                     const Integer                l      ,
                                           Status                *status )
{
    Real value = 0.0e+00 ;
    if ( self != NULL )
    {
        auto Integer  ijkl ;
        ijkl = DoubleSymmetricMatrix_Index ( i, j, k, l ) ;
        if ( ijkl < self->size ) value = self->data[ijkl] ;
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
    return value ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Increment an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DoubleSymmetricMatrix_IncrementItem ( const DoubleSymmetricMatrix *self   ,
                                           const Integer                i      ,
                                           const Integer                j      ,
                                           const Integer                k      ,
                                           const Integer                l      ,
                                           const Real                   value  ,
                                                 Status                *status )
{
    if ( self != NULL )
    {
        auto Integer  ijkl ;
        ijkl = DoubleSymmetricMatrix_Index ( i, j, k, l ) ;
        if ( ijkl < self->size ) self->data[ijkl] += value ;
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get an index into the array.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define IJIndex(i,j) ( ( i * ( i + 1 ) ) / 2 + j )
Integer  DoubleSymmetricMatrix_Index ( const Integer  i, const Integer  j, const Integer  k, const Integer  l )
{
    auto Integer  pqrs, p, pq, q, r, rs, s ;
    p  = Maximum ( i, j ) ; q  = Minimum ( i, j ) ;
    r  = Maximum ( k, l ) ; s  = Minimum ( k, l ) ;
    pq = IJIndex ( p, q ) ; rs = IJIndex ( r, s ) ;
    if ( pq >= rs ) pqrs = IJIndex ( pq, rs ) ;
    else            pqrs = IJIndex ( rs, pq ) ;
    return pqrs ;
}
# undef IJIndex

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the default iterator.
!---------------------------------------------------------------------------------------------------------------------------------*/
Iterator *DoubleSymmetricMatrix_MakeIterator ( const DoubleSymmetricMatrix *self, Status *status )
{
    Iterator *iterator = NULL ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto Iterator1D *iterator1D ; 
        iterator   = Iterator_Allocate   ( status ) ;
        iterator1D = Iterator1D_Allocate ( status ) ;
        Iterator1D_Initialize            ( iterator1D, 0, self->size, 1, status ) ;
        Iterator1D_MakeIterator          ( iterator1D, iterator ) ;
    }
    return iterator ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Printing.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DoubleSymmetricMatrix_Print ( const DoubleSymmetricMatrix *self )
{
    if ( self == NULL ) printf ( "Null double symmetric matrix.\n" ) ;
    else
    {
        auto Integer  i, ij, ijkl, j, k, kl, l, upper ;
        for ( i = ij = ijkl = 0 ; i < self->extent ; i++ )
        {
            for ( j = 0 ; j <= i ; ij++, j++ )
            {
                for ( k = kl = 0 ; k <= i ; k++ )
                {
                    if ( k == i ) upper = j ;
                    else          upper = k ;
                    for ( l = 0 ; l <= upper ; ijkl++, kl++, l++ )
                    {
                        printf ( "%5d %5d %5d %5d %5d %5d %5d %12.6f\n", ijkl, i, j, k, l, ij, kl, self->data[ijkl] ) ;
                    }
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set all items.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DoubleSymmetricMatrix_Set ( DoubleSymmetricMatrix *self, Real value )
{
    if ( self != NULL ) Utilities1D_Set ( self->size, self->data, 1, value ) ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set an item.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DoubleSymmetricMatrix_SetItem ( const DoubleSymmetricMatrix *self   ,
                                     const Integer                i      ,
                                     const Integer                j      ,
                                     const Integer                k      ,
                                     const Integer                l      ,
                                     const Real                   value  ,
                                           Status                *status )
{
    if ( self != NULL )
    {
        auto Integer  ijkl ;
        ijkl = DoubleSymmetricMatrix_Index ( i, j, k, l ) ;
        if ( ijkl < self->size ) self->data[ijkl] = value ;
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Unweight the matrix.
!---------------------------------------------------------------------------------------------------------------------------------*/
void DoubleSymmetricMatrix_Unweight ( DoubleSymmetricMatrix *self )
{
    if ( ( self != NULL ) && ( self->extent > 0 ) )
    {
        auto Integer  i, ij, ijkl, j, k, kl, l, upper ;
        auto Real     w ;
        for ( i = ij = ijkl = 0 ; i < self->extent ; i++ )
        {
            for ( j = 0 ; j <= i ; ij++, j++ )
            {
                for ( k = kl = 0 ; k <= i ; k++ )
                {
                    if ( k == i ) upper = j ;
                    else          upper = k ;
                    for ( l = 0 ; l <= upper ; ijkl++, kl++, l++ )
                    {
                        w = 0.125e+00 ;
                        if ( i  == j  ) w *= 2.0e+00 ;
                        if ( k  == l  ) w *= 2.0e+00 ;
                        if ( ij == kl ) w *= 2.0e+00 ;
                        self->data[ijkl] *= w ;
                    }
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . The view size.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer DoubleSymmetricMatrix_ViewSize ( const Integer extent )
{
    Integer size = 0 ;
    if ( extent > 0 )
    {
        auto Integer extent01 ;
        extent01 = ( extent * ( extent + 1 ) ) / 2 ;
        size     = extent01 * extent01 ;
    }
    return size ;
}

# undef _UseCBLAS
