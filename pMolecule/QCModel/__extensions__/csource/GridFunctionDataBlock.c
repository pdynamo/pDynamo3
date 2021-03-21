/*==================================================================================================================================
! . This module defines a data structure for storing function and derivative values on a grid.
!=================================================================================================================================*/

# include "GridFunctionDataBlock.h"
# include "Memory.h"
# include "NumericalMacros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
GridFunctionDataBlock *GridFunctionDataBlock_Allocate ( const Integer  numberOfFunctions ,
                                                        const Integer  numberOfPoints    ,
                                                        const Integer  order             ,
                                                              Status  *status            )
{
    GridFunctionDataBlock *self = NULL ;
    if ( Status_IsOK ( status ) )
    {
        self = Memory_AllocateType ( GridFunctionDataBlock ) ;
        if ( self != NULL )
        {
            auto Integer f, o, p ;
            f = Maximum ( numberOfFunctions   , 0 ) ; self->numberOfFunctions = f ;
            p = Maximum ( numberOfPoints      , 0 ) ; self->numberOfPoints    = p ;
            o = Minimum ( Maximum ( order, 0 ), 3 ) ; self->order             = o ;
            self->indices = IntegerArray1D_AllocateWithExtent ( f,    status ) ;
            self->f       = RealArray2D_AllocateWithExtents   ( f, p, status ) ;
            if ( o > 0 )
            {
               self->fX   = RealArray2D_AllocateWithExtents ( f, p, status ) ;
               self->fY   = RealArray2D_AllocateWithExtents ( f, p, status ) ;
               self->fZ   = RealArray2D_AllocateWithExtents ( f, p, status ) ;
            }
            else { self->fX = NULL ; self->fY = NULL ; self->fZ = NULL ; }
            if ( o > 1 )
            {
               self->fXX  = RealArray2D_AllocateWithExtents ( f, p, status ) ;
               self->fXY  = RealArray2D_AllocateWithExtents ( f, p, status ) ;
               self->fXZ  = RealArray2D_AllocateWithExtents ( f, p, status ) ;
               self->fYY  = RealArray2D_AllocateWithExtents ( f, p, status ) ;
               self->fYZ  = RealArray2D_AllocateWithExtents ( f, p, status ) ;
               self->fZZ  = RealArray2D_AllocateWithExtents ( f, p, status ) ;
            }
            else
            {
                self->fXX = NULL ; self->fXY = NULL ; self->fXZ = NULL ;
                self->fYY = NULL ; self->fYZ = NULL ; self->fZZ = NULL ;
            }
            if ( o > 2 )
            {
               self->fXXX = RealArray2D_AllocateWithExtents ( f, p, status ) ;
               self->fXXY = RealArray2D_AllocateWithExtents ( f, p, status ) ;
               self->fXXZ = RealArray2D_AllocateWithExtents ( f, p, status ) ;
               self->fXYY = RealArray2D_AllocateWithExtents ( f, p, status ) ;
               self->fXYZ = RealArray2D_AllocateWithExtents ( f, p, status ) ;
               self->fXZZ = RealArray2D_AllocateWithExtents ( f, p, status ) ;
               self->fYYY = RealArray2D_AllocateWithExtents ( f, p, status ) ;
               self->fYYZ = RealArray2D_AllocateWithExtents ( f, p, status ) ;
               self->fYZZ = RealArray2D_AllocateWithExtents ( f, p, status ) ;
               self->fZZZ = RealArray2D_AllocateWithExtents ( f, p, status ) ;
            }
            else
            {
                self->fXXX = NULL ; self->fXXY = NULL ; self->fXXZ = NULL ; self->fXYY = NULL ; self->fXYZ = NULL ;
                self->fXZZ = NULL ; self->fYYY = NULL ; self->fYYZ = NULL ; self->fYZZ = NULL ; self->fZZZ = NULL ;
            }
            if ( ! Status_IsOK ( status ) ) GridFunctionDataBlock_Deallocate ( &self ) ;
        }
    }
    if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . The size in bytes of the data block.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real GridFunctionDataBlock_ByteSize ( const GridFunctionDataBlock *self )
{
    Real size = 0.0e+00 ;
    if ( self != NULL )
    {
        auto Integer f, n, p ;
        f = View2D_Rows    ( self->f ) ;
        p = View2D_Columns ( self->f ) ;
        n = 1 ;
        if ( self->order > 0 ) n +=  3 ;
        if ( self->order > 1 ) n +=  6 ;
        if ( self->order > 2 ) n += 10 ;
        size  = sizeof ( GridFunctionDataBlock ) + ( sizeof ( IntegerArray1D ) + sizeof ( Integer ) * f ) + n * ( sizeof ( RealArray2D ) + sizeof ( Real ) * f * p ) ;
    }
    return size ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GridFunctionDataBlock_Deallocate ( GridFunctionDataBlock **self )
{
    if ( ( self != NULL ) && ( (*self) != NULL ) )
    {
        IntegerArray1D_Deallocate ( &((*self)->indices) ) ;
        RealArray2D_Deallocate    ( &((*self)->f      ) ) ;
        RealArray2D_Deallocate    ( &((*self)->fX     ) ) ;
        RealArray2D_Deallocate    ( &((*self)->fY     ) ) ;
        RealArray2D_Deallocate    ( &((*self)->fZ     ) ) ;
        RealArray2D_Deallocate    ( &((*self)->fXX    ) ) ;
        RealArray2D_Deallocate    ( &((*self)->fXY    ) ) ;
        RealArray2D_Deallocate    ( &((*self)->fXZ    ) ) ;
        RealArray2D_Deallocate    ( &((*self)->fYY    ) ) ;
        RealArray2D_Deallocate    ( &((*self)->fYZ    ) ) ;
        RealArray2D_Deallocate    ( &((*self)->fZZ    ) ) ;
        RealArray2D_Deallocate    ( &((*self)->fXXX   ) ) ;
        RealArray2D_Deallocate    ( &((*self)->fXXY   ) ) ;
        RealArray2D_Deallocate    ( &((*self)->fXXZ   ) ) ;
        RealArray2D_Deallocate    ( &((*self)->fXYY   ) ) ;
        RealArray2D_Deallocate    ( &((*self)->fXYZ   ) ) ;
        RealArray2D_Deallocate    ( &((*self)->fXZZ   ) ) ;
        RealArray2D_Deallocate    ( &((*self)->fYYY   ) ) ;
        RealArray2D_Deallocate    ( &((*self)->fYYZ   ) ) ;
        RealArray2D_Deallocate    ( &((*self)->fYZZ   ) ) ;
        RealArray2D_Deallocate    ( &((*self)->fZZZ   ) ) ;
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Filtering.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GridFunctionDataBlock_FilterValues ( GridFunctionDataBlock *self, const Integer fStart, const Real *tolerance )
{
    if ( ( self != NULL ) && ( tolerance != NULL ) && ( (*tolerance) > 0.0e+00 ) )
    {
        auto Integer     f, f0, f1, n, o ;
        auto Real        t        ;
        auto RealArray1D new, old ;
        f0 = Maximum ( fStart, 0 ) ;
        f1 = self->numberOfFunctions ;
        o  = self->order ;
        t  = (*tolerance) ;
        for ( f = n = f0 ; f < f1 ; f++ )
        {
            RealArray2D_RowView ( self->f, f, False, &old, NULL ) ;
            if ( RealArray1D_AbsoluteMaximum ( &old ) > t )
            {
                if ( f != n )
                {
                    Array1D_Item ( self->indices, n ) = Array1D_Item ( self->indices, f ) ;
                    RealArray2D_RowView ( self->f, n, False, &new, NULL ) ;
                    RealArray1D_CopyTo  ( &old, &new, NULL ) ;
                    if ( o > 0 )
                    {
                        RealArray2D_RowView ( self->fX, f, False, &old, NULL ) ; RealArray2D_RowView ( self->fX, n, False, &new, NULL ) ; RealArray1D_CopyTo  ( &old, &new, NULL ) ;
                        RealArray2D_RowView ( self->fY, f, False, &old, NULL ) ; RealArray2D_RowView ( self->fY, n, False, &new, NULL ) ; RealArray1D_CopyTo  ( &old, &new, NULL ) ;
                        RealArray2D_RowView ( self->fZ, f, False, &old, NULL ) ; RealArray2D_RowView ( self->fZ, n, False, &new, NULL ) ; RealArray1D_CopyTo  ( &old, &new, NULL ) ;
                    }
                    if ( o > 1 )
                    {
                        RealArray2D_RowView ( self->fXX, f, False, &old, NULL ) ; RealArray2D_RowView ( self->fXX, n, False, &new, NULL ) ; RealArray1D_CopyTo  ( &old, &new, NULL ) ;
                        RealArray2D_RowView ( self->fXY, f, False, &old, NULL ) ; RealArray2D_RowView ( self->fXY, n, False, &new, NULL ) ; RealArray1D_CopyTo  ( &old, &new, NULL ) ;
                        RealArray2D_RowView ( self->fXZ, f, False, &old, NULL ) ; RealArray2D_RowView ( self->fXZ, n, False, &new, NULL ) ; RealArray1D_CopyTo  ( &old, &new, NULL ) ;
                        RealArray2D_RowView ( self->fYY, f, False, &old, NULL ) ; RealArray2D_RowView ( self->fYY, n, False, &new, NULL ) ; RealArray1D_CopyTo  ( &old, &new, NULL ) ;
                        RealArray2D_RowView ( self->fYZ, f, False, &old, NULL ) ; RealArray2D_RowView ( self->fYZ, n, False, &new, NULL ) ; RealArray1D_CopyTo  ( &old, &new, NULL ) ;
                        RealArray2D_RowView ( self->fZZ, f, False, &old, NULL ) ; RealArray2D_RowView ( self->fZZ, n, False, &new, NULL ) ; RealArray1D_CopyTo  ( &old, &new, NULL ) ;
                    }
                    if ( o > 2 )
                    {
                        RealArray2D_RowView ( self->fXXX, f, False, &old, NULL ) ; RealArray2D_RowView ( self->fXXX, n, False, &new, NULL ) ; RealArray1D_CopyTo  ( &old, &new, NULL ) ;
                        RealArray2D_RowView ( self->fXXY, f, False, &old, NULL ) ; RealArray2D_RowView ( self->fXXY, n, False, &new, NULL ) ; RealArray1D_CopyTo  ( &old, &new, NULL ) ;
                        RealArray2D_RowView ( self->fXXZ, f, False, &old, NULL ) ; RealArray2D_RowView ( self->fXXZ, n, False, &new, NULL ) ; RealArray1D_CopyTo  ( &old, &new, NULL ) ;
                        RealArray2D_RowView ( self->fXYY, f, False, &old, NULL ) ; RealArray2D_RowView ( self->fXYY, n, False, &new, NULL ) ; RealArray1D_CopyTo  ( &old, &new, NULL ) ;
                        RealArray2D_RowView ( self->fXYZ, f, False, &old, NULL ) ; RealArray2D_RowView ( self->fXYZ, n, False, &new, NULL ) ; RealArray1D_CopyTo  ( &old, &new, NULL ) ;
                        RealArray2D_RowView ( self->fXZZ, f, False, &old, NULL ) ; RealArray2D_RowView ( self->fXZZ, n, False, &new, NULL ) ; RealArray1D_CopyTo  ( &old, &new, NULL ) ;
                        RealArray2D_RowView ( self->fYYY, f, False, &old, NULL ) ; RealArray2D_RowView ( self->fYYY, n, False, &new, NULL ) ; RealArray1D_CopyTo  ( &old, &new, NULL ) ;
                        RealArray2D_RowView ( self->fYYZ, f, False, &old, NULL ) ; RealArray2D_RowView ( self->fYYZ, n, False, &new, NULL ) ; RealArray1D_CopyTo  ( &old, &new, NULL ) ;
                        RealArray2D_RowView ( self->fYZZ, f, False, &old, NULL ) ; RealArray2D_RowView ( self->fYZZ, n, False, &new, NULL ) ; RealArray1D_CopyTo  ( &old, &new, NULL ) ;
                        RealArray2D_RowView ( self->fZZZ, f, False, &old, NULL ) ; RealArray2D_RowView ( self->fZZZ, n, False, &new, NULL ) ; RealArray1D_CopyTo  ( &old, &new, NULL ) ;
                    }
                }
                n++ ;
            }
        }
        self->numberOfFunctions = n ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GridFunctionDataBlock_Initialize ( GridFunctionDataBlock *self )
{
    if ( self != NULL )
    {
        auto Integer f ;
        self->numberOfFunctions = 0 ;
        for ( f = 0 ; f < View1D_Extent ( self->indices ) ; f++ ) Array1D_Item ( self->indices, f ) = f ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Resizing - always done (larger or smaller than existing).
!---------------------------------------------------------------------------------------------------------------------------------*/
void GridFunctionDataBlock_Resize ( GridFunctionDataBlock *self, const Integer numberOfFunctions, Status *status )
{
    if ( Status_IsOK ( status ) && ( self != NULL ) && ( View2D_Rows ( self->f ) != numberOfFunctions ) )
    {
        self->numberOfFunctions = numberOfFunctions ;
        IntegerArray1D_Resize ( self->indices, numberOfFunctions, status ) ;
        RealArray2D_Resize    ( self->f      , numberOfFunctions, status ) ;
        if ( self->order > 0 )
        {
            RealArray2D_Resize ( self->fX  , numberOfFunctions, status ) ;
            RealArray2D_Resize ( self->fY  , numberOfFunctions, status ) ;
            RealArray2D_Resize ( self->fZ  , numberOfFunctions, status ) ;
        }
        if ( self->order > 1 )
        {
            RealArray2D_Resize ( self->fXX , numberOfFunctions, status ) ;
            RealArray2D_Resize ( self->fXY , numberOfFunctions, status ) ;
            RealArray2D_Resize ( self->fXZ , numberOfFunctions, status ) ;
            RealArray2D_Resize ( self->fYY , numberOfFunctions, status ) ;
            RealArray2D_Resize ( self->fYZ , numberOfFunctions, status ) ;
            RealArray2D_Resize ( self->fZZ , numberOfFunctions, status ) ;
        }
        if ( self->order > 2 )
        {
            RealArray2D_Resize ( self->fXXX, numberOfFunctions, status ) ;
            RealArray2D_Resize ( self->fXXY, numberOfFunctions, status ) ;
            RealArray2D_Resize ( self->fXXZ, numberOfFunctions, status ) ;
            RealArray2D_Resize ( self->fXYY, numberOfFunctions, status ) ;
            RealArray2D_Resize ( self->fXYZ, numberOfFunctions, status ) ;
            RealArray2D_Resize ( self->fXZZ, numberOfFunctions, status ) ;
            RealArray2D_Resize ( self->fYYY, numberOfFunctions, status ) ;
            RealArray2D_Resize ( self->fYYZ, numberOfFunctions, status ) ;
            RealArray2D_Resize ( self->fYZZ, numberOfFunctions, status ) ;
            RealArray2D_Resize ( self->fZZZ, numberOfFunctions, status ) ;
        }
    }
}
