/*==================================================================================================================================
! . A container for Gaussian basis sets.
!=================================================================================================================================*/

# include "GaussianBasisContainer.h"
# include "Memory.h"
# include "NumericalMacros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
GaussianBasisContainer *GaussianBasisContainer_Allocate ( const Integer capacity, Status *status )
{
    GaussianBasisContainer *self = Memory_AllocateType ( GaussianBasisContainer ) ;
    if ( self != NULL )
    {
        self->capacity               = capacity ;
        self->centerFunctionPointers = NULL     ;
        self->entries                = NULL     ;
        self->functionCenters        = NULL     ;
        self->isOwner                = False    ;
        if ( capacity > 0 )
        {
            self->entries = Memory_AllocateArrayOfReferences ( capacity, GaussianBasis ) ;
            if ( self->entries == NULL ) GaussianBasisContainer_Deallocate ( &self ) ;
            else
            {
                auto Integer  i ;
                for ( i = 0 ; i < capacity ; i++ ) self->entries[i] = NULL ;
            }
        }
    }
    if ( self == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning (without index arrays).
!---------------------------------------------------------------------------------------------------------------------------------*/
GaussianBasisContainer *GaussianBasisContainer_Clone ( const GaussianBasisContainer *self, Status *status )
{
    GaussianBasisContainer *clone = NULL ;
    if ( self != NULL )
    {
        clone = GaussianBasisContainer_Allocate ( self->capacity, status ) ;
        if ( clone != NULL )
        {
            auto Integer  i ;
            clone->isOwner = self->isOwner ;
            if ( self->isOwner )
            {
                for ( i = 0 ; i < self->capacity ; i++ )
                {
                    clone->entries[i] = GaussianBasis_Clone ( self->entries[i], status ) ;
                    if ( clone->entries[i] == NULL ) { GaussianBasisContainer_Deallocate ( &clone ) ; break ; }
                }
            }
            else
            {
                for ( i = 0 ; i < self->capacity ; i++ ) { clone->entries[i] = self->entries[i] ; }
            }
        }
        if ( clone == NULL ) Status_Set ( status, Status_OutOfMemory ) ;
    }
    return clone ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainer_Deallocate ( GaussianBasisContainer **self )
{
    if ( (*self) != NULL )
    {
        if ( (*self)->isOwner ) /* . Always False for moment. */
        {
            auto Integer  i ;
            for ( i = 0 ; i < (*self)->capacity ; i++ ) GaussianBasis_Deallocate ( &((*self)->entries[i]) ) ;
        }
        IntegerArray1D_Deallocate ( &((*self)->centerFunctionPointers) ) ;
        IntegerArray1D_Deallocate ( &((*self)->functionCenters       ) ) ;
        Memory_Deallocate         (   (*self)->entries                 ) ;
        Memory_Deallocate         (   (*self)                          ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Largest basis.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer GaussianBasisContainer_LargestBasis ( const GaussianBasisContainer *self, const Boolean forC )
{
    Integer  n = 0 ;
    if ( self != NULL )
    {
        auto Integer        i     ;
        auto GaussianBasis *basis ;
        for ( i = 0 ; i < self->capacity ; i++ )
        {
            basis = self->entries[i] ;
            if ( forC ) n = Maximum ( n, basis->nCBF   ) ;
            else        n = Maximum ( n, basis->nBasis ) ;
        }
    }
    return n ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Largest shell.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer GaussianBasisContainer_LargestShell ( const GaussianBasisContainer *self, const Boolean forC )
{
    Integer  n = 0 ;
    if ( self != NULL )
    {
        auto Integer        i     ;
        auto GaussianBasis *basis ;
        for ( i = 0 ; i < self->capacity ; i++ )
        {
            basis = self->entries[i] ;
            n     = Maximum ( n, GaussianBasis_LargestShell ( basis, forC ) ) ;
        }
    }
    return n ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make C->S transformation for the full container. For debugging only!
! . The transformation is always made even if it's the identity!
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifdef _MakeC2S_
void GaussianBasisContainer_MakeC2S ( const GaussianBasisContainer *self, RealArray2D *T, Status *status )
{
    if ( ( self != NULL ) && ( T != NULL ) && Status_IsOK ( status ) )
    {
        if ( ( View2D_Rows    ( T ) >= GaussianBasisContainer_NumberOfWorkFunctions ( self ) ) &&
             ( View2D_Columns ( T ) >= GaussianBasisContainer_NumberOfFunctions     ( self ) ) )
        {
            auto Integer        c, C, C0 = 0, b, r, R, R0 = 0, s ;
            auto GaussianBasis *basis ;
            auto RealArray2D   *c2s   ;
            RealArray2D_Set ( T, 0.0e+00 ) ;
            /* . Loop over bases and shells. */
            for ( b = 0 ; b < self->capacity ; b++ )
            {
                basis = self->entries[b] ;
                for ( s = 0 ; s < basis->nShells ; s++ )
                {
                    R   = basis->shells[s].nCBF   ;
                    C   = basis->shells[s].nBasis ;
                    c2s = basis->shells[s].c2s    ;
                    if ( c2s == NULL )
                    {
                        for ( c = 0 ; c < C ; c++ ) Array2D_Item ( T, c+R0, c+C0 ) = 1.0e+00 ; /* . C = R for Cartesian shells but C0 != R0 for mixed containers. */
                    }
                    else
                    {
                        for ( r = 0 ; r < R ; r++ )
                        {
                            for ( c = 0 ; c < C ; c++ ) { Array2D_Item ( T, r+R0, c+C0 ) = Array2D_Item ( c2s, r, c ) ; }
                        }
                    }
                    C0 += C ;
                    R0 += R ;
                }
            }

        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}
# endif

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make index arrays.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------------------------------------------------------------
!   The arrays are:
!   - centerFunctionPointers pointers to the functions on each center such that for center c the first and
!                            last functions are cFP[c] and cFP[c+1]-1, respectively.
!   - functionCenters        the center to which the function belongs
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainer_MakeIndexArrays ( GaussianBasisContainer *self                   ,
                                              IntegerArray1D         *centerFunctionPointers ,
                                              IntegerArray1D         *functionCenters        ,
                                              Status                 *status                 )
{
    if ( ( self                   != NULL ) &&
         ( centerFunctionPointers != NULL ) &&
         ( functionCenters        != NULL ) &&
         Status_IsOK ( status ) )
    {
        if ( ( View1D_Extent ( centerFunctionPointers ) == ( self->capacity + 1 ) ) &&
             ( View1D_Extent ( functionCenters        ) ==  GaussianBasisContainer_NumberOfFunctions ( self ) ) )
        {
            auto Integer        c, f, n, t ;
            auto GaussianBasis *cBasis ;
            Array1D_Item ( centerFunctionPointers, 0 ) = 0 ;
            for ( c = t = 0 ; c < self->capacity ; c++ )
            {
                cBasis = self->entries[c] ;
                n      = cBasis->nBasis ;
                for ( f = 0 ; f < n ; f++, t++ ) Array1D_Item ( functionCenters, t ) = c ;
                Array1D_Item ( centerFunctionPointers, c+1 ) = t ;
            }
            self->centerFunctionPointers = IntegerArray1D_CloneShallow ( centerFunctionPointers, status ) ;
            self->functionCenters        = IntegerArray1D_CloneShallow ( functionCenters       , status ) ;
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Number of functions in the container.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer GaussianBasisContainer_NumberOfFunctions ( const GaussianBasisContainer *self )
{
    Integer t = 0 ;
    if ( self != NULL )
    {
        auto Integer        c ;
        auto GaussianBasis *cBasis ;
        for ( c = 0 ; c < self->capacity ; c++ )
        {
            cBasis = self->entries[c] ;
            t     += cBasis->nBasis ;
        }
    }
    return t ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Number of work functions in the container.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer GaussianBasisContainer_NumberOfWorkFunctions ( const GaussianBasisContainer *self )
{
    Integer t = 0 ;
    if ( self != NULL )
    {
        auto Integer        c ;
        auto GaussianBasis *cBasis ;
        for ( c = 0 ; c < self->capacity ; c++ )
        {
            cBasis = self->entries[c] ;
            t     += cBasis->nCBF ;
        }
    }
    return t ;
}
