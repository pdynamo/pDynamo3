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
        self->capacity = capacity ;
        self->entries  = NULL     ;
        self->isOwner  = False    ;
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
! . Cloning.
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
                    clone->entries[i] = GaussianBasis_Clone ( self->entries[i] ) ;
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
        Memory_Deallocate ( (*self)->entries ) ;
        Memory_Deallocate ( (*self)          ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Largest basis.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer  GaussianBasisContainer_LargestBasis ( const GaussianBasisContainer *self, const Integer  doWork )
{
    Integer  n = 0 ;
    if ( self != NULL )
    {
        auto Integer        i     ;
        auto GaussianBasis *basis ;
        for ( i = 0 ; i < self->capacity ; i++ )
        {
            basis = self->entries[i] ;
            if ( doWork ) n = Maximum ( n, basis->nbasisw ) ;
            else          n = Maximum ( n, basis->nbasis  ) ;
        }
    }
    return n ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Basis atom indices.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainer_MakeBasisAtomIndices ( const GaussianBasisContainer *self    ,
                                                   const Boolean                 doWork  ,
                                                         IntegerArray1D         *indices ,
                                                         Status                 *status  )
{
    if ( ( self != NULL ) && ( indices != NULL ) && Status_IsOK ( status ) )
    {
        if ( View1D_Extent ( indices ) == GaussianBasisContainer_NumberOfBasisFunctions ( self, doWork ) )
        {
            auto Integer        c, f, n, t ;
            auto GaussianBasis *cBasis  ;
            IntegerArray1D_Set ( indices, 0 ) ;
            for ( c = t = 0 ; c < self->capacity ; c++ )
            {
                cBasis = self->entries[c] ;
                if ( doWork ) n = cBasis->nbasisw ;
                else          n = cBasis->nbasis  ;
                for ( f = 0 ; f < n ; f++, t++ ) Array1D_Item ( indices, t ) = c ;
            }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Basis indices.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainer_MakeBasisIndices ( const GaussianBasisContainer *self    ,
                                               const Boolean                 doWork  ,
                                                     IntegerArray1D         *indices ,
                                                     Status                 *status  )
{
    if ( ( self != NULL ) && ( indices != NULL ) && Status_IsOK ( status ) )
    {
        if ( View1D_Extent ( indices ) == ( self->capacity + 1 ) )
        {
            auto Integer        i, n = 0 ;
            auto GaussianBasis *iBasis   ;
            IntegerArray1D_Set ( indices, 0 ) ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                Array1D_Item ( indices, i ) = n ;
                iBasis = self->entries[i] ;
                if ( doWork ) n += iBasis->nbasisw ;
                else          n += iBasis->nbasis  ;
            }
            Array1D_Item ( indices, self->capacity ) = n ;
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the transformations from the working to actual representations.
! . These transformations are block diagonal and, hence, sparse.
! . However, for ease, they are treated as dense here.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasisContainer_MakeFunctionTransformations ( const GaussianBasisContainer *self   ,
                                                                RealArray2D            *c2o    ,
                                                                RealArray2D            *o2c    ,
                                                                Status                 *status )
{
    if ( ( self != NULL ) &&
         ( c2o  != NULL ) &&
         ( o2c  != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer  nC, nS ;
        nC = GaussianBasisContainer_NumberOfBasisFunctions ( self, True  ) ;
        nS = GaussianBasisContainer_NumberOfBasisFunctions ( self, False ) ;
        if ( ( View2D_Rows ( c2o ) == nC ) && ( View2D_Columns ( c2o ) == nS ) &&  
             ( View2D_Rows ( o2c ) == nC ) && ( View2D_Columns ( o2c ) == nS ) )
        {
            auto Integer        c, c0 = 0, i, s, s0 = 0 ;
            auto GaussianBasis *iBasis ;
            auto RealArray2D   *m ;
            RealArray2D_Set ( c2o, 0.0e+00 ) ;
            RealArray2D_Set ( o2c, 0.0e+00 ) ;
            for ( i = 0 ; i < self->capacity ; i++ )
            {
                iBasis = self->entries[i] ;
                m      = iBasis->c2o ;
                for ( c = 0 ; c < View2D_Rows ( m ) ; c++ )
                {
                    for ( s = 0 ; s < View2D_Columns ( m ) ; s++ ) Array2D_Item ( c2o, c + c0, s + s0 ) =  Array2D_Item ( m, c, s ) ;
                }
                m = iBasis->o2c ;
                for ( c = 0 ; c < View2D_Rows ( m ) ; c++ )
                {
                    for ( s = 0 ; s < View2D_Columns ( m ) ; s++ ) Array2D_Item ( o2c, c + c0, s + s0 ) =  Array2D_Item ( m, c, s ) ;
                }
                c0 += iBasis->nbasisw ;
                s0 += iBasis->nbasis  ;
            }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Number of basis functions in the container.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer GaussianBasisContainer_NumberOfBasisFunctions ( const GaussianBasisContainer *self   ,
                                                        const Boolean                 doWork )
{
    Integer t = 0 ;
    if ( self != NULL )
    {
        auto Integer        c ;
        auto GaussianBasis *cBasis ;
        for ( c = 0 ; c < self->capacity ; c++ )
        {
            cBasis = self->entries[c] ;
            if ( doWork ) t += cBasis->nbasisw ;
            else          t += cBasis->nbasis  ;
        }
    }
    return t ;
}
