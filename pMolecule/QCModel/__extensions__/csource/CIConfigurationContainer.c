/*==================================================================================================================================
! . A CI configuration container.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "CIConfigurationContainer.h"
# include "Memory.h"
# include "NumericalMacros.h"

/*==================================================================================================================================
! . Local procedures declarations.
!=================================================================================================================================*/
static IntegerArray2D  *CISetup_MakePermutations ( const Integer  m, const Integer  n, Status *status ) ;
static void             CISetup_MakeSPQR         ( CIConfigurationContainer *self, Status *status ) ;
static Integer          Factorial                ( const Integer  n ) ;

static void CIConfiguration_AllocateAlphasBetas ( CIConfiguration *self, const Integer  nActive, Status *status ) ;
static void CIConfiguration_Deallocate          ( CIConfiguration *self ) ;
static void CIConfiguration_Initialize          ( CIConfiguration *self ) ;

/*==================================================================================================================================
! . CI configuration container procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation - basic.
!---------------------------------------------------------------------------------------------------------------------------------*/
CIConfigurationContainer *CIConfigurationContainer_Allocate ( const Integer nActive, const Integer nConfigurations, Status *status )
{
    CIConfigurationContainer *self        = NULL      ;
    Status                    localStatus = Status_OK ;
    if ( ( nActive < 0 ) || ( nConfigurations < 0 ) ) Status_Set ( &localStatus, Status_InvalidArgument ) ;
    else
    {
        self = Memory_AllocateType ( CIConfigurationContainer ) ;
        if ( self != NULL )
        {
            self->nActive         = nActive         ;
            self->nConfigurations = nConfigurations ;
            self->nElectrons      = 0               ;
            if ( nConfigurations > 0 )
            {
                self->configurations = Memory_AllocateArrayOfTypes ( nConfigurations, CIConfiguration ) ;
                if ( self->configurations == NULL ) Status_Set ( &localStatus, Status_OutOfMemory ) ;
                else
                {
                    auto Integer  i ;
                    for ( i = 0 ; i < nConfigurations ; i++ )
                    {
                        CIConfiguration_Initialize          ( &(self->configurations[i]) ) ;
                    }
                    for ( i = 0 ; i < nConfigurations ; i++ )
                    {
                        CIConfiguration_AllocateAlphasBetas ( &(self->configurations[i]), self->nActive, &localStatus ) ;
                        if ( ! Status_IsValueOK ( localStatus ) ) break ;
                    }
                }
            }
        }
    }
    if ( ! Status_IsValueOK ( localStatus ) ) CIConfigurationContainer_Deallocate ( &self ) ;
    Status_Set ( status, localStatus ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
CIConfigurationContainer *CIConfigurationContainer_Clone ( const CIConfigurationContainer *self, Status *status )
{
    CIConfigurationContainer *new = NULL ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        new = CIConfigurationContainer_Allocate ( self->nActive, self->nConfigurations, status ) ;
        if ( new != NULL )
        {
            auto Integer  i ;
            new->nElectrons = self->nElectrons ;
            for ( i = 0 ; i < self->nConfigurations ; i++ )
            {
                IntegerArray1D_CopyTo ( self->configurations[i].alphas, new->configurations[i].alphas, status ) ;
                IntegerArray1D_CopyTo ( self->configurations[i].betas , new->configurations[i].betas , status ) ;
            }
            CISetup_MakeSPQR ( new, status ) ;
        }
        if ( ! Status_IsOK ( status ) ) CIConfigurationContainer_Deallocate ( &new ) ;
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CIConfigurationContainer_Deallocate ( CIConfigurationContainer **self )
{
    if ( (*self) != NULL )
    {
        if ( (*self)->configurations != NULL )
        {
            auto Integer  i ;
            for ( i = 0 ; i < (*self)->nConfigurations ; i++ ) CIConfiguration_Deallocate ( &((*self)->configurations[i]) ) ;
            Memory_Deallocate ( (*self)->configurations ) ;
        }
        Memory_Deallocate ( (*self) ) ;
        (*self) = NULL   ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get an estimate of the sparsity of the CI matrix.
! . "nonZero"  is useful for sparse matrix allocation.
! . "sparsity" is an underestimate as some of the "non-zero" elements could, in fact, be zero.
!---------------------------------------------------------------------------------------------------------------------------------*/
void CIConfigurationContainer_GetCIMatrixSparsity ( CIConfigurationContainer *self, Integer  *nonZero, Real *sparsity )
{
    if ( nonZero  != NULL ) (*nonZero ) = 0       ;
    if ( sparsity != NULL ) (*sparsity) = 0.0e+00 ;
    if ( self != NULL )
    {
        auto Integer          i, j, k, na, nai, naj, nActive, nb, nOff = 0 ;
        auto IntegerArray1D  *iAlphas, *iBetas, *jAlphas, *jBetas ;
        nActive = self->nActive ;
        for ( i = 0 ; i < self->nConfigurations ; i++ )
        {
            nai     = self->configurations[i].nAlphas ;
            iAlphas = self->configurations[i].alphas  ;
            iBetas  = self->configurations[i].betas   ;
            for ( j = 0 ; j < i ; j++ )
            {
                naj     = self->configurations[j].nAlphas ;
                jAlphas = self->configurations[j].alphas  ;
                jBetas  = self->configurations[j].betas   ;
                if ( nai != naj ) continue ;
                for ( k = na = nb = 0 ; k < nActive ; k++ )
                {
                    na += abs ( Array1D_Item ( iAlphas, k ) - Array1D_Item ( jAlphas, k ) ) ;
                    nb += abs ( Array1D_Item ( iBetas,  k ) - Array1D_Item ( jBetas,  k ) ) ;
                }
                if ( ( na + nb ) <= 4 ) nOff++ ;
            }
        }
        if ( nonZero  != NULL ) (*nonZero ) = ( self->nConfigurations + nOff ) ;
        if ( sparsity != NULL ) (*sparsity) = 100.0e+00 * ( 1.0e+00 - ( ( Real ) ( self->nConfigurations + 2 * nOff ) ) / ( ( Real ) ( self->nConfigurations * self->nConfigurations ) ) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Generate all configurations consistent with a given number of up and down electrons.
!---------------------------------------------------------------------------------------------------------------------------------*/
CIConfigurationContainer *CIConfigurationContainer_MakeFull ( const Integer nActive, const Integer nUp, const Integer nDown, Status *status )
{
    CIConfigurationContainer *self = NULL ;
    IntegerArray2D           *aPermutations = NULL, *bPermutations = NULL ;
    /* . Create the permutations. */
    aPermutations = CISetup_MakePermutations ( nUp  , nActive, status ) ;
    bPermutations = CISetup_MakePermutations ( nDown, nActive, status ) ;
    if ( ( aPermutations != NULL ) && ( bPermutations != NULL ) )
    {
        auto Integer a, b, n, na, nb, nConfigurations ;
        /* . Set some counters. */
        na = View2D_Rows ( aPermutations ) ;
        nb = View2D_Rows ( bPermutations ) ;
        nConfigurations = ( na * nb ) ;
        /* . Set up the configurations. */
        self = CIConfigurationContainer_Allocate ( nActive, nConfigurations, status ) ;
        if ( self != NULL )
        {
            auto IntegerArray1D aRow, bRow ;
            self->nElectrons = ( nUp + nDown ) ;
            for ( a = n = 0 ; a < na ; a++ )
            {
                IntegerArray2D_RowView ( aPermutations, a, False, &aRow, status ) ;
                for ( b = 0 ; b < nb ; b++, n++ )
                {
                    IntegerArray2D_RowView ( bPermutations, b, False, &bRow, status ) ;
                    IntegerArray1D_CopyTo  ( &aRow, self->configurations[n].alphas, status ) ;
                    IntegerArray1D_CopyTo  ( &bRow, self->configurations[n].betas , status ) ;
                }
            }
            /* . Make remaining configuration data.*/
            CISetup_MakeSPQR ( self, status ) ;
        }
        else CIConfigurationContainer_Deallocate ( &self ) ;
    }
    /* . Finish up. */
    IntegerArray2D_Deallocate ( &aPermutations ) ;
    IntegerArray2D_Deallocate ( &bPermutations ) ;
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Generate a combination of singles and doubles configurations.
!-----------------------------------------------------------------------------------------------------------------------------------
!
! . Easier to treat in terms of numbers of alpha and betas. Use A = (C+O), B = C, VA = V, VB = (O+V).
!
! . Singles - all:
!
!   A changes: A * VA + A * VB
!   B changes: B * VB + B * VA
!
! . Singles - preserve numbers of alpha and beta:
!
!   A changes: A * VA
!   B changes: B * VB
!
! . Doubles - all:
!
!   2A  -> 2A : A(A-1)/2 * VA*(VA-1)/2
!   2B  -> 2B : B(B-1)/2 * VB*(VB-1)/2
!   A,B -> A,B: A * B * VA * VB
!   2A  -> 2B : A(A-1)/2 * VB*(VB-1)/2
!   2B  -> 2A : B(B-1)/2 * VA*(VA-1)/2
!   2A  -> A,B: A(A-1)/2 * VA*VB
!   2B  -> A,B: B(B-1)/2 * VA*VB
!   A,B -> 2A : A * B * VA(VA-1)/2
!   A.B -> 2B : A * B * VB(VB-1)/2
!
! . Doubles - preserve numbers of alpha and beta:
!
!   2A  -> 2A : A(A-1)/2 * VA*(VA-1)/2
!   2B  -> 2B : B(B-1)/2 * VB*(VB-1)/2
!   A,B -> A,B: A * B * VA * VB
!
! . Extra configurations are needed for open shell cases to ensure that the spin wavefunction is correct.
!
!---------------------------------------------------------------------------------------------------------------------------------*/
CIConfigurationContainer *CIConfigurationContainer_MakeSinglesDoubles ( const Boolean  doSingles ,
                                                                        const Boolean  doDoubles ,
                                                                        const Integer  nActive   ,
                                                                        const Integer  nclosed   ,
                                                                        const Integer  nopen     ,
                                                                              Status  *status    )
{
    auto Integer                   na, naa, nab, nb, nbb, nConfigurations, nvirtual, va, vaa, vab, vb, vbb ;
    auto CIConfigurationContainer *self = NULL ;
    /* . Set up some counters. */
    nvirtual = nActive - ( nclosed + nopen ) ;
    na  = nclosed + nopen ;
    nb  = nclosed ;
    va  = nActive - na ;
    vb  = nActive - nb ;
    naa = ( na * ( na - 1 ) ) / 2 ;
    nab = na * nb ;
    nbb = ( nb * ( nb - 1 ) ) / 2 ;
    vaa = ( va * ( va - 1 ) ) / 2 ;
    vab = va * vb ;
    vbb = ( vb * ( vb - 1 ) ) / 2 ;
    /* . Determine the number of configurations - including the ground state! */
    nConfigurations = 1 ;
    if ( doSingles ) nConfigurations += ( na * va + nb * vb ) ;
    if ( doDoubles ) nConfigurations += ( naa * vaa + nbb * vbb + nab * vab ) ;
    /* . Extra configurations needed for open-shell systems. */
    if ( nopen > 0 )
    {
        if ( doSingles && ( ! doDoubles ) ) nConfigurations += ( nopen * nclosed * nvirtual ) ;
        if ( doDoubles )
        {
            nConfigurations += nopen * ( 4 * nbb * vaa + nopen * ( nbb * nvirtual + vaa * nclosed ) ) + ( nopen * ( nopen - 1 ) ) / 2 * nbb * vaa ;
            if ( ! doSingles ) nConfigurations += ( 2 * nclosed * nvirtual ) ;
        }
    }
    /* . Allocate space. */
    self = CIConfigurationContainer_Allocate ( nActive, nConfigurations, status ) ;
    /* . Set up the configurations. */
    if ( self != NULL )
    {
        auto Integer         i, j, n = 1, o, p, u, v ;
        auto IntegerArray1D  view ;
        /* . Initialize all configurations to the ground state. */
        self->nElectrons = ( na + nb ) ;
        for ( i = 0 ; i < nConfigurations ; i++ )
        {
            IntegerArray1D_Set ( self->configurations[i].alphas, 0 ) ; IntegerArray1D_View ( self->configurations[i].alphas, 0, na, 1, False, &view, NULL ) ; IntegerArray1D_Set ( &view, 1 ) ;
            IntegerArray1D_Set ( self->configurations[i].betas , 0 ) ; IntegerArray1D_View ( self->configurations[i].betas , 0, nb, 1, False, &view, NULL ) ; IntegerArray1D_Set ( &view, 1 ) ;
        }
        /* . Singles. */
        if ( doSingles )
        {
            for ( i = 0 ; i < na ; i++ )
            {
                for ( u = na ; u < nActive ; n++, u++ ) { Array1D_Item ( self->configurations[n].alphas, i ) = 0 ; Array1D_Item ( self->configurations[n].alphas, u ) = 1 ; }
            }
            for ( i = 0 ; i < nb ; i++ )
            {
                for ( u = nb ; u < nActive ; n++, u++ ) { Array1D_Item ( self->configurations[n].betas, i ) = 0 ; Array1D_Item ( self->configurations[n].betas, u ) = 1 ; }
            }
            /* . Add in the extra doubles that are necessary to have a consistent spin wavefunction. */
            if ( ! doDoubles && ( nopen > 0 ) )
            {
                for ( i = 0 ; i < nb ; i++ )
                {
                    for ( o = nb ; o < na ; o++ )
                    {
                        for ( u = na ; u < nActive ; n++, u++ ) { Array1D_Item ( self->configurations[n].betas, i ) = 0 ; Array1D_Item ( self->configurations[n].alphas, o ) = 0 ;
                                                                  Array1D_Item ( self->configurations[n].betas, o ) = 1 ; Array1D_Item ( self->configurations[n].alphas, u ) = 1 ; }
                    }
                }
            }
        }
        /* . Doubles. */
        if ( doDoubles )
        {
            for ( i = 1 ; i < na ; i++ )
            {
                for ( j = 0 ; j < i ; j++ )
                {
                    for ( u = (na+1) ; u < nActive ; u++ )
                    {
                        for ( v = na ; v < u ; n++, v++ ) { Array1D_Item ( self->configurations[n].alphas, i ) = 0 ; Array1D_Item ( self->configurations[n].alphas, j ) = 0 ;
                                                            Array1D_Item ( self->configurations[n].alphas, u ) = 1 ; Array1D_Item ( self->configurations[n].alphas, v ) = 1 ; }
                    }
                }
            }
            for ( i = 1 ; i < nb ; i++ )
            {
                for ( j = 0 ; j < i ; j++ )
                {
                    for ( u = (nb+1) ; u < nActive ; u++ )
                    {
                        for ( v = nb ; v < u ; n++, v++ ) { Array1D_Item ( self->configurations[n].betas, i ) = 0 ; Array1D_Item ( self->configurations[n].betas, j ) = 0 ;
                                                            Array1D_Item ( self->configurations[n].betas, u ) = 1 ; Array1D_Item ( self->configurations[n].betas, v ) = 1 ; }
                    }
                }
            }
            for ( i = 0 ; i < na ; i++ )
            {
                for ( j = 0 ; j < nb ; j++ )
                {
                    for ( u = na ; u < nActive ; u++ )
                    {
                        for ( v = nb ; v < nActive ; n++, v++ ) { Array1D_Item ( self->configurations[n].alphas, i ) = 0 ; Array1D_Item ( self->configurations[n].betas, j ) = 0 ;
                                                                  Array1D_Item ( self->configurations[n].alphas, u ) = 1 ; Array1D_Item ( self->configurations[n].betas, v ) = 1 ; }
                    }
                }
            }
            /* . Add in the extra configurations that are necessary to have a consistent spin wavefunction. */
            if ( nopen > 0 )
            {
                for ( i = 1 ; i < nb ; i++ )
                {
                    for ( j = 0 ; j < i ; j++ )
                    {
                        for ( u = (na+1) ; u < nActive ; u++ )
                        {
                            for ( v = na ; v < u ; v++ )
                            {
                                for ( o = nb ; o < na ; n += 4, o++ )
                                {
                                     /* . 2 closed alpha -> 2 virtual alpha. */
                                     Array1D_Item ( self->configurations[n  ].alphas, i ) = 0 ; Array1D_Item ( self->configurations[n  ].betas , j ) = 0 ;
                                     Array1D_Item ( self->configurations[n  ].alphas, o ) = 0 ; Array1D_Item ( self->configurations[n  ].betas , o ) = 1 ;
                                     Array1D_Item ( self->configurations[n  ].alphas, u ) = 1 ; Array1D_Item ( self->configurations[n  ].alphas, v ) = 1 ;
                                     Array1D_Item ( self->configurations[n+1].alphas, j ) = 0 ; Array1D_Item ( self->configurations[n+1].betas , i ) = 0 ;
                                     Array1D_Item ( self->configurations[n+1].alphas, o ) = 0 ; Array1D_Item ( self->configurations[n+1].betas , o ) = 1 ;
                                     Array1D_Item ( self->configurations[n+1].alphas, u ) = 1 ; Array1D_Item ( self->configurations[n+1].alphas, v ) = 1 ;
                                     /* . 2 closed beta -> 2 virtual beta. */
                                     Array1D_Item ( self->configurations[n+2].betas , i ) = 0 ; Array1D_Item ( self->configurations[n+2].betas , j ) = 0 ;
                                     Array1D_Item ( self->configurations[n+2].alphas, o ) = 0 ; Array1D_Item ( self->configurations[n+2].betas , o ) = 1 ;
                                     Array1D_Item ( self->configurations[n+2].alphas, u ) = 1 ; Array1D_Item ( self->configurations[n+2].betas , v ) = 1 ;
                                     Array1D_Item ( self->configurations[n+3].betas , i ) = 0 ; Array1D_Item ( self->configurations[n+3].betas , j ) = 0 ;
                                     Array1D_Item ( self->configurations[n+3].alphas, o ) = 0 ; Array1D_Item ( self->configurations[n+3].betas , o ) = 1 ;
                                     Array1D_Item ( self->configurations[n+3].alphas, v ) = 1 ; Array1D_Item ( self->configurations[n+3].betas , u ) = 1 ;
                                }
                            }
                        }
                    }
                }
                for ( o = nb ; o < na ; o++ )
                {
                    for ( p = nb ; p < na ; p++ )
                    {
                        if ( p != o )
                        {
                            /* . 1 closed alpha, 1 open alpha -> 2 virtual alpha. */
                            for ( i = 0 ; i < nb ; i++ )
                            {
                                for ( u = (na+1) ; u < nActive ; u++ )
                                {
                                    for ( v = na ; v < u ; n++, v++ )
                                    {
                                         Array1D_Item ( self->configurations[n].betas , i ) = 0 ; Array1D_Item ( self->configurations[n].alphas, o ) = 0 ;
                                         Array1D_Item ( self->configurations[n].alphas, p ) = 0 ; Array1D_Item ( self->configurations[n].betas , p ) = 1 ;
                                         Array1D_Item ( self->configurations[n].alphas, u ) = 1 ; Array1D_Item ( self->configurations[n].alphas, v ) = 1 ;
                                    }
                                }
                            }
                            /* . 2 closed beta -> 1 open beta, 1 virtual beta (or 1 open alpha, 1 closed beta -> 1 virtual alpha, 1 virtual beta). */
                            for ( i = 1 ; i < nb ; i++ )
                            {
                                for ( j = 0 ; j < i ; j++ )
                                {
                                    for ( u = na ; u < nActive ; n++, u++ )
                                    {
                                         Array1D_Item ( self->configurations[n].betas , i ) = 0 ; Array1D_Item ( self->configurations[n].betas , j ) = 0 ;
                                         Array1D_Item ( self->configurations[n].betas , o ) = 1 ; Array1D_Item ( self->configurations[n].alphas, p ) = 0 ;
                                         Array1D_Item ( self->configurations[n].betas , p ) = 1 ; Array1D_Item ( self->configurations[n].alphas, u ) = 1 ;
                                    }
                                }
                            }
                        }
                    }
                }
                for ( o = nb ; o < na ; o++ )
                {
                    /* . 1 closed alpha, 1 closed beta (same closed) -> 1 virtual alpha, 1 virtual beta (different virtual). */
                    for ( i = 0 ; i < nb ; i++ )
                    {
                        for ( u = (na+1) ; u < nActive ; u++ )
                        {
                            for ( v = na ; v < u ; n++, v++ )
                            {
                                 Array1D_Item ( self->configurations[n].alphas, i ) = 0 ; Array1D_Item ( self->configurations[n].betas , i ) = 0 ;
                                 Array1D_Item ( self->configurations[n].alphas, o ) = 0 ; Array1D_Item ( self->configurations[n].betas , o ) = 1 ;
                                 Array1D_Item ( self->configurations[n].alphas, u ) = 1 ; Array1D_Item ( self->configurations[n].alphas, v ) = 1 ;
                            }
                        }
                    }
                    /* . 1 closed alpha, 1 closed beta (different closed) -> 1 virtual alpha, 1 virtual beta (same virtual). */
                    for ( i = 1 ; i < nb ; i++ )
                    {
                        for ( j = 0 ; j < i ; j++ )
                        {
                            for ( u = na ; u < nActive ; n++, u++ )
                            {
                                 Array1D_Item ( self->configurations[n].betas , i ) = 0 ; Array1D_Item ( self->configurations[n].betas , j ) = 0 ;
                                 Array1D_Item ( self->configurations[n].alphas, o ) = 0 ; Array1D_Item ( self->configurations[n].betas , o ) = 1 ;
                                 Array1D_Item ( self->configurations[n].alphas, u ) = 1 ; Array1D_Item ( self->configurations[n].betas , u ) = 1 ;
                            }
                        }
                    }
                }
                /* . Quadruples arising from when two open shell orbitals have beta spin. */
                for ( o = (nb+1) ; o < na ; o++ )
                {
                    for ( p = nb ; p < o ; p++ )
                    {
                        for ( i = 1 ; i < nb ; i++ )
                        {
                            for ( j = 0 ; j < i ; j++ )
                            {
                                for ( u = (na+1) ; u < nActive ; u++ )
                                {
                                    for ( v = na ; v < u ; n++, v++ )
                                    {
                                         Array1D_Item ( self->configurations[n].betas , i ) = 0 ; Array1D_Item ( self->configurations[n].betas , j ) = 0 ;
                                         Array1D_Item ( self->configurations[n].alphas, o ) = 0 ; Array1D_Item ( self->configurations[n].betas , o ) = 1 ;
                                         Array1D_Item ( self->configurations[n].alphas, p ) = 0 ; Array1D_Item ( self->configurations[n].betas , p ) = 1 ;
                                         Array1D_Item ( self->configurations[n].alphas, u ) = 1 ; Array1D_Item ( self->configurations[n].alphas, v ) = 1 ;
                                    }
                                }
                            }
                        }
                    }
                }
                if ( ! doSingles )
                {
                    for ( i = 0 ; i < nb ; i++ )
                    {
                        for ( u = na ; u < nActive ; n += 2, u++ )
                        {
                            Array1D_Item ( self->configurations[n  ].alphas, i ) = 0 ; Array1D_Item ( self->configurations[n  ].alphas, u ) = 1 ;
                            Array1D_Item ( self->configurations[n+1].betas , i ) = 0 ; Array1D_Item ( self->configurations[n+1].betas , u ) = 1 ;
                        }
                    }
                }
            }
        }
        /* . Make remaining configuration data.*/
        CISetup_MakeSPQR ( self, status ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Generate a state given a set of user-specified microStates.
!---------------------------------------------------------------------------------------------------------------------------------*/
CIConfigurationContainer *CIConfigurationContainer_MakeUserSpecified ( const IntegerArray2D  *microStates     ,
                                                                       const Integer          activeOrbitals  ,
                                                                       const Integer          activeElectrons ,
                                                                             Status          *status          )
{
    CIConfigurationContainer *self = NULL ;
    if ( ( microStates != NULL ) && Status_IsOK ( status ) )
    {
        auto Boolean         isOK = True ;
        auto Integer         i, nConfigurations ;
        auto IntegerArray1D  state ;
        /* . Basic checks. */
        nConfigurations = View2D_Rows ( microStates ) ;
        for ( i = 0 ; i < nConfigurations ; i++ )
        {
            IntegerArray2D_RowView  ( microStates, i, False, &state, NULL ) ;
            if ( IntegerArray1D_Sum ( &state ) != activeElectrons ) { isOK = False ; break ; }
        }
        /* . Basic checks. */
        if ( isOK && ( View2D_Columns ( microStates ) == ( 2 * activeOrbitals ) ) )
        {
            /* . Set up the configurations. */
            self = CIConfigurationContainer_Allocate ( activeOrbitals, nConfigurations, status ) ;
            if ( self != NULL )
            {
                auto IntegerArray1D  alphas, betas ;
                self->nElectrons = activeElectrons ;
                for ( i = 0 ; i < nConfigurations ; i++ )
                {
                    IntegerArray2D_View1D ( microStates, 1, i, 0             , activeOrbitals, 1, False, &alphas, status ) ;
                    IntegerArray2D_View1D ( microStates, 1, i, activeOrbitals, activeOrbitals, 1, False, &betas , status ) ;
                    IntegerArray1D_CopyTo ( &alphas, self->configurations[i].alphas, status ) ;
                    IntegerArray1D_CopyTo ( &betas , self->configurations[i].betas , status ) ;
                }
                /* . Make remaining configuration data.*/
                CISetup_MakeSPQR ( self, status ) ;
            }
        }
        else Status_Set ( status, Status_InvalidArgument ) ;
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Numbers characterizing the container.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer  CIConfigurationContainer_NumberOfActiveElectrons ( const CIConfigurationContainer *self ) { return ( ( self == NULL ) ? 0 : self->nElectrons      ) ; }
Integer  CIConfigurationContainer_NumberOfActiveOrbitals  ( const CIConfigurationContainer *self ) { return ( ( self == NULL ) ? 0 : self->nActive         ) ; }
Integer  CIConfigurationContainer_NumberOfConfigurations  ( const CIConfigurationContainer *self ) { return ( ( self == NULL ) ? 0 : self->nConfigurations ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find the spins of a set of CI state vectors (in column major format).
!---------------------------------------------------------------------------------------------------------------------------------*/
void CIConfigurationContainer_StateSpins ( const CIConfigurationContainer *self    ,
                                           const RealArray2D              *vectors ,
                                                 RealArray1D              *spins   ,
                                                 Status                   *status  )
{
    if ( ( self    != NULL ) &&
         ( vectors != NULL ) &&
         ( spins   != NULL ) &&
         Status_IsOK ( status ) )
    {
        auto Integer  c, s ;
        c = self->nConfigurations ;
        s = View1D_Extent ( spins ) ;
        if ( ( c == View2D_Columns ( vectors ) ) &&
             ( s <= View2D_Rows    ( vectors ) ) )
        {
            auto Integer  i, j, k, m ;
            auto Real     e = ( 0.5e+00 * ( Real ) self->nElectrons ), p, sFactor ;
            RealArray1D_Set ( spins, 0.0e+00 ) ;
            for ( i = 0 ; i < s ; i++ )
            {
                sFactor = e ;
                for ( j = 0 ; j < c ; j++ )
                {
                    sFactor -= ( pow ( Array2D_Item ( vectors, i, j ), 2 ) * self->configurations[j].spin / 4.0e+00 ) ;
                    for ( k = 0 ; k < self->configurations[j].nSPQR ; k++ )
                    {
                        if ( Array1D_Item ( self->configurations[j].parity, k ) ) p = -2.0e+00 ;
                        else                                                      p =  2.0e+00 ;
                        m        = Array1D_Item ( self->configurations[j].spqr, k ) ;
                        sFactor += ( p * Array2D_Item ( vectors, i, j ) * Array2D_Item ( vectors, i, m ) ) ;
                    }
                }
                Array1D_Item ( spins, i ) = sFactor ;
            }
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
}

/*==================================================================================================================================
! . CI setup procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Generate all possible permutations (m in n) - small numbers only!
!---------------------------------------------------------------------------------------------------------------------------------*/
# define MAXIMUMFACTORIAL 12
static Integer Factorial ( const Integer  n )
{
    auto Integer f = 0, i ;
    if ( n <= MAXIMUMFACTORIAL )
    {
        f = 1 ;
        for ( i = 2 ; i <= n ; i++ ) f *= i ;
    }
    return f ;
}

static IntegerArray2D *CISetup_MakePermutations ( const Integer m, const Integer n, Status *status )
{
    IntegerArray2D *permutations = NULL ;
    if ( Status_IsOK ( status ) )
    {
        /* . Find the number of permutations. */
        auto Integer nPermutations = Factorial ( n ) / ( Factorial ( m ) * Factorial ( n - m ) ) ;
        /* . Allocate space. */
        permutations = IntegerArray2D_AllocateWithExtents ( nPermutations, n, status ) ;
        /* . Generate the permutations. */
        if ( permutations != NULL )
        {
            auto Boolean          isOK = True ;
            auto Integer          i, j, k ;
            auto IntegerArray1D  *indices = NULL ;
            /* . Initialization. */
            IntegerArray2D_Set ( permutations, 0 ) ;
            /* . Define the initial set of indices. */
            indices = IntegerArray1D_AllocateWithExtent ( m, NULL ) ;
            for ( i = 0 ; i < m ; i++ ) Array1D_Item ( indices, i ) = i ;
            /* . First combination. */
            for ( i = 0 ; i < m ; i++ ) Array2D_Item ( permutations, 0, i ) = 1 ;
            /* . Subsequent combinations. */
            for ( k = 1 ; k < nPermutations ; k++ )
            {
                i = m - 1 ;
                while ( ( i > 0 ) && ( Array1D_Item ( indices, i ) == n - m + i ) ) i -= 1 ;
                if ( ( i == 0 ) && ( Array1D_Item ( indices, i ) == n - m + i ) ) { isOK = False ; break ; }
                Array1D_Item ( indices, i ) += 1 ;
                for ( j = i ; j < m-1 ; j++ ) Array1D_Item ( indices, j+1 ) = Array1D_Item ( indices, j ) + 1 ;
                for ( i = 0 ; i < m ; i++ ) Array2D_Item ( permutations, k, Array1D_Item ( indices, i ) ) = 1 ;
            }
            if ( ! isOK )
            {
                IntegerArray2D_Deallocate ( &permutations ) ;
                Status_Set ( status, Status_AlgorithmError ) ;
            }
            /* . Finish up. */
            IntegerArray1D_Deallocate ( &indices ) ;
        }
    }
    return permutations ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the SPQR data which is necessary for the calculation of state spins.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CISetup_MakeSPQR ( CIConfigurationContainer *self, Status *status )
{
    if ( ( self != NULL ) && ( Status_IsOK ( status ) ) )
    {
        auto Integer  i, na, nab, nb ;
        /* . Make basic spin data. */
        for ( i = 0 ; i < self->nConfigurations ; i++ )
        {
            na  = IntegerArray1D_Sum ( self->configurations[i].alphas ) ;
            nb  = IntegerArray1D_Sum ( self->configurations[i].betas  ) ;
            nab = IntegerArray1D_Dot ( self->configurations[i].alphas, self->configurations[i].betas, NULL ) ;
            self->configurations[i].nAlphas = na ;
            self->configurations[i].spin    = 4.0e+00 * ( ( Real ) nab ) - pow ( ( Real ) ( na - nb ), 2 ) ;
        }
        /* . Make SPQR. */
        if ( self->nConfigurations > 1 )
        {
            auto Integer         j, k, n, nai, naj, nSPQR = 0, p, q, r, s, t ;
            auto BooleanArray1D  bView, *parity ;
            auto IntegerArray1D  cView, *iAlphas, *iBetas, *jAlphas, *jBetas, *spqr ;
            /* . Temporary space. */
            parity = BooleanArray1D_AllocateWithExtent ( self->nConfigurations - 1, status ) ;
            spqr   = IntegerArray1D_AllocateWithExtent ( self->nConfigurations - 1, status ) ;
            /* . Double loop over configurations. */
            if ( ( parity != NULL ) && ( spqr != NULL ) )
            {
                for ( i = 1 ; i < self->nConfigurations ; i++ )
                {
                    nai     = self->configurations[i].nAlphas ;
                    iAlphas = self->configurations[i].alphas  ;
                    iBetas  = self->configurations[i].betas   ;
                    nSPQR   = 0 ;
                    for ( j = 0 ; j < i ; j++ )
                    {
                        naj     = self->configurations[j].nAlphas ;
                        jAlphas = self->configurations[j].alphas  ;
                        jBetas  = self->configurations[j].betas   ;
                        /* . Skip if there are different numbers of alpha orbitals in the two configurations. */
                        if ( nai != naj ) continue ;
                        /* . Find the differences in the numbers of alpha and beta orbitals including positional information. */
                        for ( k = na = nb = 0 ; k < self->nActive ; k++ )
                        {
                            na += abs ( Array1D_Item ( iAlphas, k ) - Array1D_Item ( jAlphas, k ) ) ;
                            nb += abs ( Array1D_Item ( iBetas,  k ) - Array1D_Item ( jBetas,  k ) ) ;
                        }
                        /* . States differing by one alpha and one beta orbital. */
                        if ( ( na == 2 ) && ( nb == 2 ) )
                        {
                            /* . Find the alpha orbitals that differ (i and j with j > i). */
                            for ( p = -1, n = 0   ; n < self->nActive ; n++ ) { if ( Array1D_Item ( iAlphas, n ) != Array1D_Item ( jAlphas, n ) ) { p = n ; break ; } }
                            for ( q = -1, n = p+1 ; n < self->nActive ; n++ ) { if ( Array1D_Item ( iAlphas, n ) != Array1D_Item ( jAlphas, n ) ) { q = n ; break ; } }
                            /* . Find the beta orbitals that differ (k and l with l > k). */
                            for ( r = -1, n = 0   ; n < self->nActive ; n++ ) { if ( Array1D_Item ( iBetas , n ) != Array1D_Item ( jBetas , n ) ) { r = n ; break ; } }
                            for ( s = -1, n = r+1 ; n < self->nActive ; n++ ) { if ( Array1D_Item ( iBetas , n ) != Array1D_Item ( jBetas , n ) ) { s = n ; break ; } }
                            /* . Save configurations of type (paqb and pbqa) and calculate the parity. */
                            if ( ( p == r ) && ( q == s ) && ( Array1D_Item ( iAlphas, p ) != Array1D_Item ( iBetas, p ) ) )
                            {
                                t = -1 ;
                                for ( n = p+1 ; n < self->nActive ; n++ ) t += Array1D_Item ( iAlphas, n ) ;
                                for ( n = q+1 ; n < self->nActive ; n++ ) t += Array1D_Item ( jAlphas, n ) ;
                                for ( n = 0   ; n <= q            ; n++ ) t += Array1D_Item ( iBetas , n ) ;
                                for ( n = 0   ; n <= p            ; n++ ) t += Array1D_Item ( jBetas , n ) ;
                                Array1D_Item ( parity, nSPQR ) = IsOdd ( t ) ;
                                Array1D_Item ( spqr  , nSPQR ) = j ;
                                nSPQR++ ;
                            }
                        }
                    }
                    self->configurations[i].nSPQR = nSPQR ;
                    if ( nSPQR > 0 )
                    {
                        BooleanArray1D_View ( parity, 0, nSPQR, 1, False, &bView, NULL ) ;
                        IntegerArray1D_View ( spqr  , 0, nSPQR, 1, False, &cView, NULL ) ;
                        self->configurations[i].parity = BooleanArray1D_CloneDeep ( &bView, status ) ;
                        self->configurations[i].spqr   = IntegerArray1D_CloneDeep ( &cView, status ) ;
                        if ( ! Status_IsOK ( status ) ) break ;
                    }
                }
            }
            /* . Finish up.*/
            BooleanArray1D_Deallocate ( &parity ) ;
            IntegerArray1D_Deallocate ( &spqr   ) ;
        }
    }
}

/*==================================================================================================================================
! . CI configuration procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocate alphas and betas.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CIConfiguration_AllocateAlphasBetas ( CIConfiguration *self, const Integer  nActive, Status *status )
{
    if ( self != NULL )
    {
        self->alphas = IntegerArray1D_AllocateWithExtent ( nActive, status ) ;
        self->betas  = IntegerArray1D_AllocateWithExtent ( nActive, status ) ;
        IntegerArray1D_Set ( self->alphas, 0 ) ;
        IntegerArray1D_Set ( self->betas , 0 ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CIConfiguration_Deallocate ( CIConfiguration *self )
{
    if ( self != NULL )
    {
        BooleanArray1D_Deallocate ( &(self->parity) ) ;
        IntegerArray1D_Deallocate ( &(self->alphas) ) ;
        IntegerArray1D_Deallocate ( &(self->betas ) ) ;
        IntegerArray1D_Deallocate ( &(self->spqr  ) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Initialization.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void CIConfiguration_Initialize ( CIConfiguration *self )
{
    if ( self != NULL )
    {
        self->nAlphas = 0 ;
        self->nSPQR   = 0 ;
        self->spin    = 0.0e+00 ;
        self->alphas  = NULL ;
        self->betas   = NULL ;
        self->parity  = NULL ;
        self->spqr    = NULL ;
    }
}
