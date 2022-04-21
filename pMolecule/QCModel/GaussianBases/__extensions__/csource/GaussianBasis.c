/*==================================================================================================================================
! . This module handles Gaussian basis functions
!===================================================================================================================================
!
! . Notes:
!
!   1. Basis sets may be Cartesian or spherical harmonical.
!
!   2. Shells can be of a single angular momentum type or have multiple
!      consecutive values (from lLow to lHigh). E.g. s, p, d, f, g, sp, spd
!      spdf, df, etc.
!
!   3. The Cartesian basis function order is calculated as follows:
!
!      for l in range ( lmin, lmax + 1 ):
!          for z in range ( 0, l + 1 ):
!              for y in range ( 0, l - z + 1 ):
!                  x = l - y - z
!
!   4. The spherical harmonic basis function order for each l is:
!
!      m = 0, 1, -1, 2, -2, 3, -3, ... = 0, 1c, 1s, 2c, 2s, 3c, 3s, ...
!
!   5. Basis functions are a sum of primitive functions multiplied by an appropriate angular function:
!
!          f_i = [ sum_p c_p exp(-a_p*r^2) ] * A_i
!
!      All functions of a given shell have the same primitive expansion and only the angular functions
!      differ.
!
!      For Cartesian functions and spherical harmonic functions, the angular functions are, respectively:
!
!          A_i = x^(u_i) * y^(v_i) * z^(w_i)
!          A_i = r^l * RealSphericalHarmonic (independent of r and assumed normalized)
!
!          l is the angular momentum for the shell with l = u_i + v_i + w_i (u,v,w vary but l constant).
!
!   6. The (overlap) normalizations for Cartesian and spherical harmonic primitives are, respectively:
!
!          Nc_ip = sqrt[ ( 2 a_p / pi )^(3/2) * ( ( 4 a_p )^l / ( (2u_i-1)!! * (2v_i-1)!! * (2w_i-1)!! ) ) ]
!          Ns_p  = sqrt[ ( 2 a_p / pi )^(3/2) * ( ( 4 a_p )^l / (2l-1)!! ) ]
!
!          These two expressions are equivalent when one of the Cartesian powers = l and the others are 0.
!
!      By convention the primitives of input basis sets are often spherically harmonically normalized. To get
!      unnormalized primitives for actual computation the input coefficients must be modified by including the
!      normalization factor using: c_p(unnorm) = c_p(norm) * Ns_p. This follows from:
!
!          f_i = [ sum_p c_p(norm) * ( Ns_p * e_p ) ] * A_i = [ sum_p c_p(unnorm) * e_p ] * A_i 
!
!   7. The normalization, N_i, of a function f_i is independent of other functions as long as all coefficients
!      for that function are modified consistently, i.e. c_p for f_i -> c_p * N_i = c_ip. The final results
!      (orbitals, energies, etc.) should also be independent of normalization as the weights of each function
!      in these results are determined
!
!   8. Cartesian can be transformed into spherical harmonic functions and vice-versa. As the transformations
!      mix functions with different angular components, the transformations must correspond to the angular
!      normalization of the functions that is being used.
!
!      The transformations calculated here are those between normalized Cartesian and normalized spherical
!      harmonic functions. These transformations can be changed into those for unnormalized Cartesians
!      using:
!
!      T_c_to_s_unnorm = T_c_to_s_norm * m_i
!      T_s_to_c_unnorm = T_s_to_c_norm / m_i
!
!      where m_i = sqrt ( ( 2l-1 )!! / ( 2u_i-1 )!! ( 2v_i-1 )!! ( 2w_i-1 )!! ).
!
!      Alternatively, these m_i factors can be incorporated into the primitive coefficients for the function
!      using c_ip = c_p * m_i.
!
!   9. The Cartesian to spherical transformation procedures produce transformations between normalized Cartesian
!      and normalized spherical harmonic functions. The formulas are from:
!
!          H. B. Schlegel and M. J. Frisch
!          "Transformation Between Cartesian and Pure Spherical Harmonic Gaussians"
!          Int J Quant Chem 54, 83-87 (1995).
!
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "GaussianBasis.h"
# include "Memory.h"
# include "NumericalMacros.h"

/* . Options - one of either _NormalizePrimitives_ or _UnnormalizeTransformations_ must be set. */
# define _NormalizePrimitives_
/*# define _PrintC2S_*/
/*# define _UnnormalizeTransformations_*/

/*----------------------------------------------------------------------------------------------------------------------------------
! . Gauss-Hermite quadrature parameters (used for KE and overlap integrals).
! . If necessary extra weights are easy to calculate or to find.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The indices into the data. */
/* . Just use GHINDEX[GHMAXPT+1] and then for p in GHINDEX[gh], GHINDEX[gh+1]. */
const Integer  GHINDEX[GHMAXPT+1] = { 0, 1, 3, 6, 10, 15, 21, 28, 36, 45, 55 } ;

/* . The abscissae. */
const Real GHABSCISSAE[GHNDATA] =
    { 0.0e+00,
     -0.707106781186548e+00,  0.707106781186548e+00,
     -1.224744871391589e+00,  0.0e+00,                1.224744871391589e+00,
     -1.650680123885785e+00, -0.524647623275290e+00,  0.524647623275290e+00,  1.650680123885785e+00,
     -2.020182870456086e+00, -0.958572464613819e+00,  0.0e+00,                0.958572464613819e+00,  2.020182870456086e+00,
     -2.350604973674492e+00, -1.335849074013697e+00, -0.436077411927617e+00,  0.436077411927617e+00,  1.335849074013697e+00,  2.350604973674492e+00,
     -2.651961356835233e+00, -1.673551628767471e+00, -0.816287882858965e+00,  0.0e+00,                0.816287882858965e+00,  1.673551628767471e+00,  2.651961356835233e+00,
     -2.930637420257244e+00, -1.981656756695843e+00, -1.157193712446780e+00, -0.381186990207322e+00,  0.381186990207322e+00,  1.157193712446780e+00,  1.981656756695843e+00,  2.930637420257244e+00,
     -3.190993201781528e+00, -2.266580584531843e+00, -1.468553289216668e+00, -0.723551018752838e+00,  0.000000000000000e+00,  0.723551018752838e+00,  1.468553289216668e+00,  2.266580584531843e+00,  3.190993201781528e+00,
     -3.436159118837738e+00, -2.532731674232790e+00, -1.756683649299882e+00, -1.036610829789514e+00, -3.429013272237046e-01,  3.429013272237046e-01,  1.036610829789514e+00,  1.756683649299882e+00,  2.532731674232790e+00,  3.436159118837738e+00 } ;

/* . The weights. */
const Real GHWEIGHTS[GHNDATA] =
    { 1.77245385090552e+00,
      0.8862269254528e+00,  0.8862269254528e+00,
      0.2954089751509e+00,  1.181635900604e+00,   0.2954089751509e+00,
      8.131283544725e-02,   8.049140900055e-01,   8.049140900055e-01,   8.131283544725e-02,
      1.995324205905e-02,   3.936193231522e-01,   9.453087204829e-01,   3.936193231522e-01,   1.995324205905e-02,
      4.530009905509e-03,   1.570673203229e-01,   7.246295952244e-01,   7.246295952244e-01,   1.570673203229e-01,   4.530009905509e-03,
      9.717812450995e-04,   5.451558281913e-02,   4.256072526101e-01,   8.102646175568e-01,   4.256072526101e-01,   5.451558281913e-02,   9.717812450995e-04,
      1.996040722114e-04,   1.707798300741e-02,   2.078023258149e-01,   6.611470125582e-01,   6.611470125582e-01,   2.078023258149e-01,   1.707798300741e-02,   1.996040722114e-04,
      3.960697726326e-05,   4.943624275537e-03,   8.847452739438e-02,   4.326515590026e-01,   7.202352156061e-01,   4.326515590026e-01,   8.847452739438e-02,   4.943624275537e-03,   3.960697726326e-05,
      7.640432855233e-06,   1.343645746781e-03,   3.387439445548e-02,   2.401386110823e-01,   6.108626337353e-01,   6.108626337353e-01,   2.401386110823e-01,   3.387439445548e-02,   1.343645746781e-03,   7.640432855233e-06 } ;

/*==================================================================================================================================
! . Utility functions.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Odd factorial.
! . (-1)!!, like 0!!, is defined to be one.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real OddFactorial ( const Integer n )
{
    Real f = 0.0e+00 ;
    if ( IsOdd ( n ) && ( n > -2 ) )
    {
        f = 1.0e+00 ;
        if ( n > 1 )
        {
            auto Integer i ;
            for ( i = 3 ; i <= n ; i += 2 ) { f *= ( Real ) i ; }
        }
    }
    return f ;
}

/*==================================================================================================================================
! . Standard procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
GaussianBasis *GaussianBasis_Allocate ( const Integer nShells )
{
    GaussianBasis *self = NULL ;
    if ( nShells > 0 )
    {
        auto Integer i ;
        self = Memory_AllocateType ( GaussianBasis ) ;
        self->pNormalized  =    True ;
        self->isSpherical  =    True ;
        self->basisType    = GaussianBasisType_Orbital ; /* . Default. */
        self->atomicNumber =      -1 ;
        self->lHigh        =       0 ;
        self->nBasis       =       0 ;
        self->nCBF         =       0 ;
        self->nShells      = nShells ;
        self->cbfPowX      = NULL    ;
        self->cbfPowY      = NULL    ;
        self->cbfPowZ      = NULL    ;
        self->shells       = Memory_AllocateArrayOfTypes ( nShells, Shell ) ;
        for ( i = 0 ; i < nShells ; i++ )
        {
           self->shells[i].lHigh       = 0 ;
           self->shells[i].lLow        = 0 ;
           self->shells[i].nBasis      = 0 ;
           self->shells[i].nCBF        = 0 ;
           self->shells[i].nPrimitives = 0 ;
           self->shells[i].nStart      = 0 ;
           self->shells[i].nStartC     = 0 ;
           self->shells[i].cbfPowX     = NULL ;
           self->shells[i].cbfPowY     = NULL ;
           self->shells[i].cbfPowZ     = NULL ;
           self->shells[i].c2s         = NULL ;
           self->shells[i].s2c         = NULL ;
           self->shells[i].primitives  = NULL ;
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocate a shell.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . No checking! */
void GaussianBasis_AllocateShell (       GaussianBasis *self        ,
                                   const Integer        iShell      ,
                                   const Integer        lHigh       ,
                                   const Integer        lLow        ,
                                   const Integer        nPrimitives )
{
    if ( self != NULL )
    {
        if ( ( iShell >= 0 ) && ( iShell < self->nShells ) )
        {
            auto Integer i, l ;
            /* . Set the shell options. */
            self->shells[iShell].lHigh       = lHigh       ;
            self->shells[iShell].lLow        = lLow        ;
            self->shells[iShell].nPrimitives = nPrimitives ;
            self->shells[iShell].nCBF        = SumOfCartesians ( lHigh ) - SumOfCartesians ( lLow - 1 ) ;
            if ( self->isSpherical ) self->shells[iShell].nBasis = SumOfSphericals ( lHigh ) - SumOfSphericals ( lLow - 1 ) ;
            else                     self->shells[iShell].nBasis = self->shells[iShell].nCBF ;
            /* . Transformations. */
            if ( self->isSpherical )
            {
                self->shells[iShell].c2s = GaussianBasis_TransformationCartesianToSpherical ( self->shells[iShell].lLow, self->shells[iShell].lHigh                          , NULL ) ;
                self->shells[iShell].s2c = GaussianBasis_TransformationSphericalToCartesian ( self->shells[iShell].lLow, self->shells[iShell].lHigh, self->shells[iShell].c2s, NULL ) ;
/*
# ifdef _PrintC2S_
printf ( "\nShell %d Spherical Transformations:\n", iShell ) ; fflush ( stdout ) ;
printf ( "Element %d, AM %d - %d, primitives %d, CBF %d, BF %d\n", self->atomicNumber, lLow, lHigh, nPrimitives, self->shells[iShell].nCBF, self->shells[iShell].nBasis ) ; fflush ( stdout ) ;
printf ( "C->S:" ) ; RealArray2D_Print ( self->shells[iShell].c2s ) ; fflush ( stdout ) ;
printf ( "S->C:" ) ; RealArray2D_Print ( self->shells[iShell].s2c ) ; fflush ( stdout ) ;
# endif
*/
            }
            /* . Allocate space for the primitives. */
            self->shells[iShell].primitives = Memory_AllocateArrayOfTypes ( self->shells[iShell].nPrimitives, Primitive ) ;
            /* . Allocate space for the coefficients. */
            l = self->shells[iShell].lHigh - self->shells[iShell].lLow + 1 ;
            for ( i = 0 ; i < self->shells[iShell].nPrimitives ; i++ )
            {
                self->shells[iShell].primitives[i].coefficients0 = Memory_AllocateArrayOfTypes ( l                        , Real ) ;
                self->shells[iShell].primitives[i].coefficients  = Memory_AllocateArrayOfTypes ( l                        , Real ) ;
                self->shells[iShell].primitives[i].cCBF          = Memory_AllocateArrayOfTypes ( self->shells[iShell].nCBF, Real ) ;
                Memory_Set ( self->shells[iShell].primitives[i].coefficients0, l                        , 0.0e+00 ) ;
                Memory_Set ( self->shells[iShell].primitives[i].coefficients , l                        , 0.0e+00 ) ;
                Memory_Set ( self->shells[iShell].primitives[i].cCBF         , self->shells[iShell].nCBF, 0.0e+00 ) ;
            }
            /* . Increment the shell and basis set function indices. */
            for ( i = self->nBasis = self->nCBF = 0 ; i < self->nShells ; i++ )
            {
                self->shells[i].nStartC = self->nCBF   ;
                self->shells[i].nStart  = self->nBasis ;
                self->nBasis += self->shells[i].nBasis ;
                self->nCBF   += self->shells[i].nCBF   ;
            }
            /* . Update the maximum angular momentum for the basis. */
            if ( self->shells[iShell].lHigh > self->lHigh ) self->lHigh = self->shells[iShell].lHigh ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
GaussianBasis *GaussianBasis_Clone ( const GaussianBasis *self, Status *status )
{
    GaussianBasis *new = NULL ;
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer c, i, nC, nCC, p ;
        /* . Basis data. */
        nC  = SumOfCartesians ( self->lHigh ) ;
        new = GaussianBasis_Allocate ( self->nShells ) ;
        new->pNormalized  = self->pNormalized  ;
        new->isSpherical  = self->isSpherical  ;
        new->basisType    = self->basisType    ;
        new->atomicNumber = self->atomicNumber ;
        new->lHigh        = self->lHigh        ;
        new->nBasis       = self->nBasis       ;
        new->nCBF         = self->nCBF         ;
        new->cbfPowX      = Integer_Clone ( self->cbfPowX, nC, status ) ;
        new->cbfPowY      = Integer_Clone ( self->cbfPowY, nC, status ) ;
        new->cbfPowZ      = Integer_Clone ( self->cbfPowZ, nC, status ) ;
        /* . Shell data. */
        for ( i = 0 ; i < new->nShells ; i++ )
        {
            new->shells[i].lHigh       = self->shells[i].lHigh       ;
            new->shells[i].lLow        = self->shells[i].lLow        ;
            new->shells[i].nBasis      = self->shells[i].nBasis      ;
            new->shells[i].nCBF        = self->shells[i].nCBF        ;
            new->shells[i].nPrimitives = self->shells[i].nPrimitives ;
            new->shells[i].nStart      = self->shells[i].nStart      ;
            new->shells[i].nStartC     = self->shells[i].nStartC     ;
            /* . Cartesian powers. */
            nC                         = SumOfCartesians ( new->shells[i].lLow - 1 ) ;
            new->shells[i].cbfPowX     = &new->cbfPowX[nC] ;
            new->shells[i].cbfPowY     = &new->cbfPowY[nC] ;
            new->shells[i].cbfPowZ     = &new->cbfPowZ[nC] ;
            /* . Transformations. */
            new->shells[i].c2s         = RealArray2D_CloneDeep ( self->shells[i].c2s, status ) ;
            new->shells[i].s2c         = RealArray2D_CloneDeep ( self->shells[i].s2c, status ) ;
            /* . Counters. */
            nC  = new->shells[i].lHigh - new->shells[i].lLow + 1 ;
            nCC = new->shells[i].nCBF ;
            /* . Primitive data. */
            new->shells[i].primitives = Memory_AllocateArrayOfTypes ( new->shells[i].nPrimitives, Primitive ) ;
            for ( p = 0 ; p < new->shells[i].nPrimitives ; p++ )
            {
                new->shells[i].primitives[p].exponent  = self->shells[i].primitives[p].exponent  ;
                new->shells[i].primitives[p].exponent0 = self->shells[i].primitives[p].exponent0 ;
                /* . Coefficients. */
                new->shells[i].primitives[p].coefficients0 = Memory_AllocateArrayOfTypes ( nC , Real ) ;
                new->shells[i].primitives[p].coefficients  = Memory_AllocateArrayOfTypes ( nC , Real ) ;
                new->shells[i].primitives[p].cCBF          = Memory_AllocateArrayOfTypes ( nCC, Real ) ;
                for ( c = 0 ; c < nC ; c++ )
                {
                    new->shells[i].primitives[p].coefficients0[c] = self->shells[i].primitives[p].coefficients0[c] ;
                    new->shells[i].primitives[p].coefficients [c] = self->shells[i].primitives[p].coefficients [c] ;
                }
                for ( c = 0 ; c < nCC ; c++ ) new->shells[i].primitives[p].cCBF[c] = self->shells[i].primitives[p].cCBF[c] ;
            }
        }
    }
    return new ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Deallocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasis_Deallocate ( GaussianBasis **self )
{
    if ( (*self) != NULL )
    {
        auto Integer p, s ;
        if ( (*self)->shells != NULL )
        {
            for ( s = 0 ; s < (*self)->nShells ; s++ )
            {
                if ( (*self)->shells[s].primitives != NULL )
                {
                    for ( p = 0 ; p < (*self)->shells[s].nPrimitives ; p++ )
                    {
                        Memory_Deallocate ( (*self)->shells[s].primitives[p].cCBF          ) ;
                        Memory_Deallocate ( (*self)->shells[s].primitives[p].coefficients  ) ;
                        Memory_Deallocate ( (*self)->shells[s].primitives[p].coefficients0 ) ;
                    }
                    Memory_Deallocate ( (*self)->shells[s].primitives ) ;
                }
                /* . Cartesian powers. */
                (*self)->shells[s].cbfPowX = NULL ;
                (*self)->shells[s].cbfPowY = NULL ;
                (*self)->shells[s].cbfPowZ = NULL ;
                /* . Transformations. */
                RealArray2D_Deallocate ( &((*self)->shells[s].c2s) ) ;
                RealArray2D_Deallocate ( &((*self)->shells[s].s2c) ) ;
            }
            Integer_Deallocate ( &((*self)->cbfPowX) ) ;
            Integer_Deallocate ( &((*self)->cbfPowY) ) ;
            Integer_Deallocate ( &((*self)->cbfPowZ) ) ;
            Memory_Deallocate  ( (*self)->shells  ) ;
        }
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Fill the primitive cCBF for the basis.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasis_FillPrimitiveCCBF ( GaussianBasis *self )
{
    if ( self != NULL )
    {
        auto Integer i, iShell, l, n, p ;
        for ( iShell = 0 ; iShell < self->nShells ; iShell++ )
        {
            for ( p = 0 ; p < self->shells[iShell].nPrimitives ; p++ )
            {
                for ( l = self->shells[iShell].lLow, n = 0 ; l <= self->shells[iShell].lHigh ; l++ )
                {
                    for ( i = 0 ; i < NumberOfCartesians ( l ) ; i++, n++ )
                    {
                        self->shells[iShell].primitives[p].cCBF[n] = self->shells[iShell].primitives[p].coefficients[l-self->shells[iShell].lLow] ;
                    }
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Finalize the basis ready for use.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasis_Finalize ( GaussianBasis *self, Status *status )
{
    GaussianBasis_FillPrimitiveCCBF      ( self         ) ;
    GaussianBasis_MakeCBFPowers          ( self, status ) ;
# ifdef _NormalizePrimitives_
    GaussianBasis_NormalizePrimitiveCCBF ( self, status ) ;
# endif
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Largest shell.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer GaussianBasis_LargestShell ( const GaussianBasis *self, const Boolean forC )
{
    Integer  n = 0 ;
    if ( self != NULL )
    {
        auto Integer i ;
        for ( i = 0 ; i < self->nShells ; i++ )
        {
            if ( forC ) n = Maximum ( n, self->shells[i].nCBF   ) ;
            else        n = Maximum ( n, self->shells[i].nBasis ) ;
        }
    }
    return n ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make the CBF powers for the basis.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasis_MakeCBFPowers ( GaussianBasis *self, Status *status )
{
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer iShell, l, n, x, y, z ;
        n = SumOfCartesians ( self->lHigh ) ;
        self->cbfPowX = Integer_Allocate ( n, status ) ;
        self->cbfPowY = Integer_Allocate ( n, status ) ;
        self->cbfPowZ = Integer_Allocate ( n, status ) ;
        if ( Status_IsOK ( status ) )
        {
# ifdef _PrintC2S_
printf ( "\nCBF powers (n, l, x, y, z):\n" ) ; fflush ( stdout ) ;
# endif
            for ( l = 0, n = 0 ; l <= self->lHigh ; l++ )
            {
                for ( z = 0 ; z <= l ; z++ )
                {
                    for ( y = 0 ; y <= (l-z) ; y++, n++ )
                    {
                        x = l - y - z ;
                        self->cbfPowX[n] = x ;
                        self->cbfPowY[n] = y ;
                        self->cbfPowZ[n] = z ;
# ifdef _PrintC2S_
printf ( "\n%d %d %d %d %d\n", n, l, x, y, z ) ; fflush ( stdout ) ;
# endif
                    }
                }
            }
            for ( iShell = 0 ; iShell < self->nShells ; iShell++ )
            {
                n = SumOfCartesians ( self->shells[iShell].lLow - 1 ) ;
                self->shells[iShell].cbfPowX = &self->cbfPowX[n] ;
                self->shells[iShell].cbfPowY = &self->cbfPowY[n] ;
                self->shells[iShell].cbfPowZ = &self->cbfPowZ[n] ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Normalize the Cartesian primitive coefficients.
! . The cbfPow arrays need to exist.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasis_NormalizePrimitiveCCBF ( GaussianBasis *self, Status *status )
{
    if ( self != NULL )
    {
        auto Integer  n ;
        auto Real    *factors ;
        /* . Determine the normalization factors. */
        n       = SumOfCartesians ( self->lHigh ) ;
        factors = Real_Allocate ( n, status ) ;
        if ( Status_IsOK ( status ) )
        {
            auto Integer i, i0, iShell, l, p, x, y, z ;
            for ( i = 0 ; i < n ; i++ )
            {
                x = self->cbfPowX[i] ;
                y = self->cbfPowY[i] ;
                z = self->cbfPowZ[i] ;
                l = x + y + z ;
                factors[i] = sqrt (   OddFactorial ( 2*l-1 ) /
                                    ( OddFactorial ( 2*x-1 ) *
                                      OddFactorial ( 2*y-1 ) *
                                      OddFactorial ( 2*z-1 ) ) ) ;
            }
            for ( iShell = 0 ; iShell < self->nShells ; iShell++ )
            {
                i0 = SumOfCartesians ( self->shells[iShell].lLow - 1 ) ;
                for ( p = 0 ; p < self->shells[iShell].nPrimitives ; p++ )
                {
                    for ( i = 0 ; i < self->shells[iShell].nCBF ; i++ )
                    {
                        self->shells[iShell].primitives[p].cCBF[i] *= factors[i+i0] ;
                    }
                }
            }
        }
        Real_Deallocate ( &factors ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Scale the exponents of a shell.
! . Scaling is done from exponent0.
! . The basis should be renormalized afterwards.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasis_ScaleShellExponents ( GaussianBasis *self, const Integer index, const Real zeta, Status *status )
{
    if ( self != NULL )
    {
        if ( index < self->nShells )
        {
            auto Integer  p ;
            auto Real     zeta2 = zeta * zeta ;
            for ( p = 0 ; p < self->shells[index].nPrimitives ; p++ )
            {
                self->shells[index].primitives[p].exponent = self->shells[index].primitives[p].exponent0 * zeta2 ;
            }
        }
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Generate the Cartesian to spherical harmonic transformation.
!---------------------------------------------------------------------------------------------------------------------------------*/
# ifdef _PrintC2S_
# include "IntegerArray2D.h"
# endif
RealArray2D *GaussianBasis_TransformationCartesianToSpherical ( const Integer lLow   ,
                                                                const Integer lHigh  ,
                                                                      Status *status )
{
    RealArray2D *result = NULL ;
    if ( ( lLow >= 0 ) && ( lHigh >= lLow ) && Status_IsOK ( status ) )
    {
        auto Integer      i, j, k, kmx2, ktmx, l, m, nc, ns, x, y, z ;
        auto Real         ab, c, d, *fac, g0, gm, gp ;
        auto RealArray1D *factorial = NULL ;
# ifdef _PrintC2S_
        auto IntegerArray2D *cXYZ ;
# endif
# ifdef _UnnormalizeTransformations_
        auto Real f ;
# endif
        /* . Calculate the factorial. */
        factorial = RealArray1D_AllocateWithExtent ( 2 * lHigh + 1, status ) ; fac = factorial->data ;
        RealArray1D_Set ( factorial, 1.0e+00 ) ;
        for ( i = 1 ; i < 2 * lHigh + 1 ; i++ ) fac[i] = ( Real ) i * fac[i-1] ;
        /* . Find the size of the transformation. */
        for ( l = lLow, nc = ns = 0 ; l <= lHigh ; l++ )
        {
            nc += ( ( l + 1 ) * ( l + 2 ) ) / 2 ;
            ns += 2 * l + 1 ;
        }
        /* . Allocate space. */
        result = RealArray2D_AllocateWithExtents ( nc, ns, NULL ) ;
        RealArray2D_Set ( result, 0.0e+00 ) ;
# ifdef _PrintC2S_
cXYZ = IntegerArray2D_AllocateWithExtents ( nc, 3, NULL );
# endif
        /* . Calculate the transformation. */
        for ( l = lLow, nc = ns = 0 ; l <= lHigh ; l++ )
        {
            for ( z = 0 ; z <= l ; z++ )
            {
                for ( y = 0 ; y <= l - z ; nc++, y++ )
                {
                    x = l - y - z ;
# ifdef _UnnormalizeTransformations_
                    f = sqrt (   OddFactorial ( 2*l-1 ) /
                               ( OddFactorial ( 2*x-1 ) *
                                 OddFactorial ( 2*y-1 ) *
                                 OddFactorial ( 2*z-1 ) ) ) ;
# endif
# ifdef _PrintC2S_
Array2D_Item ( cXYZ, nc, 0 ) = x ;
Array2D_Item ( cXYZ, nc, 1 ) = y ;
Array2D_Item ( cXYZ, nc, 2 ) = z ;
# endif
                    for ( m = 0 ; m <= l ; m++ )
                    {
                        j = ( l - m - z ) / 2 ;
                        if ( ( 2 * j ) == ( l - m - z ) )
                        {
                            g0 = gm = gp = 0.0e+00 ;
                            for ( i = 0 ; i <= ( l - m ) / 2 ; i++ )
                            {
                                c = 0.0e+00 ;
                                if ( ( j >= 0 ) && ( j <= i ) )
                                {
                                    c = ( fac[l]/(fac[i]*fac[l-i]) ) * ( fac[2*l-2*i]* pow ( -1.0e+00, i ) / fac[l-m-2*i] )* ( fac[i]/(fac[j]*fac[i-j]) ) ;
                                    for ( k = 0 ; k <= j ; k++ )
                                    {
                                        if ( ( x-2*k ) >= 0 && ( x-2*k ) <= m )
                                        {
                                            d = c * (fac[j]/(fac[k]*fac[j-k]))* (fac[m]/(fac[x-2*k]*fac[m+2*k-x])) ;
                                            if ( m == 0 )
                                            {
                                                kmx2 = k - x/2 ;
                                                if ( x % 2 == 0 ) g0 += d * pow ( -1.0e+00, kmx2 ) ;
                                            }
                                            else
                                            {
                                                ktmx = ( 2*k+m-x ) / 2 ;
                                                if ( abs ( m - x ) % 2 == 0 ) gp += d * sqrt ( 2.0e+00 ) * pow ( -1.0e+00, ktmx ) ;
                                                else                          gm += d * sqrt ( 2.0e+00 ) * pow ( -1.0e+00, ktmx ) ;
                                            }
                                        }
                                    }
                                }
                            }
                            ab = sqrt ( fac[2*x]*fac[2*y]*fac[2*z]*fac[l]/ (fac[2*l]*fac[x]*fac[y]*fac[z]) ) * sqrt ( fac[l-m]/fac[l+m] ) / ( pow ( 2, l ) * fac[l] ) ;
# ifdef _UnnormalizeTransformations_
                            ab *= f ;
# endif
                            if ( m == 0 )
                            {
                                Array2D_Item ( result, nc,       ns ) = ab * g0 ;
                            }
                            else
                            {
                                Array2D_Item ( result, nc, 2*m-1+ns ) = ab * gp ;
                                Array2D_Item ( result, nc, 2*m  +ns ) = ab * gm ;
                            }
                        }
                    }
                }
            }
            ns += ( 2 * l + 1 ) ;
        }
        /* . Deallocate space. */
        RealArray1D_Deallocate ( &factorial ) ;
# ifdef _PrintC2S_
printf ( "\nC->S conversion for l_Low=%d to l_High=%d:\n", lLow, lHigh ) ; fflush ( stdout ) ;
IntegerArray2D_Print ( cXYZ ) ; fflush ( stdout ) ;
RealArray2D_Print ( result ) ; fflush ( stdout ) ;
IntegerArray2D_Deallocate ( &cXYZ ) ;
# endif
    }
    return result ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Generate the spherical harmonic to Cartesian transformation.
!---------------------------------------------------------------------------------------------------------------------------------*/
RealArray2D *GaussianBasis_TransformationSphericalToCartesian ( const Integer      lLow   ,
                                                                const Integer      lHigh  ,
                                                                const RealArray2D *c2s    ,
                                                                      Status      *status )
{
    RealArray2D *result = NULL ;
    if ( ( lLow >= 0 ) && ( lHigh >= lLow ) && ( c2s != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer      i, j, k, l, nc, ns, x, x1, x2, y, y1, y2, z, z1, z2 ;
        auto Real         a1, a2, *fac, s ;
        auto RealArray1D *factorial = NULL ;
# ifdef _UnnormalizeTransformations_
        auto Real f1, f2 ;
# endif
        /* . Calculate the factorial. */
        factorial = RealArray1D_AllocateWithExtent ( 2 * lHigh + 1, status ) ; fac = factorial->data ;
        RealArray1D_Set ( factorial, 1.0e+00 ) ;
        for ( i = 1 ; i < 2 * lHigh + 1 ; i++ ) fac[i] = ( Real ) i * fac[i-1] ;
        /* . Allocate space. */
        result = RealArray2D_AllocateWithExtents ( View2D_Rows ( c2s ), View2D_Columns ( c2s ), NULL ) ;
        RealArray2D_Set ( result, 0.0e+00 ) ;
        /* . Find the size of the transformation. */
        for ( l = lLow, nc = ns = 0 ; l <= lHigh ; l++ )
        {
            for ( i = z1 = 0 ; z1 <= l ; z1++ )
            {
                for ( y1 = 0 ; y1 <= l - z1 ; i++, y1++ )
                {
                    x1 = l - y1 - z1 ;
                    a1 = sqrt ( ( fac[x1] * fac[y1] * fac[z1] ) / ( fac[2*x1] * fac[2*y1] * fac[2*z1] ) ) ;
# ifdef _UnnormalizeTransformations_
                    f1 = sqrt (   OddFactorial ( 2*l -1 ) /
                                ( OddFactorial ( 2*x1-1 ) *
                                  OddFactorial ( 2*y1-1 ) *
                                  OddFactorial ( 2*z1-1 ) ) ) ;
# endif
                    for ( j = z2 = 0 ; z2 <= l ; z2++ )
                    {
                        for ( y2 = 0 ; y2 <= l - z2 ; j++, y2++ )
                        {
                            x2 = l - y2 - z2 ;
                            a2 = sqrt ( fac[x2] * fac[y2] * fac[z2] / ( fac[2*x2] * fac[2*y2] * fac[2*z2] ) ) ;
# ifdef _UnnormalizeTransformations_
                            f2 = sqrt (   OddFactorial ( 2*l -1 ) /
                                        ( OddFactorial ( 2*x2-1 ) *
                                          OddFactorial ( 2*y2-1 ) *
                                          OddFactorial ( 2*z2-1 ) ) ) ;
# endif
                            x  = x1 + x2 ;
                            y  = y1 + y2 ;
                            z  = z1 + z2 ;
                            if ( ( x % 2 == 0 ) && ( y % 2 == 0 ) && ( z % 2 == 0 ) )
                            {
                                s = a1 * a2 * fac[x] * fac[y] * fac[z] / ( fac[x/2] * fac[y/2] * fac[z/2] ) ;
# ifdef _UnnormalizeTransformations_
                                s /= ( f1 * f2 ) ;
# endif
                                for ( k = 0 ; k < 2 * l + 1 ; k++ ) Array2D_Item ( result, i + nc, k + ns ) += s * Array2D_Item ( c2s, j + nc, k + ns ) ;
                            }
                        }
                    }
                }
            }
            nc += ( ( l + 1 ) * ( l + 2 ) ) / 2 ;
            ns += ( 2 * l + 1 ) ;
        }
        /* . Deallocate space. */
        RealArray1D_Deallocate ( &factorial ) ;
# ifdef _PrintC2S_
printf ( "\nS->C conversion for l_Low=%d to l_High=%d:\n", lLow, lHigh ) ; fflush ( stdout ) ;
RealArray2D_Print ( result ) ; fflush ( stdout ) ;
# endif
    }
    return result ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Unnormalize the primitives of a basis if pNormalized is True.
! . Exponents0 is used to change coefficients0 into coefficients.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasis_UnnormalizePrimitives ( GaussianBasis *self )
{
    if ( ( self != NULL ) && ( self->pNormalized ) )
    {
        auto Integer i, iShell, l, p ;
        auto Real    e2, fac ;
        /* . Loop over the coefficients of each angular momentum for each shell. */
        /* . Shells have same primitive exponents but different primitive coefficients for each angular momentum. */
        for ( iShell = 0 ; iShell < self->nShells ; iShell++ )
        {
            for ( l = self->shells[iShell].lLow ; l <= self->shells[iShell].lHigh ; l++ )
            {
                for ( p = 0 ; p < self->shells[iShell].nPrimitives ; p++ )
                {
                    e2  = 2.0e+00 * self->shells[iShell].primitives[p].exponent0 ;
                    fac = PI32 / ( e2 * sqrt ( e2 ) ) ;
                    for ( i = 1 ; i <= l ; i++ ) fac *= ( Real ) ( 2*i - 1 ) / ( 2.0e+00 * e2 ) ;
                    self->shells[iShell].primitives[p].coefficients [l-self->shells[iShell].lLow] =
                  ( self->shells[iShell].primitives[p].coefficients0[l-self->shells[iShell].lLow] / sqrt ( fac ) ) ;
                }
            }
        }
    }
}
