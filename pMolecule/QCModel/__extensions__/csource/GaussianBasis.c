/*==================================================================================================================================
! . This module handles Gaussian basis functions
!===================================================================================================================================
!
! . Notes:
!
!   1. Basis sets may be Cartesian or spherical harmonical.
!
!   2. Shells can be of a single angular momentum type or have multiple
!      consecutive values (from am_low to am_high). E.g. s, p, d, f, g, sp, spd
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
!      m = 0, 1, -1, 2, -2, 3, -3, etc.
!
!   5. The normalization factors are calculated using:
!
!      Sqrt ( F0 F0 Fl / Fi Fj Fk ) = Sqrt ( ( 2l - 1 )!! / (2i-1)!! (2j-1)!! (2k-1)!! )
!
!      where Fi = Sqrt ( Pi / 2 a ) * (2i - 1)!! / ( 4^i * a^2i ) and a is the basis function exponent.
!
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "GaussianBasis.h"
# include "Memory.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cartesian basis function and shell parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The first and last Cartesian basis functions of a particular angular momentum. Equivalent to MAXCBFSUM. */
const Integer  CBFSTART[MAXAMP1] = { 0, 1, 4, 10, 20 } ;
const Integer  CBFSTOP[MAXAMP1]  = { 0, 3, 9, 19, 34 } ;

/* . The x,y and z powers of the Cartesian basis functions. */
const Integer  CBFPOWX[MAXCBFSUM] = { 0,  1, 0, 0,  2, 1, 0, 1, 0, 0,  3, 2, 1, 0, 2, 1, 0, 1, 0, 0,  4, 3, 2, 1, 0, 3, 2, 1, 0, 2, 1, 0, 1, 0, 0 } ;
const Integer  CBFPOWY[MAXCBFSUM] = { 0,  0, 1, 0,  0, 1, 2, 0, 1, 0,  0, 1, 2, 3, 0, 1, 2, 0, 1, 0,  0, 1, 2, 3, 4, 0, 1, 2, 3, 0, 1, 2, 0, 1, 0 } ;
const Integer  CBFPOWZ[MAXCBFSUM] = { 0,  0, 0, 1,  0, 0, 0, 1, 1, 2,  0, 0, 0, 0, 1, 1, 1, 2, 2, 3,  0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4 } ;

/* . The shell definitions. */
# define NSHELLTYPES 8
const ShellDefinition SHELLTYPES_CBF[NSHELLTYPES] = { { 0, 0,  0,  1,  1 } ,
                                                      { 1, 1,  1,  3,  3 } ,
                                                      { 2, 2,  4,  6,  6 } ,
                                                      { 3, 3, 10, 10, 10 } ,
                                                      { 4, 4, 20, 15, 15 } ,
                                                      { 0, 1,  0,  4,  4 } ,
                                                      { 0, 2,  0, 10, 10 } ,
                                                      { 0, 3,  0, 20, 20 } } ;
const ShellDefinition SHELLTYPES_SPH[NSHELLTYPES] = { { 0, 0,  0,  1,  1 } ,
                                                      { 1, 1,  1,  3,  3 } ,
                                                      { 2, 2,  4,  5,  6 } ,
                                                      { 3, 3, 10,  7, 10 } ,
                                                      { 4, 4, 20,  9, 15 } ,
                                                      { 0, 1,  0,  4,  4 } ,
                                                      { 0, 2,  0,  9, 10 } ,
                                                      { 0, 3,  0, 16, 20 } } ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Gauss-Hermite quadrature parameters (used for KE and overlap integrals).
!
! . For two center overlaps, the maximum number of points is (i+j+2)/2-1
! . where i is ammax+1 and j is ammax+2 (ammax = maximum angular momentum).
! . This implies 2C-overlaps (and derivatives) can be calculated up to ammax
! . = 10 (m functions).
!
! . For three center overlaps, the maximum is (i+j+k+2)/2 where i is ammax+1,
! . j is ammax+1 and k is ammax. (These are for overlaps only - not KE
! . integrals). This implies can go up to h functions.
!
! . If necessary extra weights are easy to calculate or to find.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The indices into the data. */
const Integer  GHFIRST[GHMAXPT] = { 0, 1, 3, 6, 10, 15, 21, 28, 36, 45 } ;
const Integer  GHLAST[GHMAXPT]  = { 0, 2, 5, 9, 14, 20, 27, 35, 44, 54 } ;

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

/*----------------------------------------------------------------------------------------------------------------------------------
! . Static procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static RealArray2D *CartesianToSphericalTransformation ( const Integer  amlow, const Integer  amhigh ) ;
static RealArray2D *SphericalToCartesianTransformation ( const Integer  amlow, const Integer  amhigh, const RealArray2D *c2s ) ;

/*==================================================================================================================================
! . Standard procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
GaussianBasis *GaussianBasis_Allocate ( const Integer  nshells )
{
    GaussianBasis *self = NULL ;
    if ( nshells > 0 )
    {
        auto Integer  i ;
        self = Memory_AllocateType ( GaussianBasis ) ;
        self->QNORMALIZEDPRIMITIVES   =    True ;
        self->QSPHERICAL              =    True ;
        self->QTOSPHERICAL            =   False ; /* . For the moment everything is done in Cartesians. */
        self->atomicNumber            =      -1 ;
        self->maximum_angularmomentum =       0 ;
        self->nbasis                  =       0 ;
        self->nbasisw                 =       0 ;
        self->nshells                 = nshells ;
        self->basisType               = GaussianBasisType_Orbital   ; /* . Defaults. */
        self->normalizationType       = NormalizationType_Symmetric ;
        self->c2o                     =    NULL ;
        self->o2c                     =    NULL ;
        self->shells                  = Memory_AllocateArrayOfTypes ( nshells, Shell ) ;
        for ( i = 0 ; i < nshells ; i++ )
        {
           self->shells[i].nbasisw     = 0 ;
           self->shells[i].nprimitives = 0 ;
           self->shells[i].nstart      = 0 ;
           self->shells[i].nstartw     = 0 ;
           self->shells[i].c2s         = NULL ;
           self->shells[i].s2c         = NULL ;
           self->shells[i].primitives  = NULL ;
           self->shells[i].type        = NULL ;
        }
    }
    return self ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocate a shell.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* # define DEBUGPRINT */
void GaussianBasis_AllocateShell ( GaussianBasis *self, const Integer  ishell, const Integer  nprimitives, const Integer  typeindex )
{
    if ( self != NULL )
    {
        if ( ( ishell >= 0 ) && ( ishell < self->nshells ) )
        {
            auto Integer   am, i ;
            /* . Set the shell options. */
            self->shells[ishell].nprimitives = nprimitives ;
            /* . Get the shell type. */
            if ( self->QSPHERICAL ) self->shells[ishell].type = &(SHELLTYPES_SPH[typeindex]) ;
            else                    self->shells[ishell].type = &(SHELLTYPES_CBF[typeindex]) ;
            /* . Transformations. */
            if ( self->QSPHERICAL )
            {
                self->shells[ishell].c2s = CartesianToSphericalTransformation ( self->shells[ishell].type->angularmomentum_low, self->shells[ishell].type->angularmomentum_high                           ) ;
                self->shells[ishell].s2c = SphericalToCartesianTransformation ( self->shells[ishell].type->angularmomentum_low, self->shells[ishell].type->angularmomentum_high, self->shells[ishell].c2s ) ;
            }
            /* . Allocate space for the primitives. */
            self->shells[ishell].primitives = Memory_AllocateArrayOfTypes ( self->shells[ishell].nprimitives, Primitive ) ;
            /* . Allocate space for the coefficients. */
            am = self->shells[ishell].type->angularmomentum_high - self->shells[ishell].type->angularmomentum_low + 1 ;
            for ( i = 0 ; i < self->shells[ishell].nprimitives ; i++ )
            {
                self->shells[ishell].primitives[i].coefficients0 = Memory_AllocateArrayOfTypes ( am                             , Real ) ;
                self->shells[ishell].primitives[i].coefficients  = Memory_AllocateArrayOfTypes ( am                             , Real ) ;
                self->shells[ishell].primitives[i].ccbf          = Memory_AllocateArrayOfTypes ( self->shells[ishell].type->ncbf, Real ) ;
                Memory_Set ( self->shells[ishell].primitives[i].coefficients0, am                             , 0.0e+00 ) ;
                Memory_Set ( self->shells[ishell].primitives[i].coefficients , am                             , 0.0e+00 ) ;
                Memory_Set ( self->shells[ishell].primitives[i].ccbf         , self->shells[ishell].type->ncbf, 0.0e+00 ) ;
            }
            /* . Increment the shell and basis set function indices. */
            for ( i = self->nbasis = self->nbasisw = 0 ; i < self->nshells ; i++ )
            {
                self->shells[i].nstartw = self->nbasisw ;
                self->shells[i].nstart  = self->nbasis  ;
                if ( self->shells[i].type != NULL )
                {
                    self->nbasis += self->shells[i].type->nbasis ;
                    if ( self->QTOSPHERICAL )
                    {
                        self->shells[i].nbasisw  = self->shells[i].type->nbasis ;
                        self->nbasisw           += self->shells[i].type->nbasis ;
                    }
                    else
                    {
                        self->shells[i].nbasisw  = self->shells[i].type->ncbf ;
                        self->nbasisw           += self->shells[i].type->ncbf ;
                    }
                }
            }
            /* . Update the maximum angular momentum for the basis. */
            if ( self->shells[ishell].type->angularmomentum_high > self->maximum_angularmomentum ) self->maximum_angularmomentum = self->shells[ishell].type->angularmomentum_high ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cloning.
!---------------------------------------------------------------------------------------------------------------------------------*/
GaussianBasis *GaussianBasis_Clone ( const GaussianBasis *self )
{
    GaussianBasis *new = NULL ;
    if ( self != NULL )
    {
        auto Integer  c, i, nc, ncc, p ;
        /* . Basis data. */
        new = GaussianBasis_Allocate ( self->nshells ) ;
        new->QNORMALIZEDPRIMITIVES   = self->QNORMALIZEDPRIMITIVES   ;
        new->QSPHERICAL              = self->QSPHERICAL              ;
        new->QTOSPHERICAL            = self->QTOSPHERICAL            ;
        new->atomicNumber            = self->atomicNumber            ;
        new->maximum_angularmomentum = self->maximum_angularmomentum ;
        new->nbasis                  = self->nbasis                  ;
        new->nbasisw                 = self->nbasisw                 ;
        new->basisType               = self->basisType               ;
        new->normalizationType       = self->normalizationType       ;
        /* . Shell data. */
        for ( i = 0 ; i < new->nshells ; i++ )
        {
            new->shells[i].nbasisw     = self->shells[i].nbasisw     ;
            new->shells[i].nprimitives = self->shells[i].nprimitives ;
            new->shells[i].nstart      = self->shells[i].nstart      ;
            new->shells[i].nstartw     = self->shells[i].nstartw     ;
            new->shells[i].type        = self->shells[i].type        ;
            /* . Transformations. */
            new->shells[i].c2s = RealArray2D_CloneDeep ( self->shells[i].c2s, NULL ) ;
            new->shells[i].s2c = RealArray2D_CloneDeep ( self->shells[i].s2c, NULL ) ;
            /* . Counters. */
            nc  = new->shells[i].type->angularmomentum_high - new->shells[i].type->angularmomentum_low + 1 ;
            ncc = new->shells[i].type->ncbf ;
            /* . Primitive data. */
            new->shells[i].primitives = Memory_AllocateArrayOfTypes ( new->shells[i].nprimitives, Primitive ) ;
            for ( p = 0 ; p < new->shells[i].nprimitives ; p++ )
            {
                new->shells[i].primitives[p].exponent  = self->shells[i].primitives[p].exponent  ;
                new->shells[i].primitives[p].exponent0 = self->shells[i].primitives[p].exponent0 ;
                /* . Coefficients. */
                new->shells[i].primitives[p].coefficients0 = Memory_AllocateArrayOfTypes ( nc , Real ) ;
                new->shells[i].primitives[p].coefficients  = Memory_AllocateArrayOfTypes ( nc , Real ) ;
                new->shells[i].primitives[p].ccbf          = Memory_AllocateArrayOfTypes ( ncc, Real ) ;
                for ( c = 0 ; c < nc ; c++ )
                {
                    new->shells[i].primitives[p].coefficients0[c] = self->shells[i].primitives[p].coefficients0[c] ;
                    new->shells[i].primitives[p].coefficients [c] = self->shells[i].primitives[p].coefficients [c] ;
                }
                for ( c = 0 ; c < ncc ; c++ ) new->shells[i].primitives[p].ccbf[c] = self->shells[i].primitives[p].ccbf[c] ;
            }
        }
        /* . Transformations. */
        new->c2o = RealArray2D_CloneDeep ( self->c2o, NULL ) ;
        new->o2c = RealArray2D_CloneDeep ( self->o2c, NULL ) ;
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
        auto Integer  p, s ;
        if ( (*self)->shells != NULL )
        {
            for ( s = 0 ; s < (*self)->nshells ; s++ )
            {
                if ( (*self)->shells[s].primitives != NULL )
                {
                    for ( p = 0 ; p < (*self)->shells[s].nprimitives ; p++ )
                    {
                        Memory_Deallocate ( (*self)->shells[s].primitives[p].ccbf          ) ;
                        Memory_Deallocate ( (*self)->shells[s].primitives[p].coefficients  ) ;
                        Memory_Deallocate ( (*self)->shells[s].primitives[p].coefficients0 ) ;
                    }
                    Memory_Deallocate ( (*self)->shells[s].primitives ) ;
                }
                (*self)->shells[s].type = NULL ;
                /* . Transformations. */
                RealArray2D_Deallocate ( &((*self)->shells[s].c2s) ) ;
                RealArray2D_Deallocate ( &((*self)->shells[s].s2c) ) ;
            }
            Memory_Deallocate ( (*self)->shells ) ;
            /* . Transformations. */
            RealArray2D_Deallocate ( &((*self)->c2o) ) ;
            RealArray2D_Deallocate ( &((*self)->o2c) ) ;
        }
        Memory_Deallocate ( (*self) ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Fill the primitive ccbf for the basis.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasis_FillPrimitiveCCBF ( GaussianBasis *self )
{
    if ( self != NULL )
    {
        auto Integer  am, i, ishell, n, p ;
        for ( ishell = 0 ; ishell < self->nshells ; ishell++ )
        {
            for ( p = 0 ; p < self->shells[ishell].nprimitives ; p++ )
            {
                for ( am = self->shells[ishell].type->angularmomentum_low, n = 0 ; am <= self->shells[ishell].type->angularmomentum_high ; am++ )
                {
                    for ( i = 0 ; i < NUMBEROFCARTESIANS ( am ) ; i++, n++ )
                    {
                        self->shells[ishell].primitives[p].ccbf[n] = self->shells[ishell].primitives[p].coefficients[am-self->shells[ishell].type->angularmomentum_low] ;
                    }
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Scale the exponents of a shell.
! . Scaling is done from exponent0.
! . The basis should be renormalized afterwards.
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasis_ScaleShellExponents ( GaussianBasis *self, const Integer  index, const Real zeta, Status *status )
{
    if ( self != NULL )
    {
        if ( index < self->nshells )
        {
            auto Integer  p ;
            auto Real     zeta2 = zeta * zeta ;
            for ( p = 0 ; p < self->shells[index].nprimitives ; p++ )
            {
                self->shells[index].primitives[p].exponent = self->shells[index].primitives[p].exponent0 * zeta2 ;
            }
        }
        else Status_Set ( status, Status_IndexOutOfRange ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Unnormalize the primitives of a basis if QNORMALIZEDPRIMITIVES is True.
! . Exponents0 is used to change coefficients0 into coefficients.
! . The general expression for a Gaussian * x^l y^m z^n is:
!
! . Sqrt ( 2^(2*(l+m+n)+3/2) * zeta^(l+m+n+3/2) / (2l-1)!! (2m-1)!! (2n-1)!! pi^(3/2) )
!
!---------------------------------------------------------------------------------------------------------------------------------*/
void GaussianBasis_UnnormalizePrimitives ( GaussianBasis *self )
{
    if ( ( self != NULL ) && ( self->QNORMALIZEDPRIMITIVES ) )
    {
        auto Integer  am, i, ishell, p ;
        auto Real     ex, fac ;
        /* . Loop over the coefficients of each angular momentum for each shell. */
        for ( ishell = 0 ; ishell < self->nshells ; ishell++ )
        {
            for ( am = self->shells[ishell].type->angularmomentum_low ; am <= self->shells[ishell].type->angularmomentum_high ; am++ )
            {
                for ( p = 0 ; p < self->shells[ishell].nprimitives ; p++ )
                {
                    ex  = 2.0e+00 * self->shells[ishell].primitives[p].exponent0 ;
                    fac = PI32 / ( ex * sqrt ( ex ) ) ;
                    for ( i = 1 ; i <= am ; i++ ) fac *= ( Real ) ( 2*i - 1 ) / ( 2.0e+00 * ex ) ;
                    self->shells[ishell].primitives[p].coefficients [am-self->shells[ishell].type->angularmomentum_low] =
                    self->shells[ishell].primitives[p].coefficients0[am-self->shells[ishell].type->angularmomentum_low] / sqrt ( fac ) ;
                }
            }
        }
    }
}

/*==================================================================================================================================
! . Cartesian to spherical transformation procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Generate the Cartesian to spherical transformation.
!---------------------------------------------------------------------------------------------------------------------------------*/
static RealArray2D *CartesianToSphericalTransformation ( const Integer  amlow, const Integer  amhigh )
{
    RealArray2D *result = NULL ;
    if ( ( amlow >= 0 ) && ( amhigh >= amlow ) )
    {
        auto Integer      i, j, k, kmx2, ktmx, l, m, nc, ns, x, y, z ;
        auto Real         ab, c, d, *fac, g0, gm, gp ;
        auto RealArray1D *factorial = NULL ;

        /* . Calculate the factorial. */
        factorial = RealArray1D_AllocateWithExtent ( 2 * amhigh + 1, NULL ) ; fac = factorial->data ;
        RealArray1D_Set ( factorial, 1.0e+00 ) ;
        for ( i = 1 ; i < 2 * amhigh + 1 ; i++ ) fac[i] = ( Real ) i * fac[i-1] ;

        /* . Find the size of the transformation. */
        for ( l = amlow, nc = ns = 0 ; l <= amhigh ; l++ )
        {
            nc += ( ( l + 1 ) * ( l + 2 ) ) / 2 ;
            ns += 2 * l + 1 ;
        }

        /* . Allocate space. */
        result = RealArray2D_AllocateWithExtents ( nc, ns, NULL ) ;
        RealArray2D_Set ( result, 0.0e+00 ) ;

        /* . Calculate the transformation. */
        for ( l = amlow, nc = ns = 0 ; l <= amhigh ; l++ )
        {
            for ( z = 0 ; z <= l ; z++ )
            {
                for ( y = 0 ; y <= l - z ; nc++, y++ )
                {
                    x  = l - y - z ;
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
    }
    return result ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Generate the spherical to Cartesian transformation.
!---------------------------------------------------------------------------------------------------------------------------------*/
static RealArray2D *SphericalToCartesianTransformation ( const Integer  amlow, const Integer  amhigh, const RealArray2D *c2s )
{
    RealArray2D *result = NULL ;
    if ( ( amlow >= 0 ) && ( amhigh >= amlow ) && ( c2s != NULL ) )
    {
        auto Integer      i, j, k, l, nc, ns, x, x1, x2, y, y1, y2, z, z1, z2 ;
        auto Real         a1, a2, *fac, s ;
        auto RealArray1D *factorial = NULL ;

        /* . Calculate the factorial. */
        factorial = RealArray1D_AllocateWithExtent ( 2 * amhigh + 1, NULL ) ; fac = factorial->data ;
        RealArray1D_Set ( factorial, 1.0e+00 ) ;
        for ( i = 1 ; i < 2 * amhigh + 1 ; i++ ) fac[i] = ( Real ) i * fac[i-1] ;

        /* . Allocate space. */
        result = RealArray2D_AllocateWithExtents ( View2D_Rows ( c2s ), View2D_Columns ( c2s ), NULL ) ;
        RealArray2D_Set ( result, 0.0e+00 ) ;

        /* . Find the size of the transformation. */
        for ( l = amlow, nc = ns = 0 ; l <= amhigh ; l++ )
        {
            for ( i = z1 = 0 ; z1 <= l ; z1++ )
            {
                for ( y1 = 0 ; y1 <= l - z1 ; i++, y1++ )
                {
                    x1 = l - y1 - z1 ;
                    a1 = sqrt ( ( fac[x1] * fac[y1] * fac[z1] ) / ( fac[2*x1] * fac[2*y1] * fac[2*z1] ) ) ;
                    for ( j = z2 = 0 ; z2 <= l ; z2++ )
                    {
                        for ( y2 = 0 ; y2 <= l - z2 ; j++, y2++ )
                        {
                            x2 = l - y2 - z2 ;
                            a2 = sqrt ( fac[x2] * fac[y2] * fac[z2] / ( fac[2*x2] * fac[2*y2] * fac[2*z2] ) ) ;
                            x  = x1 + x2 ;
                            y  = y1 + y2 ;
                            z  = z1 + z2 ;
                            if ( ( x % 2 == 0 ) && ( y % 2 == 0 ) && ( z % 2 == 0 ) )
                            {
                                s = a1 * a2 * fac[x] * fac[y] * fac[z] / ( fac[x/2] * fac[y/2] * fac[z/2] ) ;
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
    }
    return result ;
}

