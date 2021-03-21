/*==================================================================================================================================
! . The Gaussian basis module.
!=================================================================================================================================*/
# ifndef _GAUSSIANBASIS
# define _GAUSSIANBASIS

# include <stdio.h>

# include "Boolean.h"
# include "Integer.h"
# include "Real.h"
# include "RealArray2D.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Constants.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define PI    3.14159265358979e+00

# define PI252 3.49868366552497e+01 /* . Equivalent to 2 * power ( pi, 5./2. ). */
# define PI32  5.56832799683170e+00
# define RLN10 2.30258509299405e+00

/*----------------------------------------------------------------------------------------------------------------------------------
! . Enumerations.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Basis types. */
typedef enum {
    GaussianBasisType_Coulomb = 1 ,
    GaussianBasisType_Orbital = 2 ,
    GaussianBasisType_Poisson = 3
} GaussianBasisType ;

/* . Normalization types. */
typedef enum {
    NormalizationType_Canonical = 1 ,
    NormalizationType_Diagonal  = 2 ,
    NormalizationType_Symmetric = 3
} NormalizationType ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Miscellaneous Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Integral tolerances. */
# define INTEGRAL_TOLERANCE            1.0e-12
# define INVERSE_FIT_TOLERANCE         1.0e-5
# define PRIMITIVE_OVERLAP_TOLERANCE ( RLN10 * 18 )

/* . Currently up to g functions. */
# define MAXIMUM_ANGULAR_MOMENTUM 4
# define MAXAMP1                  ( MAXIMUM_ANGULAR_MOMENTUM + 1 )
# define MAXAMP2                  ( MAXIMUM_ANGULAR_MOMENTUM + 2 )
# define MAXAMP3                  ( MAXIMUM_ANGULAR_MOMENTUM + 3 )
# define MAXAMP4                  ( MAXIMUM_ANGULAR_MOMENTUM + 4 )

/* . The number of Cartesian basis functions for a shell of the maximum angular momentum. */
# define MAXCBF ( MAXAMP1 * ( MAXAMP1 + 1 ) ) / 2

/* . The sum of the number of Cartesian basis functions up to the maximum angular momentum. */
# define MAXCBFSUM ( MAXAMP1 * ( MAXAMP1 + 1 ) * ( MAXAMP1 + 2 ) ) / 6

/* . The maximum number of Gauss-Hermite quadrature points. */
# define GHMAXPT 10
# define GHNDATA ( GHMAXPT * ( GHMAXPT + 1 ) ) / 2

/* . The number of Cartesian and Spherical Harmonic functions for a given angular momentum. */
# define NUMBEROFCARTESIANS( l ) ( ( ( (l)+1 ) * ( (l)+2 ) ) / 2 )
# define NUMBEROFSPHERICALS( l ) ( 2 * (l) + 1 )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Basis Data.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The primitive type. */
typedef struct
{
    Real *ccbf          ;
    Real *coefficients  ;
    Real *coefficients0 ; /* . Input coefficients and exponents are always unchanged. */
    Real  exponent      ;
    Real  exponent0     ;
} Primitive ;

/* . The shell definition type. */
typedef struct {
    Integer angularmomentum_low  ;
    Integer angularmomentum_high ;
    Integer cbfindex             ;
    Integer nbasis               ;
    Integer ncbf                 ;
/* . Include pointers to transformation matrices here. */
} ShellDefinition ;

/* . The shell type. */
typedef struct
{
    Integer                nbasisw     ;
    Integer                nprimitives ;
    Integer                nstart      ;
    Integer                nstartw     ;
    RealArray2D           *c2s         ; /* . Cartesian and spherical harmonic transformations. */
    RealArray2D           *s2c         ;
    Primitive             *primitives  ;
    const ShellDefinition *type        ;
} Shell ;

/* . The basis type. */
typedef struct
{
    Boolean            QNORMALIZEDPRIMITIVES   ; /* . This flag refers to input coefficients only. Internal coefficients always correspond to unnormalized primitives. */
    Boolean            QSPHERICAL              ; /* . Cartesian or spherical basis. */
    Boolean            QTOSPHERICAL            ; /* . Flag indicating whether integrals, etc. should be transformed to a spherical harmonic basis (if the basis is spherical) or left in the default Cartesian representation. */
    Integer            atomicNumber            ;
    Integer            maximum_angularmomentum ;
    Integer            nbasis                  ;
    Integer            nbasisw                 ;
    Integer            nshells                 ;
    GaussianBasisType  basisType               ;
    NormalizationType  normalizationType       ;
    RealArray2D       *c2o                     ;  /* . Cartesian to orthogonal transformations. */
    RealArray2D       *o2c                     ;
    Shell             *shells                  ;
} GaussianBasis ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Cartesian basis function and shell data.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The x,y and z powers of the Cartesian basis functions. */
extern const Integer CBFPOWX[MAXCBFSUM] ;
extern const Integer CBFPOWY[MAXCBFSUM] ;
extern const Integer CBFPOWZ[MAXCBFSUM] ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Gauss-Hermite quadrature parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The indices into the data. */
extern const Integer GHFIRST[GHMAXPT] ;
extern const Integer GHLAST[GHMAXPT]  ;

/* . The abscissae. */
extern const Real GHABSCISSAE[GHNDATA] ;

/* . The weights. */
extern const Real GHWEIGHTS[GHNDATA] ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern GaussianBasis *GaussianBasis_Allocate              ( const Integer         nshells     ) ; 
extern void           GaussianBasis_AllocateShell         (       GaussianBasis  *self        ,   
                                                            const Integer         ishell      ,   
                                                            const Integer         nprimitives ,   
                                                            const Integer         typeindex   ) ; 
extern GaussianBasis *GaussianBasis_Clone                 ( const GaussianBasis  *self        ) ; 
extern void           GaussianBasis_Deallocate            (       GaussianBasis **self        ) ; 
extern void           GaussianBasis_FillPrimitiveCCBF     (       GaussianBasis  *self        ) ;
extern void           GaussianBasis_ScaleShellExponents   (       GaussianBasis  *self        ,
                                                            const Integer         index       ,
                                                            const Real            zeta        ,
                                                                  Status         *status      ) ;
extern void           GaussianBasis_UnnormalizePrimitives (       GaussianBasis  *self        ) ;

# endif
