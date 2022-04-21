/*==================================================================================================================================
! . The Gaussian basis module.
!=================================================================================================================================*/
# ifndef _GAUSSIANBASIS
# define _GAUSSIANBASIS

# include <stdio.h>

# include "Boolean.h"
# include "Integer.h"
# include "IntegerUtilities.h"
# include "Real.h"
# include "RealArray2D.h"
# include "RealUtilities.h"
# include "Status.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Constants.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define PI    3.14159265358979e+00
# define PI12  1.77245385090551e+00
# define PI252 3.49868366552497e+01 /* . Equivalent to 2 * power ( pi, 5./2. ). */
# define PI32  5.56832799683170e+00
# define RLN10 2.30258509299405e+00

/*----------------------------------------------------------------------------------------------------------------------------------
! . Enumerations.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Operators. */
typedef enum {
    GaussianBasisOperator_AntiCoulomb = 1 ,
    GaussianBasisOperator_Coulomb     = 2 ,
    GaussianBasisOperator_Dipole      = 3 ,
    GaussianBasisOperator_Kinetic     = 4 ,
    GaussianBasisOperator_Overlap     = 5 ,
    GaussianBasisOperator_Poisson     = 6 ,
    GaussianBasisOperator_Quadrupole  = 7
} GaussianBasisOperator ;

/* . Basis types. */
typedef enum {
    GaussianBasisType_Density = 1 , /* . Coulomb orthogonalization by default. */
    GaussianBasisType_Orbital = 2   /* . Overlap orthogonalization. */
} GaussianBasisType ;

 /* . Poisson density fit bases have Poisson orthogonalization but their <ij|f> integrals are overlaps. */

/*----------------------------------------------------------------------------------------------------------------------------------
! . Miscellaneous Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The number of Gauss-Hermite quadrature points - these may need to be increased for I functions! */
# define GHMAXPT 10
# define GHNDATA ( GHMAXPT * ( GHMAXPT + 1 ) ) / 2

/* . Integral tolerances. */
# define INTEGRAL_TOLERANCE            1.0e-12
# define INVERSE_FIT_TOLERANCE         1.0e-5
# define PRIMITIVE_OVERLAP_TOLERANCE ( RLN10 * 18 )

/* . Angular momenta. */
# define MAXIMUM_ANGULAR_MOMENTUM 6 /* . Up to I functions - should be sufficient for all (?) bases in the basis set exchange. */
# define MAXAMP1 ( MAXIMUM_ANGULAR_MOMENTUM + 1 )
# define MAXAMP2 ( MAXIMUM_ANGULAR_MOMENTUM + 2 )
# define MAXAMP3 ( MAXIMUM_ANGULAR_MOMENTUM + 3 )
# define MAXAMP4 ( MAXIMUM_ANGULAR_MOMENTUM + 4 )
# define MAXAMP5 ( MAXIMUM_ANGULAR_MOMENTUM + 5 )
# define MAXAMP6 ( MAXIMUM_ANGULAR_MOMENTUM + 6 )
# define MAXAMP7 ( MAXIMUM_ANGULAR_MOMENTUM + 7 )

/* . The number of Cartesian and Spherical Harmonic functions for a given angular momentum. */
# define NumberOfCartesians( l ) ( ( ( (l)+1 ) * ( (l)+2 ) ) / 2 )
# define NumberOfSphericals( l ) ( 2 * (l) + 1 )

/* . The sum of Cartesian and Spherical Harmonic functions up to and including a given angular momentum. */
/* . These functions work correctly for l = -1. */
# define SumOfCartesians( l ) ( ( ( (l)+1 ) * ( (l)+2 ) * ( (l)+3 ) ) / 6 )
# define SumOfSphericals( l ) ( ( (l) + 1 ) * ( (l) + 1 ) )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Basis Data.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The primitive type. */
typedef struct
{
    Real *cCBF          ;
    Real *coefficients  ;
    Real *coefficients0 ; /* . Input coefficients and exponents are always unchanged. */
    Real  exponent      ;
    Real  exponent0     ;
} Primitive ;

/* . The shell type. */
typedef struct
{
    Integer      lHigh       ;
    Integer      lLow        ;
    Integer      nBasis      ;
    Integer      nCBF        ;
    Integer      nPrimitives ;
    Integer      nStart      ;
    Integer      nStartC     ;
    Integer     *cbfPowX     ;
    Integer     *cbfPowY     ;
    Integer     *cbfPowZ     ;
    RealArray2D *c2s         ; /* . Cartesian and spherical harmonic transformations. */
    RealArray2D *s2c         ;
    Primitive   *primitives  ;
} Shell ;

/* . The basis type. */
typedef struct
{
    Boolean            isSpherical  ; /* . Cartesian or spherical basis. */
    Boolean            pNormalized  ; /* . This flag refers to input coefficients only. Internal coefficients always correspond to unnormalized primitives. */
    GaussianBasisType  basisType    ;
    Integer            atomicNumber ;
    Integer            lHigh        ;
    Integer            nBasis       ;
    Integer            nCBF         ;
    Integer            nShells      ;
    Integer           *cbfPowX      ;
    Integer           *cbfPowY      ;
    Integer           *cbfPowZ      ;
    Shell             *shells       ;
} GaussianBasis ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Gauss-Hermite quadrature parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The indices into the data. */
extern const Integer GHINDEX [GHMAXPT+1] ;

/* . The abscissae and weights. */
extern const Real GHABSCISSAE[GHNDATA] ;
extern const Real GHWEIGHTS  [GHNDATA] ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
extern GaussianBasis *GaussianBasis_Allocate                           ( const Integer         nShells     ) ;
extern void           GaussianBasis_AllocateShell                      (       GaussianBasis  *self        ,
                                                                         const Integer         iShell      ,
                                                                         const Integer         lHigh       ,
                                                                         const Integer         lLow        ,
                                                                         const Integer         nPrimitives ) ;
extern GaussianBasis *GaussianBasis_Clone                              ( const GaussianBasis  *self        ,
                                                                               Status         *status      ) ;
extern void           GaussianBasis_Deallocate                         (       GaussianBasis **self        ) ;
extern void           GaussianBasis_FillPrimitiveCCBF                  (       GaussianBasis  *self        ) ;
extern void           GaussianBasis_Finalize                           (       GaussianBasis  *self        ,
                                                                               Status         *status      ) ;
extern Integer        GaussianBasis_LargestShell                       ( const GaussianBasis  *self        ,
                                                                         const Boolean         forC        ) ;
extern void           GaussianBasis_MakeCBFPowers                      (       GaussianBasis  *self        ,
                                                                               Status         *status      ) ;
extern void           GaussianBasis_NormalizePrimitiveCCBF             (       GaussianBasis  *self        ,
                                                                               Status         *status      ) ;
extern void           GaussianBasis_ScaleShellExponents                (       GaussianBasis  *self        ,
                                                                         const Integer         index       ,
                                                                         const Real            zeta        ,
                                                                               Status         *status      ) ;
extern RealArray2D   *GaussianBasis_TransformationCartesianToSpherical ( const Integer         lLow        ,
                                                                         const Integer         lHigh       ,
                                                                               Status         *status      ) ;
extern RealArray2D   *GaussianBasis_TransformationSphericalToCartesian ( const Integer         lLow        ,
                                                                         const Integer         lHigh       ,
                                                                         const RealArray2D    *c2s         ,
                                                                               Status         *status      ) ;
extern void           GaussianBasis_UnnormalizePrimitives              (       GaussianBasis  *self        ) ;

# endif
