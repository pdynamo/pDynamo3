/*==================================================================================================================================
! . Utility procedures for calculating the integrals in a MNDO method.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>

# include "MNDODefinitions.h"
# include "MNDOIntegralDefinitions.h"
# include "MNDOIntegralUtilities.h"
# include "MNDOParameters.h"
# include "NumericalMacros.h"
# include "Units.h"

# define MNDODORBITALS

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the transformation matrices for a given atom pair (i-j).
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . Phi zero when yji = 0.
! . ca  = cos(phi)     , sa  = sin(phi)
! . cb  = cos(theta)   , sb  = sin(theta)
! . c2a = cos(2*phi)   , s2a = sin(2*phi)
! . c2b = cos(2*theta) , s2b = sin(2*phi)
!
! . There is a problem when atoms are aligned on z-axis as phi is undefined.
! . Therefore, axes are swapped.
!
! . Note that a NULL transformation implies the identity whereas a NULL derivative transformation implies a zero matrix!
*/
# define _AlignmentTolerance 0.99999999e+00
void MNDOIntegralUtilities_GetTransformationMatrices ( const Integer       ni               ,
                                                       const Integer       nj               ,
                                                       const Real          r                ,
                                                       const Real          x                ,
                                                       const Real          y                ,
                                                       const Real          z                ,
                                                             RealArray2D **iTransformation  ,
                                                             RealArray2D **jTransformation  ,
                                                             RealArray2D **iTransformationX ,
                                                             RealArray2D **iTransformationY ,
                                                             RealArray2D **iTransformationZ ,
                                                             RealArray2D **jTransformationX ,
                                                             RealArray2D **jTransformationY ,
                                                             RealArray2D **jTransformationZ )
{
    Boolean      doGradients ;
    Integer      norbitals   ;
    RealArray2D *it = NULL, *itX = NULL, *itY = NULL, *itZ = NULL, *jt = NULL, *jtX = NULL, *jtY = NULL, *jtZ = NULL ;
    /* . Initialization. */
    doGradients = ( iTransformationX != NULL ) ||
                  ( iTransformationY != NULL ) ||
                  ( iTransformationZ != NULL ) ||
                  ( jTransformationX != NULL ) ||
                  ( jTransformationY != NULL ) ||
                  ( jTransformationZ != NULL ) ;
    norbitals = Maximum ( ni, nj ) ;
    /* . Only continue if there are p orbitals or higher and the distance is big enough. */
    if ( ( norbitals > 1 ) && ( r > SMALL_RIJ ) )
    {
        auto Boolean      axesSwapped ;
        auto Integer      i, ij, j, m, mn, n ;
        auto Real         b, b2, ca, cb, c2a, c2b, sa, sb, s2a, s2b, xji, yji, zji ;
        auto Real         caX = 0.0e+00, caY = 0.0e+00, cbX = 0.0e+00, cbY = 0.0e+00, cbZ = 0.0e+00,
                          saX = 0.0e+00, saY = 0.0e+00, sbX = 0.0e+00, sbY = 0.0e+00, sbZ = 0.0e+00 ;
        auto RealArray2D *clesser         = NULL, *clesserX         = NULL, *clesserY         = NULL, *clesserZ         = NULL,
                         *cTransformation = NULL, *cTransformationX = NULL, *cTransformationY = NULL, *cTransformationZ = NULL,
                         *oTransformation = NULL, *oTransformationX = NULL, *oTransformationY = NULL, *oTransformationZ = NULL,
                         *sTransformation = NULL ;
        /* . Allocate the orbital transformation matrix and initialize to the ss case. */
        oTransformation = RealArray2D_AllocateWithExtents ( norbitals, norbitals, NULL ) ;
        RealArray2D_Set ( oTransformation, 0.0e+00 ) ;
        Array2D_Item ( oTransformation, S, S ) = 1.0e+00 ;
        if ( doGradients )
        {
            oTransformationX = RealArray2D_AllocateWithExtents ( norbitals, norbitals, NULL ) ; RealArray2D_Set ( oTransformationX, 0.0e+00 ) ;
            oTransformationY = RealArray2D_AllocateWithExtents ( norbitals, norbitals, NULL ) ; RealArray2D_Set ( oTransformationY, 0.0e+00 ) ;
            oTransformationZ = RealArray2D_AllocateWithExtents ( norbitals, norbitals, NULL ) ; RealArray2D_Set ( oTransformationZ, 0.0e+00 ) ;
        }
        /* . Check for z-axis alignment. */
        xji = x ; yji = y ; zji = z ;
        axesSwapped = ( fabs ( zji / r ) > _AlignmentTolerance ) ;
        if ( axesSwapped )
        {
            b = xji ; xji = zji ; zji = b ;
            sTransformation = RealArray2D_AllocateWithExtents ( norbitals, norbitals, NULL ) ;
            RealArray2D_Set ( sTransformation, 0.0e+00 ) ;
        }
        /* . p-orbital transformation. */
        b2 = xji * xji + yji * yji ;
        b  = sqrt ( b2 ) ;
        sb = b / r ;
        ca = xji / b ;
        sa = yji / b ;
        cb = zji / r   ;
        /* . Rotation matrix. */
        Array2D_Item ( oTransformation,  PX, PSIGMA   ) = ca*sb ;
        Array2D_Item ( oTransformation,  PX, PPIPLUS  ) = ca*cb ;
        Array2D_Item ( oTransformation,  PX, PPIMINUS ) = -sa   ;
        Array2D_Item ( oTransformation,  PY, PSIGMA   ) = sa*sb ;
        Array2D_Item ( oTransformation,  PY, PPIPLUS  ) = sa*cb ;
        Array2D_Item ( oTransformation,  PY, PPIMINUS ) = ca    ;
        Array2D_Item ( oTransformation,  PZ, PSIGMA   ) = cb    ;
        Array2D_Item ( oTransformation,  PZ, PPIPLUS  ) = -sb   ;
/*            Array2D_Item ( oTransformation,  PZ, PPIMINUS ) = 0.0 ; */
        /* . Derivatives. */
        if ( doGradients )
        {
            auto Real r2 = r * r ;
            /* . Various factors. */
            caX =   yji * yji / ( b2 * b ) ;
            caY = - xji * yji / ( b2 * b ) ;
            saX =   caY ;
            saY =   xji * xji / ( b2 * b ) ;
            cbX = - xji * zji       / (     r2 * r ) ;
            cbY = - yji * zji       / (     r2 * r ) ;
            cbZ =   b2              / (     r2 * r ) ;
            sbX =   xji * zji * zji / ( b * r2 * r ) ;
            sbY =   yji * zji * zji / ( b * r2 * r ) ;
            sbZ = - b * zji         / (     r2 * r ) ;
            /* . X. */
            Array2D_Item ( oTransformationX,  PX, PSIGMA   ) = caX*sb + ca*sbX ;
            Array2D_Item ( oTransformationX,  PX, PPIPLUS  ) = caX*cb + ca*cbX ;
            Array2D_Item ( oTransformationX,  PX, PPIMINUS ) = -saX ;
            Array2D_Item ( oTransformationX,  PY, PSIGMA   ) = saX*sb + sa*sbX ;
            Array2D_Item ( oTransformationX,  PY, PPIPLUS  ) = saX*cb + sa*cbX ;
            Array2D_Item ( oTransformationX,  PY, PPIMINUS ) = caX  ;
            Array2D_Item ( oTransformationX,  PZ, PSIGMA   ) = cbX  ;
            Array2D_Item ( oTransformationX,  PZ, PPIPLUS  ) = -sbX ;
            /* . Y. */
            Array2D_Item ( oTransformationY,  PX, PSIGMA   ) = caY*sb + ca*sbY ;
            Array2D_Item ( oTransformationY,  PX, PPIPLUS  ) = caY*cb + ca*cbY ;
            Array2D_Item ( oTransformationY,  PX, PPIMINUS ) = -saY ;
            Array2D_Item ( oTransformationY,  PY, PSIGMA   ) = saY*sb + sa*sbY ;
            Array2D_Item ( oTransformationY,  PY, PPIPLUS  ) = saY*cb + sa*cbY ;
            Array2D_Item ( oTransformationY,  PY, PPIMINUS ) = caY  ;
            Array2D_Item ( oTransformationY,  PZ, PSIGMA   ) = cbY  ;
            Array2D_Item ( oTransformationY,  PZ, PPIPLUS  ) = -sbY ;
            /* . Z. */
            Array2D_Item ( oTransformationZ,  PX, PSIGMA   ) = ca*sbZ ;
            Array2D_Item ( oTransformationZ,  PX, PPIPLUS  ) = ca*cbZ ;
            Array2D_Item ( oTransformationZ,  PY, PSIGMA   ) = sa*sbZ ;
            Array2D_Item ( oTransformationZ,  PY, PPIPLUS  ) = sa*cbZ ;
            Array2D_Item ( oTransformationZ,  PZ, PSIGMA   ) = cbZ  ;
            Array2D_Item ( oTransformationZ,  PZ, PPIPLUS  ) = -sbZ ;
        }
        /* . Swap transformation. */
        if ( axesSwapped )
        {
            Array2D_Item ( sTransformation,  S,  S ) = 1.0e+00 ;
            Array2D_Item ( sTransformation, PX, PZ ) = 1.0e+00 ;
            Array2D_Item ( sTransformation, PY, PY ) = 1.0e+00 ;
            Array2D_Item ( sTransformation, PZ, PX ) = 1.0e+00 ;
        }
        /* . d-orbital transformation. */
        if ( norbitals == 9 )
        {
            auto Real pt5sq3 ;
            /* . Initialization. */
            pt5sq3 = 0.5e+00 * sqrt ( 3.0e+00 ) ;
            c2a = 2.0e+00 * ca * ca - 1.0e+00 ;
            c2b = 2.0e+00 * cb * cb - 1.0e+00 ;
            s2a = 2.0e+00 * sa * ca ;
            s2b = 2.0e+00 * sb * cb ;
            /* . Rotation matrix. */
            Array2D_Item ( oTransformation, DX2Y2, DSIGMA      ) = pt5sq3 * c2a * sb * sb ;
            Array2D_Item ( oTransformation, DX2Y2, DPIPLUS     ) = 0.5e+00 * c2a * s2b ;
            Array2D_Item ( oTransformation, DX2Y2, DPIMINUS    ) = -s2a * sb ;
            Array2D_Item ( oTransformation, DX2Y2, DDELTAPLUS  ) = c2a * ( cb * cb + 0.5e+00 * sb * sb ) ;
            Array2D_Item ( oTransformation, DX2Y2, DDELTAMINUS ) = -s2a * cb ;
            Array2D_Item ( oTransformation, DXZ  , DSIGMA      ) = pt5sq3 * ca * s2b ;
            Array2D_Item ( oTransformation, DXZ  , DPIPLUS     ) = ca * c2b ;
            Array2D_Item ( oTransformation, DXZ  , DPIMINUS    ) = -sa * cb ;
            Array2D_Item ( oTransformation, DXZ  , DDELTAPLUS  ) = -0.5e+00 * ca * s2b ;
            Array2D_Item ( oTransformation, DXZ  , DDELTAMINUS ) = sa * sb ;
            Array2D_Item ( oTransformation, DZ2  , DSIGMA      ) = cb * cb - 0.5e+00 * sb * sb ;
            Array2D_Item ( oTransformation, DZ2  , DPIPLUS     ) = -pt5sq3 * s2b ;
/*                Array2D_Item ( oTransformation, DZ2  , DPIMINUS    ) = 0.0e+00 ; */
            Array2D_Item ( oTransformation, DZ2  , DDELTAPLUS  ) = pt5sq3 * sb * sb ;
/*                Array2D_Item ( oTransformation, DZ2  , DDELTAMINUS ) = 0.0e+00 ; */
            Array2D_Item ( oTransformation, DYZ  , DSIGMA      ) = pt5sq3 * sa * s2b ;
            Array2D_Item ( oTransformation, DYZ  , DPIPLUS     ) = sa * c2b ;
            Array2D_Item ( oTransformation, DYZ  , DPIMINUS    ) = ca * cb ;
            Array2D_Item ( oTransformation, DYZ  , DDELTAPLUS  ) = -0.5e+00 * sa * s2b ;
            Array2D_Item ( oTransformation, DYZ  , DDELTAMINUS ) = -ca * sb ;
            Array2D_Item ( oTransformation, DXY  , DSIGMA      ) = pt5sq3 * s2a * sb * sb ;
            Array2D_Item ( oTransformation, DXY  , DPIPLUS     ) = 0.5e+00 * s2a * s2b ;
            Array2D_Item ( oTransformation, DXY  , DPIMINUS    ) = c2a * sb ;
            Array2D_Item ( oTransformation, DXY  , DDELTAPLUS  ) = s2a * ( cb * cb + 0.5e+00 * sb * sb ) ;
            Array2D_Item ( oTransformation, DXY  , DDELTAMINUS ) = c2a * cb ;
            /* . Derivatives. */
            if ( doGradients )
            {
                auto Real c2aX, c2aY, c2bX, c2bY, c2bZ, s2aX, s2aY, s2bX, s2bY, s2bZ ;
                /* . Various factors. */
                c2aX = 4.0e+00 * ca * caX ;
                c2aY = 4.0e+00 * ca * caY ;
                s2aX = 2.0e+00 * ( saX * ca + sa * caX ) ;
                s2aY = 2.0e+00 * ( saY * ca + sa * caY ) ;
                c2bX = 4.0e+00 * cb * cbX ;
                c2bY = 4.0e+00 * cb * cbY ;
                c2bZ = 4.0e+00 * cb * cbZ ;
                s2bX = 2.0e+00 * ( sbX * cb + sb * cbX ) ;
                s2bY = 2.0e+00 * ( sbY * cb + sb * cbY ) ;
                s2bZ = 2.0e+00 * ( sbZ * cb + sb * cbZ ) ;
                /* . X. */
                Array2D_Item ( oTransformationX, DX2Y2, DSIGMA      ) = pt5sq3 * ( c2aX * sb * sb + 2.0e+00 * c2a * sb * sbX ) ;
                Array2D_Item ( oTransformationX, DX2Y2, DPIPLUS     ) = 0.5e+00 * ( c2aX * s2b + c2a * s2bX ) ;
                Array2D_Item ( oTransformationX, DX2Y2, DPIMINUS    ) = - s2aX * sb - s2a * sbX ;
                Array2D_Item ( oTransformationX, DX2Y2, DDELTAPLUS  ) = c2aX * ( cb * cb + 0.5e+00 * sb * sb ) + c2a * ( 2.0e+00 * cb * cbX + sb * sbX ) ;
                Array2D_Item ( oTransformationX, DX2Y2, DDELTAMINUS ) = - s2aX * cb - s2a * cbX ;
                Array2D_Item ( oTransformationX, DXZ  , DSIGMA      ) = pt5sq3 * ( caX * s2b + ca * s2bX ) ;
                Array2D_Item ( oTransformationX, DXZ  , DPIPLUS     ) = caX * c2b + ca * c2bX ;
                Array2D_Item ( oTransformationX, DXZ  , DPIMINUS    ) = - saX * cb - sa * cbX ;
                Array2D_Item ( oTransformationX, DXZ  , DDELTAPLUS  ) = - 0.5e+00 * ( caX * s2b + ca * s2bX ) ;
                Array2D_Item ( oTransformationX, DXZ  , DDELTAMINUS ) = saX * sb + sa * sbX ;
                Array2D_Item ( oTransformationX, DZ2  , DSIGMA      ) = 2.0e+00 * cb * cbX - sb * sbX ;
                Array2D_Item ( oTransformationX, DZ2  , DPIPLUS     ) = - pt5sq3 * s2bX ;
                Array2D_Item ( oTransformationX, DZ2  , DDELTAPLUS  ) = pt5sq3 * 2.0e+00 * sb * sbX ;
                Array2D_Item ( oTransformationX, DYZ  , DSIGMA      ) = pt5sq3 * ( saX * s2b + sa * s2bX ) ;
                Array2D_Item ( oTransformationX, DYZ  , DPIPLUS     ) = saX * c2b + sa * c2bX ;
                Array2D_Item ( oTransformationX, DYZ  , DPIMINUS    ) = caX * cb + ca * cbX ;
                Array2D_Item ( oTransformationX, DYZ  , DDELTAPLUS  ) = - 0.5e+00 * ( saX * s2b + sa * s2bX ) ;
                Array2D_Item ( oTransformationX, DYZ  , DDELTAMINUS ) = - caX * sb - ca * sbX ;
                Array2D_Item ( oTransformationX, DXY  , DSIGMA      ) = pt5sq3 * ( s2aX * sb * sb + 2.0e+00 * s2a * sb * sbX ) ;
                Array2D_Item ( oTransformationX, DXY  , DPIPLUS     ) = 0.5e+00 * ( s2aX * s2b + s2a * s2bX ) ;
                Array2D_Item ( oTransformationX, DXY  , DPIMINUS    ) = c2aX * sb + c2a * sbX ;
                Array2D_Item ( oTransformationX, DXY  , DDELTAPLUS  ) = s2aX * ( cb * cb + 0.5e+00 * sb * sb ) + s2a * ( 2.0e+00 * cb * cbX + sb * sbX ) ;
                Array2D_Item ( oTransformationX, DXY  , DDELTAMINUS ) = c2aX * cb + c2a * cbX ;
                /* . Y. */
                Array2D_Item ( oTransformationY, DX2Y2, DSIGMA      ) = pt5sq3 * ( c2aY * sb * sb + 2.0e+00 * c2a * sb * sbY ) ;
                Array2D_Item ( oTransformationY, DX2Y2, DPIPLUS     ) = 0.5e+00 * ( c2aY * s2b + c2a * s2bY ) ;
                Array2D_Item ( oTransformationY, DX2Y2, DPIMINUS    ) = - s2aY * sb - s2a * sbY ;
                Array2D_Item ( oTransformationY, DX2Y2, DDELTAPLUS  ) = c2aY * ( cb * cb + 0.5e+00 * sb * sb ) + c2a * ( 2.0e+00 * cb * cbY + sb * sbY ) ;
                Array2D_Item ( oTransformationY, DX2Y2, DDELTAMINUS ) = - s2aY * cb - s2a * cbY ;
                Array2D_Item ( oTransformationY, DXZ  , DSIGMA      ) = pt5sq3 * ( caY * s2b + ca * s2bY ) ;
                Array2D_Item ( oTransformationY, DXZ  , DPIPLUS     ) = caY * c2b + ca * c2bY ;
                Array2D_Item ( oTransformationY, DXZ  , DPIMINUS    ) = - saY * cb - sa * cbY ;
                Array2D_Item ( oTransformationY, DXZ  , DDELTAPLUS  ) = - 0.5e+00 * ( caY * s2b + ca * s2bY ) ;
                Array2D_Item ( oTransformationY, DXZ  , DDELTAMINUS ) = saY * sb + sa * sbY ;
                Array2D_Item ( oTransformationY, DZ2  , DSIGMA      ) = 2.0e+00 * cb * cbY - sb * sbY ;
                Array2D_Item ( oTransformationY, DZ2  , DPIPLUS     ) = - pt5sq3 * s2bY ;
                Array2D_Item ( oTransformationY, DZ2  , DDELTAPLUS  ) = pt5sq3 * 2.0e+00 * sb * sbY ;
                Array2D_Item ( oTransformationY, DYZ  , DSIGMA      ) = pt5sq3 * ( saY * s2b + sa * s2bY ) ;
                Array2D_Item ( oTransformationY, DYZ  , DPIPLUS     ) = saY * c2b + sa * c2bY ;
                Array2D_Item ( oTransformationY, DYZ  , DPIMINUS    ) = caY * cb + ca * cbY ;
                Array2D_Item ( oTransformationY, DYZ  , DDELTAPLUS  ) = - 0.5e+00 * ( saY * s2b + sa * s2bY ) ;
                Array2D_Item ( oTransformationY, DYZ  , DDELTAMINUS ) = - caY * sb - ca * sbY ;
                Array2D_Item ( oTransformationY, DXY  , DSIGMA      ) = pt5sq3 * ( s2aY * sb * sb + 2.0e+00 * s2a * sb * sbY ) ;
                Array2D_Item ( oTransformationY, DXY  , DPIPLUS     ) = 0.5e+00 * ( s2aY * s2b + s2a * s2bY ) ;
                Array2D_Item ( oTransformationY, DXY  , DPIMINUS    ) = c2aY * sb + c2a * sbY ;
                Array2D_Item ( oTransformationY, DXY  , DDELTAPLUS  ) = s2aY * ( cb * cb + 0.5e+00 * sb * sb ) + s2a * ( 2.0e+00 * cb * cbY + sb * sbY ) ;
                Array2D_Item ( oTransformationY, DXY  , DDELTAMINUS ) = c2aY * cb + c2a * cbY ;
                /* . Z. */
                Array2D_Item ( oTransformationZ, DX2Y2, DSIGMA      ) = pt5sq3 * 2.0e+00 * c2a * sb * sbZ ;
                Array2D_Item ( oTransformationZ, DX2Y2, DPIPLUS     ) = 0.5e+00 * c2a * s2bZ ;
                Array2D_Item ( oTransformationZ, DX2Y2, DPIMINUS    ) = - s2a * sbZ ;
                Array2D_Item ( oTransformationZ, DX2Y2, DDELTAPLUS  ) = c2a * ( 2.0e+00 * cb * cbZ + sb * sbZ ) ;
                Array2D_Item ( oTransformationZ, DX2Y2, DDELTAMINUS ) = - s2a * cbZ ;
                Array2D_Item ( oTransformationZ, DXZ  , DSIGMA      ) = pt5sq3 * ca * s2bZ ;
                Array2D_Item ( oTransformationZ, DXZ  , DPIPLUS     ) = ca * c2bZ ;
                Array2D_Item ( oTransformationZ, DXZ  , DPIMINUS    ) = - sa * cbZ ;
                Array2D_Item ( oTransformationZ, DXZ  , DDELTAPLUS  ) = - 0.5e+00 * ca * s2bZ ;
                Array2D_Item ( oTransformationZ, DXZ  , DDELTAMINUS ) = sa * sbZ ;
                Array2D_Item ( oTransformationZ, DZ2  , DSIGMA      ) = 2.0e+00 * cb * cbZ - sb * sbZ ;
                Array2D_Item ( oTransformationZ, DZ2  , DPIPLUS     ) = - pt5sq3 * s2bZ ;
                Array2D_Item ( oTransformationZ, DZ2  , DDELTAPLUS  ) = pt5sq3 * 2.0e+00 * sb * sbZ ;
                Array2D_Item ( oTransformationZ, DYZ  , DSIGMA      ) = pt5sq3 * sa * s2bZ ;
                Array2D_Item ( oTransformationZ, DYZ  , DPIPLUS     ) = sa * c2bZ ;
                Array2D_Item ( oTransformationZ, DYZ  , DPIMINUS    ) = ca * cbZ ;
                Array2D_Item ( oTransformationZ, DYZ  , DDELTAPLUS  ) = - 0.5e+00 * sa * s2bZ ;
                Array2D_Item ( oTransformationZ, DYZ  , DDELTAMINUS ) = - ca * sbZ ;
                Array2D_Item ( oTransformationZ, DXY  , DSIGMA      ) = pt5sq3 * 2.0e+00 * s2a * sb * sbZ ;
                Array2D_Item ( oTransformationZ, DXY  , DPIPLUS     ) = 0.5e+00 * s2a * s2bZ ;
                Array2D_Item ( oTransformationZ, DXY  , DPIMINUS    ) = c2a * sbZ ;
                Array2D_Item ( oTransformationZ, DXY  , DDELTAPLUS  ) = s2a * ( 2.0e+00 * cb * cbZ + sb * sbZ ) ;
                Array2D_Item ( oTransformationZ, DXY  , DDELTAMINUS ) = c2a * cbZ ;
            }
            /* . Swap transformation. */
            if ( axesSwapped )
            {
                Array2D_Item ( sTransformation, DXZ,   DXZ   ) =  1.0e+00 ;
                Array2D_Item ( sTransformation, DXY,   DYZ   ) =  1.0e+00 ;
                Array2D_Item ( sTransformation, DYZ,   DXY   ) =  1.0e+00 ;
                Array2D_Item ( sTransformation, DX2Y2, DX2Y2 ) =  0.5e+00 ;
                Array2D_Item ( sTransformation, DX2Y2, DZ2   ) =  pt5sq3  ;
                Array2D_Item ( sTransformation, DZ2,   DX2Y2 ) =  pt5sq3  ;
                Array2D_Item ( sTransformation, DZ2,   DZ2   ) = -0.5e+00 ;
            }
        }
        /* . Swap transformation. */
        if ( axesSwapped )
        {
            auto RealArray2D *new = NULL, *swap = NULL ;
            /* . Premultiply the existing transformation by the swap transformation. */
            new  = RealArray2D_AllocateWithExtents ( norbitals, norbitals, NULL ) ;
            RealArray2D_MatrixMultiply ( False, False, 1.0e+00, sTransformation, oTransformation, 0.0e+00, new, NULL ) ;
            swap = oTransformation ; oTransformation = new ; new = swap ;
            /* . Same for gradients but also swap X and Z matrices. */
            if ( doGradients )
            {
                RealArray2D_MatrixMultiply ( False, False, 1.0e+00, sTransformation, oTransformationX, 0.0e+00, new, NULL ) ; swap = oTransformationX ; oTransformationX = new ; new = swap ;
                RealArray2D_MatrixMultiply ( False, False, 1.0e+00, sTransformation, oTransformationY, 0.0e+00, new, NULL ) ; swap = oTransformationY ; oTransformationY = new ; new = swap ;
                RealArray2D_MatrixMultiply ( False, False, 1.0e+00, sTransformation, oTransformationZ, 0.0e+00, new, NULL ) ; swap = oTransformationZ ; oTransformationZ = new ; new = swap ;
                swap = oTransformationZ ; oTransformationZ = oTransformationX ; oTransformationX = swap ;
            }
            /* . Clear up. */
            RealArray2D_Deallocate ( &new ) ;
        }
        /* . Allocate the largest transformation matrix. */
        n = ( norbitals * ( norbitals + 1 ) ) / 2 ;
        cTransformation = RealArray2D_AllocateWithExtents ( n, n, NULL ) ;
        RealArray2D_Set ( cTransformation, 0.0e+00 ) ;
        /* . Form the matrix. */
        for ( i = ij = 0 ; i < norbitals ; i++ )
        {
            for ( j = 0 ; j <= i ; j++, ij++ )
            {
                for ( m = mn = 0 ; m < norbitals ; m++, mn++ )
                {
                    for ( n = 0 ; n < m ; n++, mn++ ) Array2D_Item ( cTransformation, ij, mn ) = Array2D_Item ( oTransformation, i, m ) * Array2D_Item ( oTransformation, j, n ) +
                                                                                                 Array2D_Item ( oTransformation, i, n ) * Array2D_Item ( oTransformation, j, m ) ;
                    Array2D_Item ( cTransformation, ij, mn ) = Array2D_Item ( oTransformation, i, m ) * Array2D_Item ( oTransformation, j, m ) ;
                }
            }
        }
        /* . Derivatives. */
        if ( doGradients )
        {
            /* . Allocation. */
            n = ( norbitals * ( norbitals + 1 ) ) / 2 ;
            cTransformationX = RealArray2D_AllocateWithExtents ( n, n, NULL ) ; RealArray2D_Set ( cTransformationX, 0.0e+00 ) ;
            cTransformationY = RealArray2D_AllocateWithExtents ( n, n, NULL ) ; RealArray2D_Set ( cTransformationY, 0.0e+00 ) ;
            cTransformationZ = RealArray2D_AllocateWithExtents ( n, n, NULL ) ; RealArray2D_Set ( cTransformationZ, 0.0e+00 ) ;
            /* . Form the matrices. */
            for ( i = ij = 0 ; i < norbitals ; i++ )
            {
                for ( j = 0 ; j <= i ; j++, ij++ )
                {
                    for ( m = mn = 0 ; m < norbitals ; m++, mn++ )
                    {
                        for ( n = 0 ; n < m ; n++, mn++ )
                        {
                            Array2D_Item ( cTransformationX, ij, mn ) = Array2D_Item ( oTransformationX, i, m ) * Array2D_Item ( oTransformation,  j, n ) +
                                                                        Array2D_Item ( oTransformation,  i, m ) * Array2D_Item ( oTransformationX, j, n ) +
                                                                        Array2D_Item ( oTransformationX, i, n ) * Array2D_Item ( oTransformation,  j, m ) +
                                                                        Array2D_Item ( oTransformation,  i, n ) * Array2D_Item ( oTransformationX, j, m ) ;
                            Array2D_Item ( cTransformationY, ij, mn ) = Array2D_Item ( oTransformationY, i, m ) * Array2D_Item ( oTransformation,  j, n ) +
                                                                        Array2D_Item ( oTransformation,  i, m ) * Array2D_Item ( oTransformationY, j, n ) +
                                                                        Array2D_Item ( oTransformationY, i, n ) * Array2D_Item ( oTransformation,  j, m ) +
                                                                        Array2D_Item ( oTransformation,  i, n ) * Array2D_Item ( oTransformationY, j, m ) ;
                            Array2D_Item ( cTransformationZ, ij, mn ) = Array2D_Item ( oTransformationZ, i, m ) * Array2D_Item ( oTransformation,  j, n ) +
                                                                        Array2D_Item ( oTransformation,  i, m ) * Array2D_Item ( oTransformationZ, j, n ) +
                                                                        Array2D_Item ( oTransformationZ, i, n ) * Array2D_Item ( oTransformation,  j, m ) +
                                                                        Array2D_Item ( oTransformation,  i, n ) * Array2D_Item ( oTransformationZ, j, m ) ;
                        }
                        Array2D_Item ( cTransformationX, ij, mn ) = Array2D_Item ( oTransformationX, i, m ) * Array2D_Item ( oTransformation,  j, m ) +
                                                                    Array2D_Item ( oTransformation,  i, m ) * Array2D_Item ( oTransformationX, j, m ) ;
                        Array2D_Item ( cTransformationY, ij, mn ) = Array2D_Item ( oTransformationY, i, m ) * Array2D_Item ( oTransformation,  j, m ) +
                                                                    Array2D_Item ( oTransformation,  i, m ) * Array2D_Item ( oTransformationY, j, m ) ;
                        Array2D_Item ( cTransformationZ, ij, mn ) = Array2D_Item ( oTransformationZ, i, m ) * Array2D_Item ( oTransformation,  j, m ) +
                                                                    Array2D_Item ( oTransformation,  i, m ) * Array2D_Item ( oTransformationZ, j, m ) ;
                    }
                }
            }
        }
        /* . Define the smaller matrix as a sub-block of the larger one. */
        /* . Should use slices or views here. */
        m = Minimum ( ni, nj ) ;
        if ( ( m > 1 ) && ( m < norbitals ) )
        {
            n = ( m * ( m + 1 ) ) / 2 ;
            clesser = RealArray2D_AllocateWithExtents ( n, n, NULL ) ;
            for ( i = 0 ; i < n ; i++ )
            {
                for ( j = 0 ; j < n ; j++ ) Array2D_Item ( clesser, i, j ) = Array2D_Item ( cTransformation, i, j ) ;
            }
            if ( doGradients )
            {
                clesserX = RealArray2D_AllocateWithExtents ( n, n, NULL ) ;
                clesserY = RealArray2D_AllocateWithExtents ( n, n, NULL ) ;
                clesserZ = RealArray2D_AllocateWithExtents ( n, n, NULL ) ;
                for ( i = 0 ; i < n ; i++ )
                {
                    for ( j = 0 ; j < n ; j++ )
                    {
                        Array2D_Item ( clesserX, i, j ) = Array2D_Item ( cTransformationX, i, j ) ;
                        Array2D_Item ( clesserY, i, j ) = Array2D_Item ( cTransformationY, i, j ) ;
                        Array2D_Item ( clesserZ, i, j ) = Array2D_Item ( cTransformationZ, i, j ) ;
                    }
                }
            }
        }
        /* . Assign the matrices to i and j.*/
        if ( ni == norbitals ) it = cTransformation ;
        else                   it = clesser         ;
        if ( nj == norbitals ) jt = cTransformation ;
        else                   jt = clesser         ;
        if ( doGradients )
        {
            if ( ni == norbitals ) { itX = cTransformationX ; itY = cTransformationY ; itZ = cTransformationZ ; }
            else                   { itX = clesserX         ; itY = clesserY         ; itZ = clesserZ         ; }
            if ( nj == norbitals ) { jtX = cTransformationX ; jtY = cTransformationY ; jtZ = cTransformationZ ; }
            else                   { jtX = clesserX         ; jtY = clesserY         ; jtZ = clesserZ         ; }
        }
        /* . Clear up. */
        RealArray2D_Deallocate ( &oTransformation ) ;
        RealArray2D_Deallocate ( &sTransformation ) ;
        if ( doGradients )
        {
            RealArray2D_Deallocate ( &oTransformationX ) ;
            RealArray2D_Deallocate ( &oTransformationY ) ;
            RealArray2D_Deallocate ( &oTransformationZ ) ;
        }
    }
    /* . Finish up. */
    if ( iTransformation  != NULL ) (*iTransformation ) = it  ;
    if ( jTransformation  != NULL ) (*jTransformation ) = jt  ;
    if ( iTransformationX != NULL ) (*iTransformationX) = itX ;
    if ( iTransformationY != NULL ) (*iTransformationY) = itY ;
    if ( iTransformationZ != NULL ) (*iTransformationZ) = itZ ;
    if ( jTransformationX != NULL ) (*jTransformationX) = jtX ;
    if ( jTransformationY != NULL ) (*jTransformationY) = jtY ;
    if ( jTransformationZ != NULL ) (*jTransformationZ) = jtZ ;
}
# undef _AlignmentTolerance

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate all two-center TEIs and optionally their derivatives in the local frame.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MNDOIntegralUtilities_LocalFrame2CTEIs ( const MNDOParameters *iData   ,
                                              const MNDOParameters *jData   ,
                                              const Real            r       ,
                                                    RealArray2D    *lfteis  ,
                                                    RealArray1D    *core1b  ,
                                                    RealArray1D    *core2a  ,
                                                    RealArray2D    *dlfteis ,
                                                    RealArray1D    *dcore1b ,
                                                    RealArray1D    *dcore2a )
{
    Boolean doGradients ;
# ifdef MNDODORBITALS
    const Integer  *negative = NULL, *positive = NULL, *unique = NULL ;
          Integer  c, i, iam = 0, ij, j, jam = 0, k, kl, l, nnegative, npositive, nunique, t ;
# else
    const Integer  *positive = NULL ;
          Integer  c, iam = 0, jam = 0, npositive, t ;
# endif

    /* . Check for gradients. */
    doGradients = ( dlfteis != NULL ) && ( dcore1b != NULL ) && ( dcore2a != NULL ) ;

    /* . Initialization. */
    RealArray2D_Set ( lfteis, 0.0e+00 ) ;
    RealArray1D_Set ( core1b, 0.0e+00 ) ;
    RealArray1D_Set ( core2a, 0.0e+00 ) ;
    if ( doGradients )
    {
        RealArray2D_Set ( dlfteis, 0.0e+00 ) ;
        RealArray1D_Set ( dcore1b, 0.0e+00 ) ;
        RealArray1D_Set ( dcore2a, 0.0e+00 ) ;
    }

    /* . Get the highest AM for each atom. */
    switch ( iData->norbitals )
    {
        case 1: iam = 0 ; break ;
        case 4: iam = 1 ; break ;
# ifdef MNDODORBITALS
        case 9: iam = 2 ; break ;
# endif
    }
    switch ( jData->norbitals )
    {
        case 1: jam = 0 ; break ;
        case 4: jam = 1 ; break ;
# ifdef MNDODORBITALS
        case 9: jam = 2 ; break ;
# endif
    }

    /* . Unique non-zero integrals in the local frame and then those related by symmetry. */
    /* . sp integrals. */
    MNDOIntegralUtilities_LocalFrame2CTEIsSP ( iData, jData       , r,        lfteis, dlfteis ) ;
    MNDOIntegralUtilities_LocalFrame2COEIsSP ( iData, jData->po[8], r, False, core1b, dcore1b ) ;
    MNDOIntegralUtilities_LocalFrame2COEIsSP ( jData, iData->po[8], r, True,  core2a, dcore2a ) ;

# ifdef MNDODORBITALS
    /* . Integrals involving d-orbitals. */
    /* . TEIs. */
    /* . Define the sets of integrals to evaluate. */
    nnegative = NNEGATIVE[iam][jam] ;
    npositive = NPOSITIVE[iam][jam] ;
    nunique   = NUNIQUE  [iam][jam] ;
    switch ( jam )
    {
        case 0:
            negative = SNEGATIVE ;
            positive = SPOSITIVE ;
            unique   = SUNIQUE   ;
            break ;
        case 1:
            negative = SPNEGATIVE ;
            positive = SPPOSITIVE ;
            unique   = SPUNIQUE   ;
            break ;
        case 2:
            negative = SPDNEGATIVE ;
            positive = SPDPOSITIVE ;
            unique   = SPDUNIQUE   ;
            break ;
    }

    /* . Integrals. */
    for ( c = t = 0 ; c < nunique ; c++, t += 4 )
    {
        i  = unique[t  ] ;
        j  = unique[t+1] ;
        k  = unique[t+2] ;
        l  = unique[t+3] ;
        ij = ( i * ( i + 1 ) ) / 2 + j ;
        kl = ( k * ( k + 1 ) ) / 2 + l ;
        Array2D_Item ( lfteis, ij, kl ) = MNDOIntegralUtilities_LocalFrame2CTEI ( &MNDOIntegralUtilities_2CChargeInteraction, iData, jData, ij, kl, ORBITALAM[i], ORBITALAM[j], ORBITALAM[k], ORBITALAM[l], 0, r ) ;
    }
    for ( c = t = 0 ; c < npositive ; c++, t+= 4 ) Array2D_Item ( lfteis, positive[t], positive[t+1] ) =   Array2D_Item ( lfteis, positive[t+2], positive[t+3] ) ;
    for ( c = t = 0 ; c < nnegative ; c++, t+= 4 ) Array2D_Item ( lfteis, negative[t], negative[t+1] ) = - Array2D_Item ( lfteis, negative[t+2], negative[t+3] ) ;

    /* . First atom electron - second atom core terms. */
    for ( c = t = 0 ; c < NCUNIQUE[iam] ; c++, t += 2 )
    {
        i  = CUNIQUE[t  ] ;
        j  = CUNIQUE[t+1] ;
        ij = ( i * ( i + 1 ) ) / 2 + j ;
        Array1D_Item ( core1b, ij ) = MNDOIntegralUtilities_LocalFrame2CTEI ( &MNDOIntegralUtilities_2CChargeInteraction, iData, jData, ij, SS, ORBITALAM[i], ORBITALAM[j], 0, 0, 2, r ) ;
    }
    for ( c = t = 0 ; c < NCPOSITIVE[iam] ; c++, t+= 2 ) Array1D_Item ( core1b, CPOSITIVE[t] ) =   Array1D_Item ( core1b, CPOSITIVE[t+1] ) ;
/*    for ( c = t = 0 ; c < NCNEGATIVE[iam] ; c++, t+= 2 ) Array1D_Item ( core1b, CNEGATIVE[t] ) = - Array1D_Item ( core1b, CNEGATIVE[t+1] ) ;*/
    RealArray1D_Scale ( core1b, - jData->zcore ) ;

    /* . Second atom electron - first atom core terms. */
    for ( c = t = 0 ; c < NCUNIQUE[jam] ; c++, t += 2 )
    {
        k  = CUNIQUE[t  ] ;
        l  = CUNIQUE[t+1] ;
        kl = ( k * ( k + 1 ) ) / 2 + l ;
        Array1D_Item ( core2a, kl ) = MNDOIntegralUtilities_LocalFrame2CTEI ( &MNDOIntegralUtilities_2CChargeInteraction, iData, jData, SS, kl, 0, 0, ORBITALAM[k], ORBITALAM[l], 1, r ) ;
    }
    for ( c = t = 0 ; c < NCPOSITIVE[jam] ; c++, t+= 2 ) Array1D_Item ( core2a, CPOSITIVE[t] ) =   Array1D_Item ( core2a, CPOSITIVE[t+1] ) ;
/*    for ( c = t = 0 ; c < NCNEGATIVE[jam] ; c++, t+= 2 ) Array1D_Item ( core2a, CNEGATIVE[t] ) = - Array1D_Item ( core2a, CNEGATIVE[t+1] ) ;*/
    RealArray1D_Scale ( core2a, - iData->zcore ) ;

    /* . Derivatives. */
    if ( doGradients )
    {
        for ( c = t = 0 ; c < nunique ; c++, t += 4 )
        {
            i  = unique[t  ] ;
            j  = unique[t+1] ;
            k  = unique[t+2] ;
            l  = unique[t+3] ;
            ij = ( i * ( i + 1 ) ) / 2 + j ;
            kl = ( k * ( k + 1 ) ) / 2 + l ;
            Array2D_Item ( dlfteis, ij, kl ) = MNDOIntegralUtilities_LocalFrame2CTEI ( &MNDOIntegralUtilities_2CChargeInteractionD, iData, jData, ij, kl, ORBITALAM[i], ORBITALAM[j], ORBITALAM[k], ORBITALAM[l], 0, r ) ;
        }
        for ( c = t = 0 ; c < npositive ; c++, t+= 4 ) Array2D_Item ( dlfteis, positive[t], positive[t+1] ) =   Array2D_Item ( dlfteis, positive[t+2], positive[t+3] ) ;
        for ( c = t = 0 ; c < nnegative ; c++, t+= 4 ) Array2D_Item ( dlfteis, negative[t], negative[t+1] ) = - Array2D_Item ( dlfteis, negative[t+2], negative[t+3] ) ;

        /* . First atom electron - second atom core terms. */
        for ( c = t = 0 ; c < NCUNIQUE[iam] ; c++, t += 2 )
        {
            i  = CUNIQUE[t  ] ;
            j  = CUNIQUE[t+1] ;
            ij = ( i * ( i + 1 ) ) / 2 + j ;
            Array1D_Item ( dcore1b, ij ) = MNDOIntegralUtilities_LocalFrame2CTEI ( &MNDOIntegralUtilities_2CChargeInteractionD, iData, jData, ij, SS, ORBITALAM[i], ORBITALAM[j], 0, 0, 2, r ) ;
        }
        for ( c = t = 0 ; c < NCPOSITIVE[iam] ; c++, t+= 2 ) Array1D_Item ( dcore1b, CPOSITIVE[t] ) =   Array1D_Item ( dcore1b, CPOSITIVE[t+1] ) ;
/*        for ( c = t = 0 ; c < NCNEGATIVE[iam] ; c++, t+= 2 ) Array1D_Item ( dcore1b, CNEGATIVE[t] ) = - Array1D_Item ( dcore1b, CNEGATIVE[t+1] ) ;*/
        RealArray1D_Scale ( dcore1b, - jData->zcore ) ;

        /* . Second atom electron - first atom core terms. */
        for ( c = t = 0 ; c < NCUNIQUE[jam] ; c++, t += 2 )
        {
            k  = CUNIQUE[t  ] ;
            l  = CUNIQUE[t+1] ;
            kl = ( k * ( k + 1 ) ) / 2 + l ;
            Array1D_Item ( dcore2a, kl ) = MNDOIntegralUtilities_LocalFrame2CTEI ( &MNDOIntegralUtilities_2CChargeInteractionD, iData, jData, SS, kl, 0, 0, ORBITALAM[k], ORBITALAM[l], 1, r ) ;
        }
        for ( c = t = 0 ; c < NCPOSITIVE[jam] ; c++, t+= 2 ) Array1D_Item ( dcore2a, CPOSITIVE[t] ) =   Array1D_Item ( dcore2a, CPOSITIVE[t+1] ) ;
/*        for ( c = t = 0 ; c < NCNEGATIVE[jam] ; c++, t+= 2 ) Array1D_Item ( dcore2a, CNEGATIVE[t] ) = - Array1D_Item ( dcore2a, CNEGATIVE[t+1] ) ;*/
        RealArray1D_Scale ( dcore2a, - iData->zcore ) ;
    }
# else
    /* . TEIs. */
    /* . Define the sets of integrals to evaluate. */
    npositive = NPOSITIVE[iam][jam] ;
    switch ( jam )
    {
        case 0:
            positive = SPOSITIVE ;
            break ;
        case 1:
            positive = SPPOSITIVE ;
            break ;
    }

    /* . Integrals. */
    for ( c = t = 0 ; c < npositive ; c++, t+= 4 ) Array2D_Item ( lfteis, positive[t], positive[t+1] ) =   Array2D_Item ( lfteis, positive[t+2], positive[t+3] ) ;

    /* . First atom electron - second atom core terms. */
    for ( c = t = 0 ; c < NCPOSITIVE[iam] ; c++, t+= 2 ) Array1D_Item ( core1b, CPOSITIVE[t] ) =   Array1D_Item ( core1b, CPOSITIVE[t+1] ) ;
    RealArray1D_Scale ( core1b, - jData->zcore ) ;

    /* . Second atom electron - first atom core terms. */
    for ( c = t = 0 ; c < NCPOSITIVE[jam] ; c++, t+= 2 ) Array1D_Item ( core2a, CPOSITIVE[t] ) =   Array1D_Item ( core2a, CPOSITIVE[t+1] ) ;
    RealArray1D_Scale ( core2a, - iData->zcore ) ;

    /* . Derivatives. */
    if ( doGradients )
    {
        for ( c = t = 0 ; c < npositive ; c++, t+= 4 ) Array2D_Item ( dlfteis, positive[t], positive[t+1] ) =   Array2D_Item ( dlfteis, positive[t+2], positive[t+3] ) ;

        /* . First atom electron - second atom core terms. */
        for ( c = t = 0 ; c < NCPOSITIVE[iam] ; c++, t+= 2 ) Array1D_Item ( dcore1b, CPOSITIVE[t] ) =   Array1D_Item ( dcore1b, CPOSITIVE[t+1] ) ;
        RealArray1D_Scale ( dcore1b, - jData->zcore ) ;

        /* . Second atom electron - first atom core terms. */
        for ( c = t = 0 ; c < NCPOSITIVE[jam] ; c++, t+= 2 ) Array1D_Item ( dcore2a, CPOSITIVE[t] ) =   Array1D_Item ( dcore2a, CPOSITIVE[t+1] ) ;
        RealArray1D_Scale ( dcore2a, - iData->zcore ) ;
    }
# endif
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the unique OEIs and TEIs involving sp orbitals in the local frame.
! . Derivatives (wrt r) can also be optionally calculated.
! . Unfortunately this code is incompatible with the d-orbital code (as the formulae are different).
! . This is identical to REPP with a few sign changes to cope with the change in transformation (1,4,7,8,9,12,13,14).
! . The OEIs and TEIs are treated separately.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define EV1 0.5e+00
# define EV2 0.25e+00
# define EV3 0.125e+00
# define EV4 0.0625e+00
# define PP  0.5e+00
# define TD  2.0e+00
void MNDOIntegralUtilities_LocalFrame2COEIsSP ( const MNDOParameters *iData   ,
                                                const Real            jPO8    ,
                                                const Real            r       ,
                                                const Boolean         swapped ,
                                                      RealArray1D    *core    ,
                                                      RealArray1D    *dcore   )
{
    Boolean doGradients ;
    Real    da, qa ;
    Real    ade, aee, aqe, dze, gdze, gqxxe, gqzze, gri0, gri2, gri3, qxxe, qzze, ri0, ri2, ri3, rsq, sqr1, sqr2, sqr3, sqr4, sqr5, sqr6, xxx ;

    /* . iData is the electronic atom, jPO8 is the core. */

    /* . Check for gradients. */
    doGradients = ( dcore != NULL ) ;

    /* . s/s - always done. */
    aee = pow ( ( iData->po[0] + jPO8 ), 2 ) ;
    rsq = r * r ;
    ri0 = 1.0e+00 / sqrt ( rsq + aee ) ;
    Array1D_Item ( core, SS ) = ri0 ;
    if ( doGradients )
    {
        gri0 = - r * ri0 * ri0 * ri0 ;
        Array1D_Item ( dcore, SS ) = gri0 ;
    }

    /* . sp/s. */
    /* . Redo ri0 with po[6]. */
    if ( iData->norbitals > 3 )
    {
        da   = iData->dd ;
        qa   = iData->qq * TD ;
        if ( swapped ) { da *= -1.0e+00 ; qa *= -1.0e+00 ; }
        aee  = pow ( ( iData->po[6] + jPO8 ), 2 ) ;
        ade  = pow ( ( iData->po[1] + jPO8 ), 2 ) ;
        aqe  = pow ( ( iData->po[2] + jPO8 ), 2 ) ;
        xxx  = r+da      ; sqr1 = 1.0e+00 / sqrt ( xxx*xxx + ade ) ;
        xxx  = r-da      ; sqr2 = 1.0e+00 / sqrt ( xxx*xxx + ade ) ;
        xxx  = r+qa      ; sqr3 = 1.0e+00 / sqrt ( xxx*xxx + aqe ) ;
        xxx  = r-qa      ; sqr4 = 1.0e+00 / sqrt ( xxx*xxx + aqe ) ;
        xxx  = rsq + aqe ; sqr5 = 1.0e+00 / sqrt ( xxx ) ;
                           sqr6 = 1.0e+00 / sqrt ( xxx + qa*qa ) ;
        dze  = - EV1 * ( sqr1 - sqr2 ) ;
        qzze = EV2 * ( sqr3 + sqr4 ) - EV1 * sqr5 ;
        qxxe = EV1 * ( sqr6 - sqr5 ) ;
        ri0  = 1.0e+00 / sqrt ( rsq + aee ) ;
        ri2  = ri0 + qzze ;
        ri3  = ri0 + qxxe ;
        Array1D_Item ( core, PZS  ) = dze ;
        Array1D_Item ( core, PZPZ ) = ri2 ;
        Array1D_Item ( core, PXPX ) = ri3 ;
        if ( doGradients )
        {
            xxx   = r * sqr5*sqr5*sqr5 ;
            gdze  = EV1 * ( (r+da) * sqr1*sqr1*sqr1 - (r-da) * sqr2*sqr2*sqr2 ) ;
            gqzze = - EV2 * ( (r+qa) * sqr3*sqr3*sqr3 + (r-qa) * sqr4*sqr4*sqr4 ) + EV1 * xxx ;
            gqxxe = - EV1 * ( r * sqr6*sqr6*sqr6 - xxx ) ;
            gri0  = - r * ri0 * ri0 * ri0 ;
            gri2  = gri0 + gqzze ;
            gri3  = gri0 + gqxxe ;
            Array1D_Item ( dcore, PZS  ) = gdze ;
            Array1D_Item ( dcore, PZPZ ) = gri2 ;
            Array1D_Item ( dcore, PXPX ) = gri3 ;
        }
    }
}

void MNDOIntegralUtilities_LocalFrame2CTEIsSP ( const MNDOParameters *iData   ,
                                                const MNDOParameters *jData   ,
                                                const Real            r       ,
                                                      RealArray2D    *lfteis  ,
                                                      RealArray2D    *dlfteis )
{
    Boolean doGradients ;
    Real    da = 0.0e+00, db = 0.0e+00, qa = 0.0e+00, qa0, qb = 0.0e+00, qb0 ;
    Real    ade, aee, aed, adq, aqe, aeq, aqd, aqq, arg35, arg38, arg39, axx, rsq, www, xxx, yyy, zzz,
            dze = 0.0e+00, edz = 0.0e+00, dxdx, dzdz, qxxe = 0.0e+00, eqxx = 0.0e+00, qzze = 0.0e+00, eqzz = 0.0e+00, dzqxx, qxxdz, dxqxz, qxzdx, dzqzz, qzzdz, qxxqxx, qxxqyy, qxxqzz, qzzqxx, qxzqxz, qzzqzz,
            gdze = 0.0e+00, gedz = 0.0e+00, gdxdx, gdzdz, gqxxe = 0.0e+00, geqxx = 0.0e+00, gqzze = 0.0e+00, geqzz = 0.0e+00, gdzqxx, gqxxdz, gdxqxz, gqxzdx, gdzqzz, gqzzdz, gqxxqxx, gqxxqyy, gqxxqzz, gqzzqxx, gqxzqxz, gqzzqzz,
            gri0 = 0.0e+00,  gri2, gri3, gri10, gri11, ri0, ri2, ri3, ri10, ri11,
            sqr1, sqr2, sqr3, sqr4, sqr5, sqr6, sqr7, sqr8, sqr9, sqr10, sqr11, sqr12, sqr13, sqr14, sqr15, sqr16, sqr17, sqr18, sqr19, sqr20,
            sqr21, sqr22, sqr23, sqr24, sqr25, sqr26, sqr27, sqr28, sqr29, sqr30, sqr31, sqr32, sqr33, sqr34, sqr35, sqr36, sqr37, sqr38, sqr39, sqr40,
            sqr41, sqr42, sqr43, sqr44, sqr45, sqr46, sqr47, sqr48, sqr49, sqr50, sqr51, sqr52, sqr53, sqr54, sqr55, sqr56, sqr57, sqr58, sqr59, sqr60,
/*            sqr61, sqr62, sqr63, */
            sqr64, sqr65, sqr66, sqr67, sqr68, sqr69, sqr70, sqr71 ;

    /* . Check for gradients. */
    doGradients = ( dlfteis != NULL ) ;

    /* . s/s - always done. */
    aee = pow ( ( iData->po[0] + jData->po[0] ), 2 ) ;
    rsq = r * r ;
    ri0 = 1.0e+00 / sqrt ( rsq + aee ) ;
    Array2D_Item ( lfteis, SS, SS ) = ri0 ;
    if ( doGradients )
    {
        gri0 = - r * ri0 * ri0 * ri0 ;
        Array2D_Item ( dlfteis, SS, SS ) = gri0 ;
    }

    /* . sp/s. */
    if ( iData->norbitals > 3 )
    {
        da   = iData->dd ;
        qa   = iData->qq * TD ;
        ade  = pow ( ( iData->po[1] + jData->po[0] ), 2 ) ;
        aqe  = pow ( ( iData->po[2] + jData->po[0] ), 2 ) ;
        xxx  = r+da      ; sqr1 = 1.0e+00 / sqrt ( xxx*xxx + ade ) ;
        xxx  = r-da      ; sqr2 = 1.0e+00 / sqrt ( xxx*xxx + ade ) ;
        xxx  = r+qa      ; sqr3 = 1.0e+00 / sqrt ( xxx*xxx + aqe ) ;
        xxx  = r-qa      ; sqr4 = 1.0e+00 / sqrt ( xxx*xxx + aqe ) ;
        xxx  = rsq + aqe ; sqr5 = 1.0e+00 / sqrt ( xxx ) ;
                           sqr6 = 1.0e+00 / sqrt ( xxx + qa*qa ) ;
        dze  = - EV1 * ( sqr1 - sqr2 ) ;
        qzze = EV2 * ( sqr3 + sqr4 ) - EV1 * sqr5 ;
        qxxe = EV1 * ( sqr6 - sqr5 ) ;
        ri2  = ri0 + qzze ;
        ri3  = ri0 + qxxe ;
        Array2D_Item ( lfteis, PZS , SS ) = dze ;
        Array2D_Item ( lfteis, PZPZ, SS ) = ri2 ;
        Array2D_Item ( lfteis, PXPX, SS ) = ri3 ;
        if ( doGradients )
        {
            xxx   = r * sqr5*sqr5*sqr5 ;
            gdze  = EV1 * ( (r+da) * sqr1*sqr1*sqr1 - (r-da) * sqr2*sqr2*sqr2 ) ;
            gqzze = - EV2 * ( (r+qa) * sqr3*sqr3*sqr3 + (r-qa) * sqr4*sqr4*sqr4 ) + EV1 * xxx ;
            gqxxe = - EV1 * ( r * sqr6*sqr6*sqr6 - xxx ) ;
            gri2  = gri0 + gqzze ;
            gri3  = gri0 + gqxxe ;
            Array2D_Item ( dlfteis, PZS , SS ) = gdze ;
            Array2D_Item ( dlfteis, PZPZ, SS ) = gri2 ;
            Array2D_Item ( dlfteis, PXPX, SS ) = gri3 ;
        }
    }

    /* . s/sp. */
    if ( jData->norbitals > 3 )
    {
        db    = jData->dd ;
        qb    = jData->qq * TD ;
        aed   = pow ( ( iData->po[0] + jData->po[1] ), 2 ) ;
        aeq   = pow ( ( iData->po[0] + jData->po[2] ), 2 ) ;
        xxx   = r-db      ; sqr7  = 1.0e+00 / sqrt ( xxx*xxx + aed ) ;
        xxx   = r+db      ; sqr8  = 1.0e+00 / sqrt ( xxx*xxx + aed ) ;
        xxx   = r-qb      ; sqr9  = 1.0e+00 / sqrt ( xxx*xxx + aeq ) ;
        xxx   = r+qb      ; sqr10 = 1.0e+00 / sqrt ( xxx*xxx + aeq ) ;
        xxx   = rsq + aeq ; sqr11 = 1.0e+00 / sqrt ( xxx ) ;
                            sqr12 = 1.0e+00 / sqrt ( xxx + qb*qb ) ;
        edz   = - EV1 * ( sqr7 - sqr8 ) ;
        eqzz  = EV2 * ( sqr9 + sqr10 ) - EV1 * sqr11 ;
        eqxx  = EV1 * ( sqr12 - sqr11 ) ;
        ri10  = ri0 + eqzz ;
        ri11  = ri0 + eqxx ;
        Array2D_Item ( lfteis, SS, PZS  ) = edz  ;
        Array2D_Item ( lfteis, SS, PZPZ ) = ri10 ;
        Array2D_Item ( lfteis, SS, PXPX ) = ri11 ;
        if ( doGradients )
        {
            xxx   = r * sqr11*sqr11*sqr11 ;
            gedz  = EV1 * ( (r-db) * sqr7*sqr7*sqr7 - (r+db) * sqr8*sqr8*sqr8 ) ;
            geqzz = - EV2 * ( (r-qb) * sqr9*sqr9*sqr9 + (r+qb) * sqr10*sqr10*sqr10 ) + EV1 * xxx ;
            geqxx = - EV1 * ( r * sqr12*sqr12*sqr12 - xxx ) ;
            gri10 = gri0 + geqzz ;
            gri11 = gri0 + geqxx ;
            Array2D_Item ( dlfteis, SS, PZS  ) = gedz  ;
            Array2D_Item ( dlfteis, SS, PZPZ ) = gri10 ;
            Array2D_Item ( dlfteis, SS, PXPX ) = gri11 ;
        }
    }

    /* . sp/sp. */
    if ( ( iData->norbitals > 3 ) && ( jData->norbitals > 3 ) )
    {
        axx   = pow ( ( iData->po[1] + jData->po[1] ), 2 ) ;
        adq   = pow ( ( iData->po[1] + jData->po[2] ), 2 ) ;
        aqd   = pow ( ( iData->po[2] + jData->po[1] ), 2 ) ;
        aqq   = pow ( ( iData->po[2] + jData->po[2] ), 2 ) ;
        xxx   = da - db         ; sqr13 = 1.0e+00 / sqrt ( rsq + axx + xxx * xxx ) ;
        xxx   = da + db         ; sqr14 = 1.0e+00 / sqrt ( rsq + axx + xxx * xxx ) ;
        xxx   = r + da - db     ; sqr15 = 1.0e+00 / sqrt ( xxx * xxx + axx ) ;
        xxx   = r - da + db     ; sqr16 = 1.0e+00 / sqrt ( xxx * xxx + axx ) ;
        xxx   = r - da - db     ; sqr17 = 1.0e+00 / sqrt ( xxx * xxx + axx ) ;
        xxx   = r + da + db     ; sqr18 = 1.0e+00 / sqrt ( xxx * xxx + axx ) ;
        xxx   = r + da          ;
        yyy   = xxx * xxx + adq ; sqr19 = 1.0e+00 / sqrt ( yyy ) ;
                                  sqr20 = 1.0e+00 / sqrt ( yyy + qb * qb ) ;
        xxx   = r - da          ;
        yyy   = xxx * xxx + adq ; sqr21 = 1.0e+00 / sqrt ( yyy ) ;
                                  sqr22 = 1.0e+00 / sqrt ( yyy + qb * qb ) ;
        xxx   = r - db          ;
        yyy   = xxx * xxx + aqd ; sqr23 = 1.0e+00 / sqrt ( yyy ) ;
                                  sqr24 = 1.0e+00 / sqrt ( yyy + qa * qa ) ;
        xxx   = r + db          ;
        yyy   = xxx * xxx + aqd ; sqr25 = 1.0e+00 / sqrt ( yyy ) ;
                                  sqr26 = 1.0e+00 / sqrt ( yyy + qa * qa ) ;
        xxx   = r + da - qb     ; sqr27 = 1.0e+00 / sqrt ( xxx * xxx + adq ) ;
        xxx   = r - da - qb     ; sqr28 = 1.0e+00 / sqrt ( xxx * xxx + adq ) ;
        xxx   = r + da + qb     ; sqr29 = 1.0e+00 / sqrt ( xxx * xxx + adq ) ;
        xxx   = r - da + qb     ; sqr30 = 1.0e+00 / sqrt ( xxx * xxx + adq ) ;
        xxx   = r + qa - db     ; sqr31 = 1.0e+00 / sqrt ( xxx * xxx + aqd ) ;
        xxx   = r + qa + db     ; sqr32 = 1.0e+00 / sqrt ( xxx * xxx + aqd ) ;
        xxx   = r - qa - db     ; sqr33 = 1.0e+00 / sqrt ( xxx * xxx + aqd ) ;
        xxx   = r - qa + db     ; sqr34 = 1.0e+00 / sqrt ( xxx * xxx + aqd ) ;
        arg35 = rsq + aqq       ; sqr35 = 1.0e+00 / sqrt ( arg35 ) ;
        xxx   = qa - qb         ; sqr36 = 1.0e+00 / sqrt ( arg35 + xxx * xxx ) ;
        xxx   = qa + qb         ; sqr37 = 1.0e+00 / sqrt ( arg35 + xxx * xxx ) ;
        xxx   = arg35 + qa * qa ; sqr38 = 1.0e+00 / sqrt ( xxx ) ;
                                  sqr39 = 1.0e+00 / sqrt ( arg35 + qb * qb ) ;
                                  sqr40 = 1.0e+00 / sqrt ( xxx + qb * qb ) ;
        xxx   = r - qb          ;
        yyy   = xxx * xxx + aqq ; sqr41 = 1.0e+00 / sqrt ( yyy ) ;
                                  sqr42 = 1.0e+00 / sqrt ( yyy + qa * qa ) ;
        xxx   = r + qb          ;
        yyy   = xxx * xxx + aqq ; sqr43 = 1.0e+00 / sqrt ( yyy ) ;
                                  sqr44 = 1.0e+00 / sqrt ( yyy + qa * qa ) ;
        xxx   = r + qa          ;
        yyy   = xxx * xxx + aqq ; sqr45 = 1.0e+00 / sqrt ( yyy ) ;
                                  sqr46 = 1.0e+00 / sqrt ( yyy + qb * qb ) ;
        xxx   = r - qa          ;
        yyy   = xxx * xxx + aqq ; sqr47 = 1.0e+00 / sqrt ( yyy ) ;
                                  sqr48 = 1.0e+00 / sqrt ( yyy + qb * qb ) ;
        xxx   = r + qa - qb     ; sqr49 = 1.0e+00 / sqrt ( xxx * xxx + aqq ) ;
        xxx   = r + qa + qb     ; sqr50 = 1.0e+00 / sqrt ( xxx * xxx + aqq ) ;
        xxx   = r - qa - qb     ; sqr51 = 1.0e+00 / sqrt ( xxx * xxx + aqq ) ;
        xxx   = r - qa + qb     ; sqr52 = 1.0e+00 / sqrt ( xxx * xxx + aqq ) ;
        qa0   = iData->qq ;
        qb0   = jData->qq ;
        xxx   = pow ( ( da - qb0 ), 2 ) ;
        yyy   = pow ( ( r  - qb0 ), 2 ) ;
        zzz   = pow ( ( da + qb0 ), 2 ) ;
        www   = pow ( ( r + qb0  ), 2 ) ;
        sqr53 = 1.0e+00 / sqrt ( xxx + yyy + adq ) ;
        sqr54 = 1.0e+00 / sqrt ( xxx + www + adq ) ;
        sqr55 = 1.0e+00 / sqrt ( zzz + yyy + adq ) ;
        sqr56 = 1.0e+00 / sqrt ( zzz + www + adq ) ;
        xxx   = pow ( ( qa0 - db ), 2 ) ;
        yyy   = pow ( ( qa0 + db ), 2 ) ;
        zzz   = pow ( ( r + qa0  ), 2 ) ;
        www   = pow ( ( r - qa0  ), 2 ) ;
        sqr57 = 1.0e+00 / sqrt ( zzz + xxx + aqd ) ;
        sqr58 = 1.0e+00 / sqrt ( www + xxx + aqd ) ;
        sqr59 = 1.0e+00 / sqrt ( zzz + yyy + aqd ) ;
        sqr60 = 1.0e+00 / sqrt ( www + yyy + aqd ) ;
        xxx   = pow ( ( qa0 - qb0 ), 2 ) ;
        yyy   = pow ( ( qa0 + qb0 ), 2 ) ;
/*
        sqr61 = 1.0e+00 / sqrt ( arg35 + TD * xxx ) ;
        sqr62 = 1.0e+00 / sqrt ( arg35 + TD * yyy ) ;
        sqr63 = 1.0e+00 / sqrt ( arg35 + TD * (qa0 * qa0 + qb0 * qb0) ) ;
*/
        zzz   = pow ( ( r + qa0 - qb0 ), 2 ) ;
        sqr64 = 1.0e+00 / sqrt ( zzz + xxx + aqq ) ;
        sqr65 = 1.0e+00 / sqrt ( zzz + yyy + aqq ) ;
        zzz   = pow ( ( r + qa0 + qb0 ), 2 ) ;
        sqr66 = 1.0e+00 / sqrt ( zzz + xxx + aqq ) ;
        sqr67 = 1.0e+00 / sqrt ( zzz + yyy + aqq ) ;
        zzz   = pow ( ( r - qa0 - qb0 ), 2 ) ;
        sqr68 = 1.0e+00 / sqrt ( zzz + xxx + aqq ) ;
        sqr69 = 1.0e+00 / sqrt ( zzz + yyy + aqq ) ;
        zzz   = pow ( ( r - qa0 + qb0 ), 2 ) ;
        sqr70 = 1.0e+00 / sqrt ( zzz + xxx + aqq ) ;
        sqr71 = 1.0e+00 / sqrt ( zzz + yyy + aqq ) ;
        dxdx   =   EV1 * sqr13 - EV1 * sqr14 ;
        dzdz   =   EV2 * sqr15 + EV2 * sqr16 - EV2 * sqr17 - EV2 * sqr18 ;
        dzqxx  =   EV2 * sqr19 - EV2 * sqr20 - EV2 * sqr21 + EV2 * sqr22 ;
        qxxdz  =   EV2 * sqr23 - EV2 * sqr24 - EV2 * sqr25 + EV2 * sqr26 ;
        dzqzz  = - EV3 * sqr27 + EV3 * sqr28 - EV3 * sqr29 + EV3 * sqr30 - EV2 * sqr21 + EV2 * sqr19 ;
        qzzdz  = - EV3 * sqr31 + EV3 * sqr32 - EV3 * sqr33 + EV3 * sqr34 + EV2 * sqr23 - EV2 * sqr25 ;
        qxxqxx =   EV3 * sqr36 + EV3 * sqr37 - EV2 * sqr38 - EV2 * sqr39 + EV2 * sqr35 ;
        qxxqyy =   EV2 * sqr40 - EV2 * sqr38 - EV2 * sqr39 + EV2 * sqr35 ;
        qxxqzz =   EV3 * sqr42 + EV3 * sqr44 - EV3 * sqr41 - EV3 * sqr43 - EV2 * sqr38 + EV2 * sqr35 ;
        qzzqxx =   EV3 * sqr46 + EV3 * sqr48 - EV3 * sqr45 - EV3 * sqr47 - EV2 * sqr39 + EV2 * sqr35 ;
        qzzqzz =   EV4 * sqr49 + EV4 * sqr50 + EV4 * sqr51 + EV4 * sqr52 - EV3 * sqr47 - EV3 * sqr45 - EV3 * sqr41 - EV3 * sqr43 + EV2 * sqr35 ;
        dxqxz  = - EV2 * sqr53 + EV2 * sqr54 + EV2 * sqr55 - EV2 * sqr56 ;
        qxzdx  = - EV2 * sqr57 + EV2 * sqr58 + EV2 * sqr59 - EV2 * sqr60 ;
        qxzqxz =   EV3 * sqr64 - EV3 * sqr66 - EV3 * sqr68 + EV3 * sqr70 - EV3 * sqr65 + EV3 * sqr67 + EV3 * sqr69 - EV3 * sqr71 ;
        Array2D_Item ( lfteis, PZS , PZS  ) = dzdz ;
        Array2D_Item ( lfteis, PXS , PXS  ) = dxdx ;
        Array2D_Item ( lfteis, PZPZ, PZS  ) = edz + qzzdz ;
        Array2D_Item ( lfteis, PXPX, PZS  ) = edz + qxxdz ;
        Array2D_Item ( lfteis, PXPZ, PXS  ) = qxzdx ;
        Array2D_Item ( lfteis, PZS , PZPZ ) = dze + dzqzz ;
        Array2D_Item ( lfteis, PZS , PXPX ) = dze + dzqxx ;
        Array2D_Item ( lfteis, PXS , PXPZ ) = dxqxz ;
        Array2D_Item ( lfteis, PZPZ, PZPZ ) = ri0 + eqzz + qzze + qzzqzz ;
        Array2D_Item ( lfteis, PXPX, PZPZ ) = ri0 + eqzz + qxxe + qxxqzz ;
        Array2D_Item ( lfteis, PZPZ, PXPX ) = ri0 + eqxx + qzze + qzzqxx ;
        Array2D_Item ( lfteis, PXPX, PXPX ) = ri0 + eqxx + qxxe + qxxqxx ;
        Array2D_Item ( lfteis, PXPZ, PXPZ ) = qxzqxz ;
        Array2D_Item ( lfteis, PXPX, PYPY ) = ri0 + eqxx + qxxe + qxxqyy ;
        Array2D_Item ( lfteis, PYPX, PYPX ) = PP * ( qxxqxx - qxxqyy ) ;
        if ( doGradients )
        {
            gdxdx   = - EV1 * r * ( sqr13*sqr13*sqr13 - sqr14*sqr14*sqr14 ) ;
            gdzdz   = - EV2 * ( (r+da-db)*sqr15*sqr15*sqr15 + (r-da+db)*sqr16*sqr16*sqr16 - (r-da-db)*sqr17*sqr17*sqr17 - (r+da+db)*sqr18*sqr18*sqr18 ) ;
            www     = (r+da)*sqr19*sqr19*sqr19 ;
            xxx     = (r-da)*sqr21*sqr21*sqr21 ;
            gdzqxx  = - EV2 * ( www - (r+da)*sqr20*sqr20*sqr20 - xxx + (r-da)*sqr22*sqr22*sqr22 ) ;
            yyy     = (r-db)*sqr23*sqr23*sqr23 ;
            zzz     = (r+db)*sqr25*sqr25*sqr25 ;
            gqxxdz  = - EV2 * ( yyy - (r-db)*sqr24*sqr24*sqr24 - zzz + (r+db)*sqr26*sqr26*sqr26 ) ;
            gdzqzz  = - EV3 * ( - (r+da-qb)*sqr27*sqr27*sqr27 + (r-da-qb)*sqr28*sqr28*sqr28 - (r+da+qb)*sqr29*sqr29*sqr29 + (r-da+qb)*sqr30*sqr30*sqr30 ) + EV2 * ( xxx - www ) ;
            gqzzdz  = - EV3 * ( - (r+qa-db)*sqr31*sqr31*sqr31 + (r+qa+db)*sqr32*sqr32*sqr32 - (r-qa-db)*sqr33*sqr33*sqr33 + (r-qa+db)*sqr34*sqr34*sqr34 ) - EV2 * ( yyy - zzz ) ;
            arg35   = r*sqr35*sqr35*sqr35 ;
            arg38   = r*sqr38*sqr38*sqr38 ;
            arg39   = r*sqr39*sqr39*sqr39 ;
            gqxxqxx = - EV3 * ( r*sqr36*sqr36*sqr36 + r*sqr37*sqr37*sqr37 ) + EV2 * ( arg38 + arg39 - arg35 ) ;
            gqxxqyy = - EV2 * ( r*sqr40*sqr40*sqr40 - arg38 - arg39 + arg35 ) ;
            www     = (r-qb)*sqr41*sqr41*sqr41 ;
            xxx     = (r+qb)*sqr43*sqr43*sqr43 ;
            yyy     = (r+qa)*sqr45*sqr45*sqr45 ;
            zzz     = (r-qa)*sqr47*sqr47*sqr47 ;
            gqxxqzz = - EV3 * ( (r-qb)*sqr42*sqr42*sqr42 + (r+qb)*sqr44*sqr44*sqr44 - www - xxx ) + EV2 * ( arg38 - arg35 ) ;
            gqzzqxx = - EV3 * ( (r+qa)*sqr46*sqr46*sqr46 + (r-qa)*sqr48*sqr48*sqr48 - yyy - zzz ) + EV2 * ( arg39 - arg35 ) ;
            gqzzqzz = - EV4 * ( (r+qa-qb)*sqr49*sqr49*sqr49 + (r+qa+qb)*sqr50*sqr50*sqr50 + (r-qa-qb)*sqr51*sqr51*sqr51 + (r-qa+qb)*sqr52*sqr52*sqr52 ) + EV3 * ( zzz + yyy + www + xxx ) - EV2 * arg35 ;
            gdxqxz  = - EV2 * ( - (r-qb0)*sqr53*sqr53*sqr53 + (r+qb0)*sqr54*sqr54*sqr54 + (r-qb0)*sqr55*sqr55*sqr55 - (r+qb0)*sqr56*sqr56*sqr56 ) ;
            gqxzdx  = - EV2 * ( - (r+qa0)*sqr57*sqr57*sqr57 + (r-qa0)*sqr58*sqr58*sqr58 + (r+qa0)*sqr59*sqr59*sqr59 - (r-qa0)*sqr60*sqr60*sqr60 ) ;
            gqxzqxz = - EV3 * ( (r+qa0-qb0)*sqr64*sqr64*sqr64 - (r+qa0+qb0)*sqr66*sqr66*sqr66 - (r-qa0-qb0)*sqr68*sqr68*sqr68 + (r-qa0+qb0)*sqr70*sqr70*sqr70 -
                                (r+qa0-qb0)*sqr65*sqr65*sqr65 + (r+qa0+qb0)*sqr67*sqr67*sqr67 + (r-qa0-qb0)*sqr69*sqr69*sqr69 - (r-qa0+qb0)*sqr71*sqr71*sqr71 ) ;
            Array2D_Item ( dlfteis, PZS , PZS  ) = gdzdz ;
            Array2D_Item ( dlfteis, PXS , PXS  ) = gdxdx ;
            Array2D_Item ( dlfteis, PZPZ, PZS  ) = gedz + gqzzdz ;
            Array2D_Item ( dlfteis, PXPX, PZS  ) = gedz + gqxxdz ;
            Array2D_Item ( dlfteis, PXPZ, PXS  ) = gqxzdx ;
            Array2D_Item ( dlfteis, PZS , PZPZ ) = gdze + gdzqzz ;
            Array2D_Item ( dlfteis, PZS , PXPX ) = gdze + gdzqxx ;
            Array2D_Item ( dlfteis, PXS , PXPZ ) = gdxqxz ;
            Array2D_Item ( dlfteis, PZPZ, PZPZ ) = gri0 + geqzz + gqzze + gqzzqzz ;
            Array2D_Item ( dlfteis, PXPX, PZPZ ) = gri0 + geqzz + gqxxe + gqxxqzz ;
            Array2D_Item ( dlfteis, PZPZ, PXPX ) = gri0 + geqxx + gqzze + gqzzqxx ;
            Array2D_Item ( dlfteis, PXPX, PXPX ) = gri0 + geqxx + gqxxe + gqxxqxx ;
            Array2D_Item ( dlfteis, PXPZ, PXPZ ) = gqxzqxz ;
            Array2D_Item ( dlfteis, PXPX, PYPY ) = gri0 + geqxx + gqxxe + gqxxqyy ;
            Array2D_Item ( dlfteis, PYPX, PYPX ) = PP * ( gqxxqxx - gqxxqyy ) ;
        }
    }
}
# undef EV1
# undef EV2
# undef EV3
# undef EV4
# undef PP
# undef TD

# ifdef MNDODORBITALS
/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate a two-center TEI or its derivative in the local frame.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . i, j, k, l - either 0, 1, or 2 with i >= j and k >= l.
! . c - either 0, 1 or 2.
*/
Real MNDOIntegralUtilities_LocalFrame2CTEI ( const ChargeInteractionFunction Evaluate ,
                                             const MNDOParameters           *iData    ,
                                             const MNDOParameters           *jData    ,
                                             const Integer                   ij       ,
                                             const Integer                   kl       ,
                                             const Integer                   i        ,
                                             const Integer                   j        ,
                                             const Integer                   k        ,
                                             const Integer                   l        ,
                                             const Integer                   c        ,
                                             const Real                      r        )
{
    Real integral = 0.0e+00 ;
    if ( ( iData != NULL ) && ( jData != NULL ) && ( NCHTERMS[ij] > 0 ) && ( NCHTERMS[kl] > 0 ) )
    {
        auto Real add, chijkl[MAXCHTERMS], dij = 0.0e+00, dkl = 0.0e+00, pij = 0.0e+00, pkl = 0.0e+00 ;
        auto Integer     lij, lkl, lmin, lm1, lm2, l1, l1max, l1min, l1offset, l1terms[MAXCHTERMS], l2, l2max, l2min, l2offset, l2terms[MAXCHTERMS], m, mterms[MAXCHTERMS], nterms, t ;
        /* . Get loop indices. */
        /* . Possibilites for ( i, j, l1min, l1max ) are ( 0, 0, 0, 0 ), ( 1, 0, 1, 1 ), ( 1, 1, 0, 2 ), ( 2, 0, 2, 2 ), ( 2, 1, 1, 2 ) and ( 2, 2, 0, 2 ). */
        /* . Special cases (l1 or l2 == 0) are the diagonal ones; i.e. i == j with i = 0, 1 or 2. */
        l1min = i - j ;
        l1max = Minimum ( i + j, 2 ) ;
        lij   = ( i * ( i + 1 ) ) / 2 + j ;
        l2min = k - l ;
        l2max = Minimum ( k + l, 2 ) ;
        lkl   = ( k * ( k + 1 ) ) / 2 + l ;
        /* . Preprocessing for number of terms. */
        for ( l1 = l1min, nterms = 0 ; l1 <= l1max ; l1++ )
        {
            l1offset = ij * CHINCREMENT1 + l1 * CHINCREMENT2 + CHINCREMENT3 ;
            for ( l2 = l2min ; l2 <= l2max ; l2++ )
            {
                l2offset = kl * CHINCREMENT1 + l2 * CHINCREMENT2 + CHINCREMENT3 ;
                lmin     = Minimum ( l1, l2 ) ;
                for ( m = -lmin ; m <= lmin ; m++ )
                {
                    lm1 = CHINDICES[l1offset+m] ;
                    lm2 = CHINDICES[l2offset+m] ;
                    if ( ( lm1 > 0 ) && ( lm2 > 0 ) )
                    {
                        l1terms[nterms] = l1 ;
                        l2terms[nterms] = l2 ;
                        mterms [nterms] = abs ( m ) ;
                        chijkl [nterms] = CHTERMS[lm1-1] * CHTERMS[lm2-1] ;
                        nterms ++ ;
                    }
                }
            }
        }
        /* . Calculate the terms. */
        for ( t = 0, integral = 0.0e+00 ; t < nterms ; t++ )
        {
            l1 = l1terms[t] ;
            l2 = l2terms[t] ;
            m  = mterms [t] ;
            if ( l1 == 0 )
            {
                dij = 0.0e+00 ;
                switch ( i )
                {
                    case 0:
                        pij = iData->po[0] ;
                        if ( c == 1 ) pij = iData->po[8] ;
                        break ;
                    case 1:
                        pij = iData->po[6] ;
                        break ;
                    case 2:
                        pij = iData->po[7] ;
                        break ;
                }
            }
            else
            {
                dij = iData->ddp[lij] ;
                pij = iData->po [lij] ;
            }
            if ( l2 == 0 )
            {
                dkl = 0.0e+00 ;
                switch ( k )
                {
                    case 0:
                        pkl = jData->po[0] ;
                        if ( c == 2 ) pkl = jData->po[8] ;
                        break ;
                    case 1:
                        pkl = jData->po[6] ;
                        break ;
                    case 2:
                        pkl = jData->po[7] ;
                        break ;
                }
            }
            else
            {
                dkl = jData->ddp[lkl] ;
                pkl = jData->po [lkl] ;
            }
            add = pow ( ( pij + pkl ), 2 ) ;
            integral += chijkl[t] * Evaluate ( r, l1, l2, m, dij, dkl, add ) ;
        }
    }
    return integral ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the interaction of two point-charge configurations.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . r      distance.
! . l1,m   quantum numbers for multipole of configuration 1.
! . l2,m   quantum numbers for multipole of configuration 2.
! . da     charge separation of configuration 1.
! . db     charge separation of configuration 2.
! . add    additive term.
*/
Real MNDOIntegralUtilities_2CChargeInteraction ( const Real r, const Integer  l1, const Integer  l2, const Integer  m, const Real da, const Real db, const Real add )
{
    Real aa, ab, charg, dxdx, dxqxz, dzdz, dzqzz, qqzz, qxzdx, qxzqxz, qzzdz, qzzq, xyxy, zzzz ;
    charg = 0.0e+00 ;
    switch ( l1 )
    {
        case 0:
            switch ( l2 )
            {
                case 0:
                    charg = 1.0e+00 / sqrt( r*r + add ) ;
                break ;
                case 1:
                    charg = 1.0e+00 / sqrt(pow(r + db,2) + add)  - 1.0e+00/sqrt(pow(r - db,2) + add) ; charg /= 2.0e+00 ;
                break ;
                case 2:
                    qqzz  = 1.0e+00/sqrt(pow(r - db,2) + add) - 2.0e+00/sqrt(r*r + db*db + add) + 1.0e+00/sqrt(pow(r + db,2) + add) ;
                    charg = qqzz/4.0e+00 ;
                break ;
            }
            break ;
        case 1:
            switch ( l2 )
            {
                case 0:
                    charg = (-1.0e+00/sqrt(pow(r + da,2) + add)) + 1.0e+00/sqrt(pow(r - da,2) + add) ; charg /= 2.0e+00 ;
                break ;
                case 1:
                if ( m == 0 )
                {
                    dzdz  = 1.0e+00/sqrt(pow(r + da - db,2) + add) + 1.0e+00/sqrt(pow(r - da + db,2) + add) - 1.0e+00/sqrt(pow(r - da - db,2) + add) - 1.0e+00/sqrt(pow(r + da + db,2) + add) ;
                    charg = dzdz/4.0e+00 ;
                }
                else if ( m == 1 )
                {
                    dxdx  = 2.0e+00/sqrt(r*r + pow(da - db,2) + add) - 2.0e+00/sqrt(r*r + pow(da + db,2) + add) ;
                    charg = dxdx/4.0e+00 ;
                }
                break ;
                case 2:
                if ( m == 0 )
                {
                    dzqzz = 1.0e+00/sqrt(pow(r - da - db,2) + add) - 2.0e+00/sqrt(pow(r - da,2) + db*db + add) + 1.0e+00/sqrt(pow(r + db - da,2) + add) -
                            1.0e+00/sqrt(pow(r - db + da,2) + add) + 2.0e+00/sqrt(pow(r + da,2) + db*db + add) - 1.0e+00/sqrt(pow(r + da + db,2) + add) ;
                    charg = dzqzz/8.0e+00 ;
                }
                else if ( m == 1 )
                {
                    ab = db/sqrt(2.0e+00) ;
                    dxqxz = (-2.0e+00/sqrt(pow(r - ab,2) + pow(da - ab,2) + add)) + 2.0e+00/sqrt(pow(r + ab,2) + pow(da - ab,2) + add) + 2.0e+00/sqrt(pow(r - ab,2) + pow(da + ab,2) + add) - 2.0e+00/sqrt(pow(r + ab,2) + pow(da + ab,2) + add) ;
                    charg = dxqxz/8.0e+00 ;
                }
                break ;
            }
            break ;
        case 2:
            switch ( l2 )
            {
                case 0:
                    qzzq  = 1.0e+00/sqrt(pow(r - da,2) + add) - 2.0e+00/sqrt(r*r + da*da + add) + 1.0e+00/sqrt(pow(r + da,2) + add) ;
                    charg = qzzq/4.0e+00 ;
                break ;
                case 1:
                if ( m == 0 )
                {
                    qzzdz = (-1.0e+00/sqrt(pow(r - da - db,2) + add)) + 2.0e+00/sqrt(pow(r - db,2) + da*da + add) - 1.0e+00/sqrt(pow(r + da - db,2) + add) +
                              1.0e+00/sqrt(pow(r - da + db,2) + add) - 2.0e+00/sqrt(pow(r + db,2) + da*da + add) + 1.0e+00/sqrt(pow(r + da + db,2) + add) ;
                    charg = qzzdz/8.0e+00 ;
                }
                else if ( m == 1 )
                {
                    aa = da/sqrt(2.0e+00) ;
                    qxzdx = (-2.0e+00/sqrt(pow(r + aa,2) + pow(aa - db,2) + add)) + 2.0e+00/sqrt(pow(r - aa,2) + pow(aa - db,2) + add) + 2.0e+00/sqrt(pow(r + aa,2) + pow(aa + db,2) + add) - 2.0e+00/sqrt(pow(r - aa,2) + pow(aa + db,2) + add) ;
                    charg = qxzdx/8.0e+00 ;
                }
                break ;
                case 2:
                if ( m == 0 )
                {
                    zzzz  = 1.0e+00/sqrt(pow(r - da - db,2) + add) + 1.0e+00/sqrt(pow(r + da + db,2) + add) + 1.0e+00/sqrt(pow(r - da + db,2) + add) + 1.0e+00/sqrt(pow(r + da - db,2) + add) -
                                                          2.0e+00/sqrt(pow(r - da,2) + db*db + add) - 2.0e+00/sqrt(pow(r - db,2) + da*da + add) - 2.0e+00/sqrt(pow(r + da,2) + db*db + add) -
                                                          2.0e+00/sqrt(pow(r + db,2) + da*da + add) + 2.0e+00/sqrt(r*r + pow(da - db,2) + add) + 2.0e+00/sqrt(r*r + pow(da + db,2) + add) ;
                    xyxy  = 4.0e+00/sqrt(r*r + pow(da - db,2) + add) + 4.0e+00/sqrt(r*r + pow(da + db,2) + add) - 8.0e+00/sqrt(r*r + da*da + db*db + add) ;
                    charg = zzzz/16.0e+00 - xyxy/64.0e+00 ;
                }
                else if ( m == 1 )
                {
                    aa = da/sqrt(2.0e+00) ;
                    ab = db/sqrt(2.0e+00) ;
                    qxzqxz = 2.0e+00/sqrt(pow(r + aa - ab,2) + pow(aa - ab,2) + add) - 2.0e+00/sqrt(pow(r + aa + ab,2) + pow(aa - ab,2) + add) - 2.0e+00/sqrt(pow(r - aa - ab,2) + pow(aa - ab,2) + add) +
                             2.0e+00/sqrt(pow(r - aa + ab,2) + pow(aa - ab,2) + add) - 2.0e+00/sqrt(pow(r + aa - ab,2) + pow(aa + ab,2) + add) + 2.0e+00/sqrt(pow(r + aa + ab,2) + pow(aa + ab,2) + add) +
                             2.0e+00/sqrt(pow(r - aa - ab,2) + pow(aa + ab,2) + add) - 2.0e+00/sqrt(pow(r - aa + ab,2) + pow(aa + ab,2) + add) ;
                    charg  = qxzqxz/16.0e+00 ;
                }
                else if ( m == 2 )
                {
                    xyxy  = 4.0e+00/sqrt(r*r + pow(da - db,2) + add) + 4.0e+00/sqrt(r*r + pow(da + db,2) + add) - 8.0e+00/sqrt(r*r + da*da + db*db + add) ;
                    charg = xyxy/16.0e+00 ;
                }
                break ;
            }
            break ;
    }
    return charg ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the derivative of the interaction of two point-charge configurations.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . r      distance.
! . l1,m   quantum numbers for multipole of configuration 1.
! . l2,m   quantum numbers for multipole of configuration 2.
! . da     charge separation of configuration 1.
! . db     charge separation of configuration 2.
! . add    additive term.
*/
Real MNDOIntegralUtilities_2CChargeInteractionD ( const Real r, const Integer  l1, const Integer  l2, const Integer  m, const Real da, const Real db, const Real add )
{
    Real aa, ab, dcharg, dxdx, dxqxz, dzdz, dzqzz, fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8, fac9, fac10, qqzz, qxzdx, qxzqxz, qzzdz, qzzq, xyxy, zzzz ;
    dcharg = 0.0e+00 ;
    switch ( l1 )
    {
        case 0:
            switch ( l2 )
            {
                case 0:
                    fac1   = r*r + add ;
                    dcharg = - r / ( fac1 * sqrt ( fac1 ) ) ;
                break ;
                case 1:
                    fac1    = pow(r + db,2) + add ;
                    fac2    = pow(r - db,2) + add ;
                    dcharg  = ( r + db ) / ( fac1 * sqrt ( fac1 ) ) - ( r - db ) / ( fac2 * sqrt ( fac2 ) ) ;
                    dcharg *= - 0.5e+00 ;
                break ;
                case 2:
                    fac1   = pow(r - db,2) + add ;
                    fac2   = r*r + db*db + add ;
                    fac3   = pow(r + db,2) + add ;
                    qqzz   = ( r - db ) / ( fac1 * sqrt ( fac1 ) ) - 2.0e+00 * r / ( fac2 * sqrt ( fac2 ) ) + ( r + db ) / ( fac3 * sqrt ( fac3 ) ) ;
                    dcharg = - qqzz / 4.0e+00 ;
                break ;
            }
            break ;
        case 1:
            switch ( l2 )
            {
                case 0:
                    fac1    = pow(r + da,2) + add ;
                    fac2    = pow(r - da,2) + add ;
                    dcharg  = - ( r + da ) / ( fac1 * sqrt ( fac1 ) ) + ( r - da ) / ( fac2 * sqrt ( fac2 ) ) ;
                    dcharg *= - 0.5e+00 ;
                break ;
                case 1:
                if ( m == 0 )
                {
                    fac1   = pow(r + da - db,2) + add ;
                    fac2   = pow(r - da + db,2) + add ;
                    fac3   = pow(r - da - db,2) + add ;
                    fac4   = pow(r + da + db,2) + add ;
                    dzdz   = ( r + da - db ) / ( fac1 * sqrt ( fac1 ) ) + ( r - da + db ) / ( fac2 * sqrt ( fac2 ) ) - ( r - da - db ) / ( fac3 * sqrt ( fac3 ) ) - ( r + da + db ) / ( fac4 * sqrt ( fac4 ) ) ;
                    dcharg = - dzdz / 4.0e+00 ;
                }
                else if ( m == 1 )
                {
                    fac1   = r*r + pow(da - db,2) + add ;
                    fac2   = r*r + pow(da + db,2) + add ;
                    dxdx   = 2.0e+00 * r / ( fac1 * sqrt ( fac1 ) ) - 2.0e+00 * r / ( fac2 * sqrt ( fac2 ) ) ;
                    dcharg = - dxdx / 4.0e+00 ;
                }
                break ;
                case 2:
                if ( m == 0 )
                {
                    fac1   = pow(r - da - db,2) + add ;
                    fac2   = pow(r - da,2) + db*db + add ;
                    fac3   = pow(r + db - da,2) + add ;
                    fac4   = pow(r - db + da,2) + add ;
                    fac5   = pow(r + da,2) + db*db + add ;
                    fac6   = pow(r + da + db,2) + add ;
                    dzqzz  = ( r - da - db ) / ( fac1 * sqrt ( fac1 ) ) - 2.0e+00 * ( r - da ) / ( fac2 * sqrt ( fac2 ) ) + ( r + db - da ) / ( fac3 * sqrt ( fac3 ) ) -
                             ( r - db + da ) / ( fac4 * sqrt ( fac4 ) ) + 2.0e+00 * ( r + da ) / ( fac5 * sqrt ( fac5 ) ) - ( r + da + db ) / ( fac6 * sqrt ( fac6 ) ) ;
                    dcharg = - dzqzz / 8.0e+00 ;
                }
                else if ( m == 1 )
                {
                    ab     = db / sqrt ( 2.0e+00 ) ;
                    fac1   = pow(r - ab,2) + pow(da - ab,2) + add ;
                    fac2   = pow(r + ab,2) + pow(da - ab,2) + add ;
                    fac3   = pow(r - ab,2) + pow(da + ab,2) + add ;
                    fac4   = pow(r + ab,2) + pow(da + ab,2) + add ;
                    dxqxz  = - 2.0e+00 * ( r - ab ) / ( fac1 * sqrt ( fac1 ) ) + 2.0e+00 * ( r + ab ) / ( fac2 * sqrt ( fac2 ) ) + 2.0e+00 * ( r - ab ) / ( fac3 * sqrt ( fac3 ) ) - 2.0e+00 * ( r + ab ) / ( fac4 * sqrt ( fac4 ) ) ;
                    dcharg = - dxqxz / 8.0e+00 ;
                }
                break ;
            }
            break ;
        case 2:
            switch ( l2 )
            {
                case 0:
                    fac1   = pow(r - da,2) + add ;
                    fac2   = r*r + da*da   + add ;
                    fac3   = pow(r + da,2) + add ;
                    qzzq   = ( r - da ) / ( fac1 * sqrt ( fac1 ) ) - 2.0e+00 * r / ( fac2 * sqrt ( fac2 ) ) + ( r + da ) / ( fac3 * sqrt ( fac3 ) ) ;
                    dcharg = - qzzq / 4.0e+00 ;
                break ;
                case 1:
                if ( m == 0 )
                {
                    fac1   = pow(r - da - db,2) + add ;
                    fac2   = pow(r - db,2) + da*da + add ;
                    fac3   = pow(r + da - db,2) + add ;
                    fac4   = pow(r - da + db,2) + add ;
                    fac5   = pow(r + db,2) + da*da + add ;
                    fac6   = pow(r + da + db,2) + add ;
                    qzzdz  = - ( r - da - db ) / ( fac1 * sqrt ( fac1 ) ) + 2.0e+00 * ( r - db ) / ( fac2 * sqrt ( fac2 ) ) - ( r + da - db ) / ( fac3 * sqrt ( fac3 ) ) +
                               ( r - da + db ) / ( fac4 * sqrt ( fac4 ) ) - 2.0e+00 * ( r + db ) / ( fac5 * sqrt ( fac5 ) ) + ( r + da + db ) / ( fac6 * sqrt ( fac6 ) ) ;
                    dcharg = - qzzdz / 8.0e+00 ;
                }
                else if ( m == 1 )
                {
                    aa     = da / sqrt ( 2.0e+00 ) ;
                    fac1   = pow(r + aa,2) + pow(aa - db,2) + add ;
                    fac2   = pow(r - aa,2) + pow(aa - db,2) + add ;
                    fac3   = pow(r + aa,2) + pow(aa + db,2) + add ;
                    fac4   = pow(r - aa,2) + pow(aa + db,2) + add ;
                    qxzdx  = -2.0e+00 * ( r + aa ) / ( fac1 * sqrt ( fac1 ) ) + 2.0e+00 * ( r - aa ) / ( fac2 * sqrt ( fac2 ) ) + 2.0e+00 * ( r + aa ) / ( fac3 * sqrt ( fac3 ) ) - 2.0e+00 * ( r - aa ) / ( fac4 * sqrt ( fac4 ) ) ;
                    dcharg = - qxzdx / 8.0e+00 ;
                }
                break ;
                case 2:
                if ( m == 0 )
                {
                    fac1   = pow(r - da - db,2) + add ;
                    fac2   = pow(r + da + db,2) + add ;
                    fac3   = pow(r - da + db,2) + add ;
                    fac4   = pow(r + da - db,2) + add ;
                    fac5   = pow(r - da,2) + db*db + add ;
                    fac6   = pow(r - db,2) + da*da + add ;
                    fac7   = pow(r + da,2) + db*db + add ;
                    fac8   = pow(r + db,2) + da*da + add ;
                    fac9   = r*r + pow(da - db,2) + add ;
                    fac10  = r*r + pow(da + db,2) + add ;
                    zzzz   = ( r - da - db ) / ( fac1 * sqrt ( fac1 ) ) + ( r + da + db ) / ( fac2 * sqrt ( fac2 ) ) +
                             ( r - da + db ) / ( fac3 * sqrt ( fac3 ) ) + ( r + da - db ) / ( fac4 * sqrt ( fac4 ) ) -
                             2.0e+00 * ( r - da ) / ( fac5 * sqrt ( fac5 ) ) - 2.0e+00 * ( r - db ) / ( fac6 * sqrt ( fac6 ) ) -
                             2.0e+00 * ( r + da ) / ( fac7 * sqrt ( fac7 ) ) - 2.0e+00 * ( r + db ) / ( fac8 * sqrt ( fac8 ) ) +
                             2.0e+00 * r / ( fac9 * sqrt ( fac9 ) )          + 2.0e+00 * r / ( fac10 * sqrt ( fac10 ) ) ;
                    fac1   = r*r + pow(da - db,2) + add ;
                    fac2   = r*r + pow(da + db,2) + add ;
                    fac3   = r*r + da*da + db*db  + add ;
                    xyxy   = 4.0e+00 * r / ( fac1 * sqrt ( fac1 ) ) + 4.0e+00 * r / ( fac2 * sqrt ( fac2 ) ) - 8.0e+00 * r / ( fac3 * sqrt ( fac3 ) ) ;
                    dcharg = - ( zzzz / 16.0e+00 - xyxy / 64.0e+00 ) ;
                }
                else if ( m == 1 )
                {
                    aa      = da / sqrt ( 2.0e+00 ) ;
                    ab      = db / sqrt ( 2.0e+00 ) ;
                    fac1    = pow(r + aa - ab,2) + pow(aa - ab,2) + add ;
                    fac2    = pow(r + aa + ab,2) + pow(aa - ab,2) + add ;
                    fac3    = pow(r - aa - ab,2) + pow(aa - ab,2) + add ;
                    fac4    = pow(r - aa + ab,2) + pow(aa - ab,2) + add ;
                    fac5    = pow(r + aa - ab,2) + pow(aa + ab,2) + add ;
                    fac6    = pow(r + aa + ab,2) + pow(aa + ab,2) + add ;
                    fac7    = pow(r - aa - ab,2) + pow(aa + ab,2) + add ;
                    fac8    = pow(r - aa + ab,2) + pow(aa + ab,2) + add ;
                    qxzqxz  = 2.0e+00 * ( r + aa - ab ) / ( fac1 * sqrt ( fac1 ) ) - 2.0e+00 * ( r + aa + ab ) / ( fac2 * sqrt ( fac2 ) ) - 2.0e+00 * ( r - aa - ab ) / ( fac3 * sqrt ( fac3 ) ) +
                              2.0e+00 * ( r - aa + ab ) / ( fac4 * sqrt ( fac4 ) ) - 2.0e+00 * ( r + aa - ab ) / ( fac5 * sqrt ( fac5 ) ) + 2.0e+00 * ( r + aa + ab ) / ( fac6 * sqrt ( fac6 ) ) +
                              2.0e+00 * ( r - aa - ab ) / ( fac7 * sqrt ( fac7 ) ) - 2.0e+00 * ( r - aa + ab ) / ( fac8 * sqrt ( fac8 ) ) ;
                    dcharg  = - qxzqxz / 16.0e+00 ;
                }
                else if ( m == 2 )
                {
                    fac1   = r*r + pow(da - db,2) + add ;
                    fac2   = r*r + pow(da + db,2) + add ;
                    fac3   = r*r + da*da + db*db  + add ;
                    xyxy   = 4.0e+00 * r / ( fac1 * sqrt ( fac1 ) ) + 4.0e+00 * r / ( fac2 * sqrt ( fac2 ) ) - 8.0e+00 * r / ( fac3 * sqrt ( fac3 ) ) ;
                    dcharg = - xyxy / 16.0e+00 ;
                }
                break ;
            }
            break ;
    }
    return dcharg ;
}
# endif

# undef MNDODORBITALS

