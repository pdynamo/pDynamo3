/*==================================================================================================================================
! . Procedures for 3 x 3 matrices.
!=================================================================================================================================*/

# include <math.h>

# include "Matrix33.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
Matrix33 *Matrix33_Allocate ( void ) { return RealArray2D_AllocateWithExtents ( 3, 3, NULL ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Apply to a vector3.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Matrix33_ApplyToVector3 ( const Matrix33 *self, Vector3 *vector3 )
{
    if ( ( self != NULL ) && ( vector3 != NULL ) )
    {
        auto Real x, y, z ;
        x = Vector3_Item ( vector3, 0 ) ;
        y = Vector3_Item ( vector3, 1 ) ;
        z = Vector3_Item ( vector3, 2 ) ;
        Vector3_Item ( vector3, 0 ) = x * Matrix33_Item ( self, 0, 0 ) +
                                      y * Matrix33_Item ( self, 0, 1 ) +
                                      z * Matrix33_Item ( self, 0, 2 ) ;
        Vector3_Item ( vector3, 1 ) = x * Matrix33_Item ( self, 1, 0 ) +
                                      y * Matrix33_Item ( self, 1, 1 ) +
                                      z * Matrix33_Item ( self, 1, 2 ) ;
        Vector3_Item ( vector3, 2 ) = x * Matrix33_Item ( self, 2, 0 ) +
                                      y * Matrix33_Item ( self, 2, 1 ) +
                                      z * Matrix33_Item ( self, 2, 2 ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Find the angle and axis for a given rotation.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Use Rodriguez's identity: R = I + N * sinTheta + N^2 * ( 1 - cosTheta ) with N skew matrix of axis vector elements. */

# define Tolerance 1.0e-10 /* . Should have separate tolerances? */
void Matrix33_AngleAxisFromRotation ( const Matrix33 *self, const Real *tolerance, Real *angle, Vector3 *axis, Status *status )
{
    if ( ( self != NULL ) && ( axis != NULL ) )
    {
        /* . Is the matrix orthogonal? */
        if ( Matrix33_IsOrthogonal ( self ) )
        {
            auto Real cosTheta, sinTheta, theta, tol ;

            /* . Get tolerances. */
            if ( tolerance == NULL ) tol = Tolerance ;
            else                     tol = (*tolerance)   ;

            /* . Initial guess at the axis. */
            Vector3_Item ( axis, 0 ) = Matrix33_Item ( self, 2, 1 ) - Matrix33_Item ( self, 1, 2 ) ;
            Vector3_Item ( axis, 1 ) = Matrix33_Item ( self, 0, 2 ) - Matrix33_Item ( self, 2, 0 ) ;
            Vector3_Item ( axis, 2 ) = Matrix33_Item ( self, 1, 0 ) - Matrix33_Item ( self, 0, 1 ) ;
            Vector3_Scale ( axis, 0.5e+00 ) ;

            /* . Get the angle data. */
            cosTheta = 0.5e+00 * ( Matrix33_Item ( self, 0, 0 ) +
                                   Matrix33_Item ( self, 1, 1 ) +
                                   Matrix33_Item ( self, 2, 2 ) - 1.0e+00 ) ;
            sinTheta = Vector3_Norm2 ( axis ) ;

            /* . Theta is close to zero or 180. */
            if ( fabs ( sinTheta ) <= tol )
            {
                auto Integer c, i ;

                /* . Get the diagonal elements of B = ( R + I ) / 2 which are ni*ni. */
                for ( i = 0 ; i < 3 ; i++ ) Vector3_Item ( axis, i ) = ( Matrix33_Item ( self, i, i ) + 1.0e+00 ) / 2.0e+00 ;

                /* . Set the angle using the trace of B (0 => tr(B) = tr(I) = 3, pi => tr(B) = tr(axis.dyadic.axis) = 1. */
                if ( Vector3_Sum ( axis ) > 2.0e+00 ) theta = 0.0e+00 ;
                else                                  theta = M_PI    ;

                /* . Find the index with the maximum magnitude and then use this column as the axis. */
                c = Vector3_AbsoluteMaximumIndex ( axis ) ;
                for ( i = 0 ; i < 3 ; i++ ) Vector3_Item ( axis, i ) = Matrix33_Item ( self, i, c ) ;
                Vector3_Item ( axis, c ) += 1.0e+00 ;
                Vector3_Normalize ( axis, &tol, NULL ) ;
            }
            /* . Other theta. */
            else
            {
                theta = atan2 ( sinTheta, cosTheta ) ;
                Vector3_Scale ( axis, 1.0e+00 / sinTheta ) ;
            }

            /* . Set the angle. */
            if ( angle != NULL ) (*angle) = theta ;
        }
        else Status_Set ( status, Status_InvalidArgument ) ;
    }
}
# undef Tolerance

/*----------------------------------------------------------------------------------------------------------------------------------
! . Determinant.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Matrix33_Determinant ( const Matrix33 *self )
{
    Real det = 0.0e+00 ;
    if ( self != NULL )
    {
        det = Matrix33_Item ( self, 0, 0 ) * ( Matrix33_Item ( self, 1, 1 ) * Matrix33_Item ( self, 2, 2 ) -
                                               Matrix33_Item ( self, 2, 1 ) * Matrix33_Item ( self, 1, 2 ) ) -
              Matrix33_Item ( self, 0, 1 ) * ( Matrix33_Item ( self, 1, 0 ) * Matrix33_Item ( self, 2, 2 ) -
                                               Matrix33_Item ( self, 2, 0 ) * Matrix33_Item ( self, 1, 2 ) ) +
              Matrix33_Item ( self, 0, 2 ) * ( Matrix33_Item ( self, 1, 0 ) * Matrix33_Item ( self, 2, 1 ) -
                                               Matrix33_Item ( self, 2, 0 ) * Matrix33_Item ( self, 1, 1 ) ) ;
    }
    return det ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Get the derivative of a matrix's inverse with respect to its ijth component.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define ITEM Matrix33_Item

void Matrix33_InverseDerivative ( const Matrix33 *self, const Integer i, const Integer j, Matrix33 *other )
{
    Matrix33_Set ( other, 0.0e+00 ) ;
    if ( ( self != NULL ) && ( other != NULL ) && ( i >= 0 ) && ( i < 3 ) && ( j >= 0 ) && ( j < 3 ) )
    {
        auto Real ddet = 0.0e+00, det, m00, m01, m02, m10, m11, m12, m20, m21, m22 ;

        /* . Local variables. */
        m00  = Matrix33_Item ( self, 0, 0 ) ;
        m01  = Matrix33_Item ( self, 0, 1 ) ;
        m02  = Matrix33_Item ( self, 0, 2 ) ;
        m10  = Matrix33_Item ( self, 1, 0 ) ;
        m11  = Matrix33_Item ( self, 1, 1 ) ;
        m12  = Matrix33_Item ( self, 1, 2 ) ;
        m20  = Matrix33_Item ( self, 2, 0 ) ;
        m21  = Matrix33_Item ( self, 2, 1 ) ;
        m22  = Matrix33_Item ( self, 2, 2 ) ;

        /* . Get the inverse. */
        Matrix33_Invert ( other, self ) ;

        /* . Determinant derivative factor. */
        switch ( 3 * i + j )
        {
            case 0: ddet = - m12 * m21 + m11 * m22 ; break ;
            case 1: ddet =   m12 * m20 - m10 * m22 ; break ;
            case 2: ddet = - m11 * m20 + m10 * m21 ; break ;
            case 3: ddet =   m02 * m21 - m01 * m22 ; break ;
            case 4: ddet = - m02 * m20 + m00 * m22 ; break ;
            case 5: ddet =   m01 * m20 - m00 * m21 ; break ;
            case 6: ddet = - m02 * m11 + m01 * m12 ; break ;
            case 7: ddet =   m02 * m10 - m00 * m12 ; break ;
            case 8: ddet = - m01 * m10 + m00 * m11 ; break ;
        }
        Matrix33_Scale ( other, - ddet ) ;

        /* . Matrix derivative factor. */
        switch ( 3 * i + j )
        {
            case 0: ITEM ( other, 1, 1 ) += m22 ; ITEM ( other, 1, 2 ) -= m12 ; ITEM ( other, 2, 1 ) -= m21 ; ITEM ( other, 2, 2 ) += m11 ; break ;
            case 1: ITEM ( other, 0, 1 ) -= m22 ; ITEM ( other, 0, 2 ) += m12 ; ITEM ( other, 2, 1 ) += m20 ; ITEM ( other, 2, 2 ) -= m10 ; break ;
            case 2: ITEM ( other, 0, 1 ) += m21 ; ITEM ( other, 0, 2 ) -= m11 ; ITEM ( other, 1, 1 ) -= m20 ; ITEM ( other, 1, 2 ) += m10 ; break ;
            case 3: ITEM ( other, 1, 0 ) -= m22 ; ITEM ( other, 1, 2 ) += m02 ; ITEM ( other, 2, 0 ) += m21 ; ITEM ( other, 2, 2 ) -= m01 ; break ;
            case 4: ITEM ( other, 0, 0 ) += m22 ; ITEM ( other, 0, 2 ) -= m02 ; ITEM ( other, 2, 0 ) -= m20 ; ITEM ( other, 2, 2 ) += m00 ; break ;
            case 5: ITEM ( other, 0, 0 ) -= m21 ; ITEM ( other, 0, 2 ) += m01 ; ITEM ( other, 1, 0 ) += m20 ; ITEM ( other, 1, 2 ) -= m00 ; break ;
            case 6: ITEM ( other, 1, 0 ) += m12 ; ITEM ( other, 1, 1 ) -= m02 ; ITEM ( other, 2, 0 ) -= m11 ; ITEM ( other, 2, 1 ) += m01 ; break ;
            case 7: ITEM ( other, 0, 0 ) -= m12 ; ITEM ( other, 0, 1 ) += m02 ; ITEM ( other, 2, 0 ) += m10 ; ITEM ( other, 2, 1 ) -= m00 ; break ;
            case 8: ITEM ( other, 0, 0 ) += m11 ; ITEM ( other, 0, 1 ) -= m01 ; ITEM ( other, 1, 0 ) -= m10 ; ITEM ( other, 1, 1 ) += m00 ; break ;
        }

        /* . Determinant factor. */
        det = Matrix33_Determinant ( self ) ;
        Matrix33_Scale ( other, 1.0e+00 / det ) ;
    }
}
# undef ADD

/*----------------------------------------------------------------------------------------------------------------------------------
! . Inversion.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Matrix33_Invert ( Matrix33 *self, const Matrix33 *other )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Real a00, a01, a02, a10, a11, a12, a20, a21, a22, det ;
        /* . Elements. */
        a00 = Matrix33_Item ( other, 0, 0 ) ;
        a01 = Matrix33_Item ( other, 0, 1 ) ;
        a02 = Matrix33_Item ( other, 0, 2 ) ;
        a10 = Matrix33_Item ( other, 1, 0 ) ;
        a11 = Matrix33_Item ( other, 1, 1 ) ;
        a12 = Matrix33_Item ( other, 1, 2 ) ;
        a20 = Matrix33_Item ( other, 2, 0 ) ;
        a21 = Matrix33_Item ( other, 2, 1 ) ;
        a22 = Matrix33_Item ( other, 2, 2 ) ;
        /* . Co-factors. */
        Matrix33_Item ( self, 0, 0 ) = ( a11 * a22 - a12 * a21 ) ;
        Matrix33_Item ( self, 0, 1 ) = ( a02 * a21 - a22 * a01 ) ;
        Matrix33_Item ( self, 0, 2 ) = ( a01 * a12 - a11 * a02 ) ;
        Matrix33_Item ( self, 1, 0 ) = ( a12 * a20 - a10 * a22 ) ;
        Matrix33_Item ( self, 1, 1 ) = ( a00 * a22 - a02 * a20 ) ;
        Matrix33_Item ( self, 1, 2 ) = ( a02 * a10 - a00 * a12 ) ;
        Matrix33_Item ( self, 2, 0 ) = ( a10 * a21 - a11 * a20 ) ;
        Matrix33_Item ( self, 2, 1 ) = ( a01 * a20 - a00 * a21 ) ;
        Matrix33_Item ( self, 2, 2 ) = ( a00 * a11 - a10 * a01 ) ;
        /* . Determinant factor. */
        det = Matrix33_Determinant ( other ) ;
        Matrix33_Scale ( self, 1.0e+00 / det ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Test for equality of contents.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define TOLERANCE 1.0e-6

Boolean Matrix33_IsEqual ( const Matrix33 *self, const Matrix33 *other )
{
    Boolean QEQUAL = False ;
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Integer i, j ;
        QEQUAL = True ;
        for ( i = 0 ; i < 3 ; i++ )
        {
            for ( j = 0 ; j < 3 ; j++ )
            {
                if ( fabs ( Matrix33_Item ( self, i, j ) - Matrix33_Item ( other, i, j ) ) > TOLERANCE )
                {
                    QEQUAL = False ;
                    break ;
                }
            }
            if ( ! QEQUAL ) break ;
        }
    }
    return QEQUAL ;
}

# undef TOLERANCE

/*----------------------------------------------------------------------------------------------------------------------------------
! . Test for identity.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define TOLERANCE 1.0e-6

Boolean Matrix33_IsIdentity ( const Matrix33 *self )
{
    Boolean QIDENTITY = False ;
    if ( self != NULL )
    {
        auto Real target ;
        auto Integer    i, j   ;
        QIDENTITY = True ;
        for ( i = 0 ; i < 3 ; i++ )
        {
            for ( j = 0 ; j < 3 ; j++ )
            {
                if ( i == j ) target = 1.0e+00 ;
                else          target = 0.0e+00 ;
                if ( fabs ( Matrix33_Item ( self, i, j ) - target ) > TOLERANCE )
                {
                    QIDENTITY = False ;
                    break ;
                }
            }
            if ( ! QIDENTITY ) break ;
        }
    }
    return QIDENTITY ;
}

# undef TOLERANCE

/*----------------------------------------------------------------------------------------------------------------------------------
! . Test for improper and proper rotations.
! . First check the determinant and then for orthogonality.
!---------------------------------------------------------------------------------------------------------------------------------*/
# define TOLERANCE  1.0e-06

# define TARGETDET -1.0e+00
Boolean Matrix33_IsImproperRotation ( const Matrix33 *self )
{
    Boolean QTEST = False ;
    if ( self != NULL )
    {
        auto Real det ;
        det = Matrix33_Determinant ( self ) ;
        if ( fabs ( det - TARGETDET ) < TOLERANCE )
        {
            QTEST = Matrix33_IsOrthogonal ( self ) ;
        }
    }
    return QTEST ;
}
# undef TARGETDET

# define TARGETDET 1.0e+00
Boolean Matrix33_IsProperRotation ( const Matrix33 *self )
{
    Boolean QTEST = False ;
    if ( self != NULL )
    {
        auto Real det ;
        det = Matrix33_Determinant ( self ) ;
        if ( fabs ( det - TARGETDET ) < TOLERANCE )
        {
            QTEST = Matrix33_IsOrthogonal ( self ) ;
        }
    }
    return QTEST ;
}
# undef TARGETDET

# undef TOLERANCE

/*----------------------------------------------------------------------------------------------------------------------------------
! . Test for orthogonality.
!---------------------------------------------------------------------------------------------------------------------------------*/
Boolean Matrix33_IsOrthogonal ( const Matrix33 *self )
{
    Boolean QTEST = False ;
    if ( self != NULL )
    {
        auto Matrix33 *other ;
        other = Matrix33_Allocate ( ) ;
        Matrix33_CopyTo         ( self, other, NULL ) ;
        Matrix33_Transpose      ( other, NULL ) ;
        Matrix33_PostMultiplyBy ( other, self ) ;
        QTEST = Matrix33_IsIdentity ( other ) ;
        Matrix33_Deallocate ( &other ) ;
    }
    return QTEST ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Post-multiply by.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Matrix33_PostMultiplyBy ( Matrix33 *self, const Matrix33 *other )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Real o0, o1, o2, r0, r1, r2, s00, s01, s02, s10, s11, s12, s20, s21, s22 ;
        auto Integer    i ;
        s00 = Matrix33_Item ( self, 0, 0 ) ;
        s01 = Matrix33_Item ( self, 0, 1 ) ;
        s02 = Matrix33_Item ( self, 0, 2 ) ;
        s10 = Matrix33_Item ( self, 1, 0 ) ;
        s11 = Matrix33_Item ( self, 1, 1 ) ;
        s12 = Matrix33_Item ( self, 1, 2 ) ;
        s20 = Matrix33_Item ( self, 2, 0 ) ;
        s21 = Matrix33_Item ( self, 2, 1 ) ;
        s22 = Matrix33_Item ( self, 2, 2 ) ;
        /* . Multiply. */
        for ( i = 0 ; i < other->extent1 ; i++ )
        {
            o0 = Matrix33_Item ( other, 0, i ) ;
            o1 = Matrix33_Item ( other, 1, i ) ;
            o2 = Matrix33_Item ( other, 2, i ) ;
            r0 = s00 * o0 + s01 * o1 + s02 * o2 ;
            r1 = s10 * o0 + s11 * o1 + s12 * o2 ;
            r2 = s20 * o0 + s21 * o1 + s22 * o2 ;
            Matrix33_Item ( self, 0, i ) = r0 ;
            Matrix33_Item ( self, 1, i ) = r1 ;
            Matrix33_Item ( self, 2, i ) = r2 ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Pre-multiply by.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Matrix33_PreMultiplyBy ( Matrix33 *self, const Matrix33 *other )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Real o0, o1, o2, r0, r1, r2, s00, s01, s02, s10, s11, s12, s20, s21, s22 ;
        auto Integer    i ;
        s00 = Matrix33_Item ( self, 0, 0 ) ;
        s01 = Matrix33_Item ( self, 0, 1 ) ;
        s02 = Matrix33_Item ( self, 0, 2 ) ;
        s10 = Matrix33_Item ( self, 1, 0 ) ;
        s11 = Matrix33_Item ( self, 1, 1 ) ;
        s12 = Matrix33_Item ( self, 1, 2 ) ;
        s20 = Matrix33_Item ( self, 2, 0 ) ;
        s21 = Matrix33_Item ( self, 2, 1 ) ;
        s22 = Matrix33_Item ( self, 2, 2 ) ;
        /* . Multiply. */
        for ( i = 0 ; i < other->extent0 ; i++ )
        {
            Matrix33_GetRow ( other, i, o0, o1, o2 ) ;
            r0 = s00 * o0 + s10 * o1 + s20 * o2 ;
            r1 = s01 * o0 + s11 * o1 + s21 * o2 ;
            r2 = s02 * o0 + s12 * o1 + s22 * o2 ;
            Matrix33_SetRow ( self, i, r0, r1, r2 ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Determine a matrix33 for a reflection in a plane through the origin.
! . Normal is the normal to the plane and is assumed normalized.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status Matrix33_Reflection ( Matrix33 **self, const Vector3 *normal )
{
    Status status = Status_OK ;
    if ( self != NULL )
    {
        auto Matrix33 *work = NULL ;
        /* . Allocate space if necessary. */
        if ( (*self) == NULL ) { work = Matrix33_Allocate ( ) ; (*self) = work ; }
        else                     work = (*self) ;
        /* . Calculate the rotation. */
        if ( work != NULL )
        {
            auto Real a, b, c ;
            Matrix33_Set ( work, 0.0e+00 ) ;
            a = Vector3_Item ( normal, 0 ) ;
            b = Vector3_Item ( normal, 1 ) ;
            c = Vector3_Item ( normal, 2 ) ;
            Matrix33_Item ( work, 0, 0 ) = 1.0e+00 - 2.0e+00 * a * a ;
            Matrix33_Item ( work, 0, 1 ) =         - 2.0e+00 * a * b ;
            Matrix33_Item ( work, 0, 2 ) =         - 2.0e+00 * a * c ;
            Matrix33_Item ( work, 1, 0 ) =         - 2.0e+00 * b * a ;
            Matrix33_Item ( work, 1, 1 ) = 1.0e+00 - 2.0e+00 * b * b ;
            Matrix33_Item ( work, 1, 2 ) =         - 2.0e+00 * b * c ;
            Matrix33_Item ( work, 2, 0 ) =         - 2.0e+00 * c * a ;
            Matrix33_Item ( work, 2, 1 ) =         - 2.0e+00 * c * b ;
            Matrix33_Item ( work, 2, 2 ) = 1.0e+00 - 2.0e+00 * c * c ;
        }
        else status = Status_OutOfMemory ;
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Determine a matrix33 for a rotation about an axis (assumed normalized).
! . The rotation is anticlockwise when looking in the direction the axis is
! . travelling and clockwise when looking towards where the axis has come.
! . Thus (1,0,0) -> (0,-1,0) for a pi/2 rotation about (0,0,1).
!---------------------------------------------------------------------------------------------------------------------------------*/
Status Matrix33_RotationAboutAxis ( Matrix33 **self, const Real angle, const Real x, const Real y, const Real z )
{
    Status status = Status_OK ;
    if ( self != NULL )
    {
        auto Matrix33 *work = NULL ;
        /* . Allocate space if necessary. */
        if ( (*self) == NULL ) { work = Matrix33_Allocate ( ) ; (*self) = work ; }
        else                     work = (*self) ;
        /* . Calculate the rotation. */
        if ( work != NULL )
        {
            auto Real q0, q1, q2, q3 ;
            /* . Determine some factors. */
            q0 =     cos ( angle / 2.0e+00 ) ;
            q1 = x * sin ( angle / 2.0e+00 ) ;
            q2 = y * sin ( angle / 2.0e+00 ) ;
            q3 = z * sin ( angle / 2.0e+00 ) ;
            /* . Create the matrix. */
            status = Matrix33_RotationFromQuaternion ( self, q0, q1, q2, q3 ) ;
        }
        else status = Status_OutOfMemory ;
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Determine a matrix33 of a rotation from a quaternion (assumed normalized).
!---------------------------------------------------------------------------------------------------------------------------------*/
Status Matrix33_RotationFromQuaternion ( Matrix33 **self, const Real q0, const Real q1, const Real q2, const Real q3 )
{
    Status status = Status_OK ;
    if ( self != NULL )
    {
        auto Matrix33 *work = NULL ;
        /* . Allocate space if necessary. */
        if ( (*self) == NULL ) { work = Matrix33_Allocate ( ) ; (*self) = work ; }
        else                     work = (*self) ;
        /* . Calculate the rotation. */
        if ( work != NULL )
        {
            Matrix33_Item ( work, 0, 0 ) = q0*q0 + q1*q1 - q2*q2 - q3*q3 ;
            Matrix33_Item ( work, 1, 0 ) = 2.0e+00 * ( - q0*q3 + q1*q2 ) ;
            Matrix33_Item ( work, 2, 0 ) = 2.0e+00 * (   q0*q2 + q1*q3 ) ;
            Matrix33_Item ( work, 0, 1 ) = 2.0e+00 * (   q0*q3 + q1*q2 ) ;
            Matrix33_Item ( work, 1, 1 ) = q0*q0 - q1*q1 + q2*q2 - q3*q3 ;
            Matrix33_Item ( work, 2, 1 ) = 2.0e+00 * ( - q0*q1 + q2*q3 ) ;
            Matrix33_Item ( work, 0, 2 ) = 2.0e+00 * ( - q0*q2 + q1*q3 ) ;
            Matrix33_Item ( work, 1, 2 ) = 2.0e+00 * (   q0*q1 + q2*q3 ) ;
            Matrix33_Item ( work, 2, 2 ) = q0*q0 - q1*q1 - q2*q2 + q3*q3 ;
        }
        else status = Status_OutOfMemory ;
    }
    return status ;
}
