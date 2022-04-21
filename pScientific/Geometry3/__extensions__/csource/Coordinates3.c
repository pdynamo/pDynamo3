/*==================================================================================================================================
! . Procedures for 3-D arrays of coordinates.
!=================================================================================================================================*/

# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "BooleanBlock.h"
# include "Coordinates3.h"
# include "DenseEigenvalueSolvers.h"
# include "IntegerArray1D.h"
# include "NumericalMacros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . The value to return if there are problems with a calculation. */
# define BadValue 1.0e+30

/* . Debugging. */
# define DEBUG
/* # define DEBUGPRINTING */

/* . Local unit conversion until this module is sorted out. */
/* . Angle. */
# define UNITS_ANGLE_DEGREES_TO_RADIANS ( M_PI / 180.0e+00 )
# define UNITS_ANGLE_RADIANS_TO_DEGREES ( 180.0e+00 / M_PI )

/*----------------------------------------------------------------------------------------------------------------------------------
! . Allocation.
!---------------------------------------------------------------------------------------------------------------------------------*/
Coordinates3 *Coordinates3_Allocate ( const Integer extent, Status *status ) { return RealArray2D_AllocateWithExtents ( extent, 3, status ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate an angle between three points.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Coordinates3_Angle ( const Coordinates3 *self, const Integer i, const Integer j, const Integer k )
{
    Real costheta, rij, rkj, sintheta2, theta = 0.0e+00, xij, yij, zij, xkj, ykj, zkj ;
    if ( self != NULL )
    {
        /* . Displacement vectors. */
        Coordinates3_DifferenceRow ( self, i, j, xij, yij, zij ) ;
        rij = sqrt ( xij * xij + yij * yij + zij * zij ) ;
        xij /= rij ; yij /= rij ; zij /= rij ;
        Coordinates3_DifferenceRow ( self, k, j, xkj, ykj, zkj ) ;
        rkj = sqrt ( xkj * xkj + ykj * ykj + zkj * zkj ) ;
        xkj /= rkj ; ykj /= rkj ; zkj /= rkj ;
        /* . Angle terms. */
        costheta  = xij * xkj + yij * ykj + zij * zkj ;
        sintheta2 = pow ( yij * zkj - zij * ykj, 2 ) + pow ( zij * xkj - xij * zkj, 2 ) + pow ( xij * ykj - yij * xkj, 2 ) ;
        costheta /= sqrt ( costheta * costheta + sintheta2 ) ;
        theta     = UNITS_ANGLE_RADIANS_TO_DEGREES * acos ( costheta ) ;
    }
    return theta ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Build a point |i| as point |j| + |r| * |direction|.
! . |direction| should be normalized.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status Coordinates3_BuildPointFromDistance ( Coordinates3 *self, const Integer i, const Integer j, const Real r, const Vector3 *direction )
{
    Status status = Status_OK ;
    if ( ( self != NULL ) && ( direction != NULL ) )
    {
        if ( ( i >= 0 ) && ( i < self->extent0 ) && ( j >= 0 ) && ( j < self->extent0 ) )
        {
            auto Integer c ;
            for ( c = 0 ; c < 3 ; c++ )
            {
                Coordinates3_Item ( self, i, c ) = Coordinates3_Item ( self, j, c ) + r * ( Vector3_Item ( direction, c ) ) ;
            }
        }
        else status = Status_IndexOutOfRange ;
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Build a point |i| such that it is at a distance |r| from |j| and the angle
! . i-j-k is |theta|.
! . |direction| does not have to be normalized.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status Coordinates3_BuildPointFromDistanceAngle ( Coordinates3 *self, const Integer i, const Integer j, const Integer k, const Real r, const Real theta, const Vector3 *direction )
{
    Status status = Status_OK ;
    if ( ( self != NULL ) && ( direction != NULL ) )
    {
        if ( ( i >= 0 ) && ( i < self->extent0 ) && ( j >= 0 ) && ( j < self->extent0 ) && ( k >= 0 ) && ( k < self->extent0 ) )
        {
            auto Vector3 *dra = NULL, *drkj = NULL ;
            dra  = Vector3_Allocate ( ) ;
            drkj = Vector3_Allocate ( ) ;
            if ( ( dra != NULL ) && ( drkj != NULL ) )
            {
                auto Real dx, dy, dz ;
                /* . Get the normalized j->k displacement. */
                Coordinates3_DifferenceRow ( self, k, j, dx, dy, dz ) ;
                Vector3_Item ( drkj, 0 ) = dx ;
                Vector3_Item ( drkj, 1 ) = dy ;
                Vector3_Item ( drkj, 2 ) = dz ;
                Vector3_Normalize ( drkj, NULL, &status ) ;
                if ( status == Status_OK )
                {
                    /* . Get the normalized vector perpendicular to drkj and direction. */
                    Vector3_CopyTo ( direction, dra, NULL ) ;
                    Vector3_CrossProduct ( dra, drkj ) ;
                    Vector3_Normalize ( dra, NULL, &status ) ;
                    if ( status == Status_OK )
                    {
                        auto Real wa, wb ;
                        auto Integer c ;
                        /* . Calculate the cross product of drkj and dra (keep in dra). */
                        Vector3_CrossProduct ( dra, drkj ) ;
                        Vector3_Scale        ( dra, -1.0 ) ;
                        /* . Calculate the coordinate displacements (spherical polars). */
                        wa = r * cos ( UNITS_ANGLE_DEGREES_TO_RADIANS * theta ) ;
                        wb = r * sin ( UNITS_ANGLE_DEGREES_TO_RADIANS * theta ) ;
                        /* . Calculate the coordinates of i. */
                        for ( c = 0 ; c < 3 ; c++ )
                        {
                            Coordinates3_Item ( self, i, c ) = Coordinates3_Item ( self, j, c ) + wa * ( Vector3_Item ( drkj, c ) ) + wb * ( Vector3_Item ( dra, c ) ) ;
                        }
                    }
                }
            }
            else status = Status_OutOfMemory ;
            /* . Finish up. */
            Vector3_Deallocate ( &dra  ) ;
            Vector3_Deallocate ( &drkj ) ;
        }
        else status = Status_IndexOutOfRange ;
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Build a point |i| such that it is at a distance |r| from |j|, the angle
! . i-j-k is |theta| and the dihedral i-j-k-l is |phi|.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status Coordinates3_BuildPointFromDistanceAngleDihedral ( Coordinates3 *self, const Integer i, const Integer j, const Integer k, const Integer l, const Real r, const Real theta, const Real phi )
{
    Status status = Status_OK ;
    if ( self != NULL )
    {
        if ( ( i >= 0 ) && ( i < self->extent0 ) && ( j >= 0 ) && ( j < self->extent0 ) && ( k >= 0 ) && ( k < self->extent0 ) && ( l >= 0 ) && ( l < self->extent0 ) )
        {
            auto Vector3 *dra = NULL, *drkj = NULL, *drlk = NULL ;
            dra  = Vector3_Allocate ( ) ;
            drkj = Vector3_Allocate ( ) ;
            drlk = Vector3_Allocate ( ) ;
            if ( ( dra != NULL ) && ( drkj != NULL ) && ( drlk != NULL ) )
            {
                auto Real dx, dy, dz ;
                /* . Get the normalized j->k displacement. */
                Coordinates3_DifferenceRow ( self, k, j, dx, dy, dz ) ;
                Vector3_Item ( drkj, 0 ) = dx ;
                Vector3_Item ( drkj, 1 ) = dy ;
                Vector3_Item ( drkj, 2 ) = dz ;
                Vector3_Normalize ( drkj, NULL, &status ) ;
                if ( status == Status_OK )
                {
                    /* . Get the k->l displacement. */
                    Coordinates3_DifferenceRow ( self, l, k, dx, dy, dz ) ;
                    Vector3_Item ( drlk, 0 ) = dx ;
                    Vector3_Item ( drlk, 1 ) = dy ;
                    Vector3_Item ( drlk, 2 ) = dz ;
                    /* . Get the normalized vector perpendicular to drlk and drkj. */
                    Vector3_CopyTo ( drlk, dra, NULL ) ;
                    Vector3_CrossProduct ( dra, drkj ) ;
                    Vector3_Normalize ( dra, NULL, &status ) ;
                    if ( status == Status_OK )
                    {
                        auto Real sint, wa, wb, wc ;
                        auto Integer c ;
                        /* . Calculate the cross product of drkj and dra (put in drlk). */
                        Vector3_CopyTo ( drkj, drlk, NULL ) ;
                        Vector3_CrossProduct ( drlk, dra ) ;
                        /* . Calculate the coordinate displacements (spherical polars). */
                        sint =            sin ( UNITS_ANGLE_DEGREES_TO_RADIANS * theta ) ;
                        wa   = r *        cos ( UNITS_ANGLE_DEGREES_TO_RADIANS * theta ) ;
                        wb   = r * sint * cos ( UNITS_ANGLE_DEGREES_TO_RADIANS * phi   ) ;
                        wc   = r * sint * sin ( UNITS_ANGLE_DEGREES_TO_RADIANS * phi   ) ;
                        /* . Calculate the coordinates of i. */
                        for ( c = 0 ; c < 3 ; c++ )
                        {
                            Coordinates3_Item ( self, i, c ) = Coordinates3_Item ( self, j, c ) + wa * ( Vector3_Item ( drkj, c ) ) + wb * ( Vector3_Item ( drlk, c ) ) + wc * ( Vector3_Item ( dra, c ) ) ;
                        }
                    }
                }
            }
            else status = Status_OutOfMemory ;
            /* . Finish up. */
            Vector3_Deallocate ( &dra  ) ;
            Vector3_Deallocate ( &drkj ) ;
            Vector3_Deallocate ( &drlk ) ;
        }
        else status = Status_IndexOutOfRange ;
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Build a point |i| such that it is at a distance |r| from |j|. The direction
! . is defined by the bisector of j-i-k and the vector perpendicular to this
! . and the plane j-i-k. The |planeangle| is the angle from the bisector.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status Coordinates3_BuildPointFromDistancePlaneAngle ( Coordinates3 *self, const Integer i, const Integer j, const Integer k, const Integer l, const Real r, const Real planeangle )
{
    Status status = Status_OK ;
    if ( self != NULL )
    {
        if ( ( i >= 0 ) && ( i < self->extent0 ) && ( j >= 0 ) && ( j < self->extent0 ) && ( k >= 0 ) && ( k < self->extent0 ) && ( l >= 0 ) && ( l < self->extent0 ) )
        {
            auto Vector3 *dra = NULL, *drb = NULL ;
            dra = Vector3_Allocate ( ) ;
            drb = Vector3_Allocate ( ) ;
            if ( ( dra != NULL ) && ( drb != NULL ) )
            {
                auto Real dx, dy, dz ;
                /* . Get the normalized j->k displacement. */
                Coordinates3_DifferenceRow ( self, k, j, dx, dy, dz ) ;
                Vector3_Item ( dra, 0 ) = dx ;
                Vector3_Item ( dra, 1 ) = dy ;
                Vector3_Item ( dra, 2 ) = dz ;
                Vector3_Normalize ( dra, NULL, &status ) ;
                if ( status == Status_OK )
                {
                    /* . Get the normalized j->l displacement. */
                    Coordinates3_DifferenceRow ( self, l, j, dx, dy, dz ) ;
                    Vector3_Item ( drb, 0 ) = dx ;
                    Vector3_Item ( drb, 1 ) = dy ;
                    Vector3_Item ( drb, 2 ) = dz ;
                    Vector3_Normalize ( drb, NULL, &status ) ;
                    if ( status == Status_OK )
                    {
                        auto Real a, b ;
                        auto Integer    c    ;
                        /* . Get the bisector of j->k and j->l and the vector perpendicular to it. */
                        for ( c = 0 ; c < 3 ; c++ )
                        {
                            a = Vector3_Item ( dra, c ) ;
                            b = Vector3_Item ( drb, c ) ;
                            Vector3_Item ( dra, c ) = a + b ;
                            Vector3_Item ( drb, c ) = b - a ;
                        }
                        Vector3_Normalize ( dra, NULL, &status ) ;
                        if ( status == Status_OK )
                        {
                            Vector3_Normalize ( drb, NULL, &status ) ;
                            if ( status == Status_OK )
                            {
                                /* . Get the normalized vector perpendicular to dra and drb. */
                                Vector3_CrossProduct ( drb, dra ) ;
                                Vector3_Normalize ( drb, NULL, &status ) ;
                                if ( status == Status_OK )
                                {
                                    /* . Calculate the angle factors. */
                                    a = r * cos ( UNITS_ANGLE_DEGREES_TO_RADIANS * planeangle ) ;
                                    b = r * sin ( UNITS_ANGLE_DEGREES_TO_RADIANS * planeangle ) ;
                                    /* . Calculate the coordinates of i. */
                                    for ( c = 0 ; c < 3 ; c++ )
                                    {
                                        Coordinates3_Item ( self, i, c ) = Coordinates3_Item ( self, j, c ) + a * ( Vector3_Item ( dra, c ) ) + b * ( Vector3_Item ( drb, c ) ) ;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else status = Status_OutOfMemory ;
            /* . Finish up. */
            Vector3_Deallocate ( &dra ) ;
            Vector3_Deallocate ( &drb ) ;
        }
        else status = Status_IndexOutOfRange ;
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Build a point |i| such that it is at a distance |r| from |j|. The direction
! . is defined by the vector from the mid-point of the k-l-m plane to |j|.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status Coordinates3_BuildPointFromDistanceTetrahedralTripod ( Coordinates3 *self, const Integer i, const Integer j, const Integer k, const Integer l, const Integer m, const Real r )
{
    Status status = Status_OK ;
    if ( self != NULL )
    {
        if ( ( i >= 0 ) && ( i < self->extent0 ) && ( j >= 0 ) && ( j < self->extent0 ) && ( k >= 0 ) && ( k < self->extent0 ) && ( l >= 0 ) && ( l < self->extent0 ) && ( m >= 0 ) && ( m < self->extent0 ) )
        {
            auto Vector3 *dra = NULL, *drt = NULL ;
            dra = Vector3_Allocate ( ) ;
            drt = Vector3_Allocate ( ) ;
            if ( ( dra != NULL ) && ( drt != NULL ) )
            {
                auto Real dx, dy, dz ;
                /* . Get the normalized j->k displacement. */
                Coordinates3_DifferenceRow ( self, k, j, dx, dy, dz ) ;
                Vector3_Item ( dra, 0 ) = dx ;
                Vector3_Item ( dra, 1 ) = dy ;
                Vector3_Item ( dra, 2 ) = dz ;
                Vector3_Normalize ( dra, NULL, &status ) ;
                if ( status == Status_OK )
                {
                    /* . Get the normalized j->l displacement. */
                    Coordinates3_DifferenceRow ( self, l, j, dx, dy, dz ) ;
                    Vector3_Item ( drt, 0 ) = dx ;
                    Vector3_Item ( drt, 1 ) = dy ;
                    Vector3_Item ( drt, 2 ) = dz ;
                    Vector3_Normalize ( drt, NULL, &status ) ;
                    if ( status == Status_OK )
                    {
                        /* . Increment dra. */
                        Vector3_Add ( dra, 1.0e+00, drt, NULL ) ;
                        /* . Get the normalized j->m displacement. */
                        Coordinates3_DifferenceRow ( self, m, j, dx, dy, dz ) ;
                        Vector3_Item ( drt, 0 ) = dx ;
                        Vector3_Item ( drt, 1 ) = dy ;
                        Vector3_Item ( drt, 2 ) = dz ;
                        Vector3_Normalize ( drt, NULL, &status ) ;
                        if ( status == Status_OK )
                        {
                            /* . Increment dra and normalize it. */
                            Vector3_Add ( dra, 1.0e+00, drt, NULL ) ;
                            Vector3_Normalize ( drt, NULL, &status ) ;
                            if ( status == Status_OK )
                            {
                                auto Integer c ;
                                /* . Calculate the coordinates of i. */
                                for ( c = 0 ; c < 3 ; c++ )
                                {
                                    Coordinates3_Item ( self, i, c ) = Coordinates3_Item ( self, j, c ) - r * ( Vector3_Item ( dra, c ) ) ;
                                }
                            }
                        }
                    }
                }
            }
            else status = Status_OutOfMemory ;
            /* . Finish up. */
            Vector3_Deallocate ( &dra ) ;
            Vector3_Deallocate ( &drt ) ;
        }
        else status = Status_IndexOutOfRange ;
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate a center.
!---------------------------------------------------------------------------------------------------------------------------------*/
Status Coordinates3_Center ( const Coordinates3 *self, const Selection *selection, const RealArray1D *weights, Vector3 **center )
{
    Status status = Status_OK ;
    if ( self != NULL )
    {
        auto Vector3 *work = NULL ;
        if ( (*center) == NULL ) { work = Vector3_Allocate ( ) ; (*center) = work ; }
        else                       work = (*center) ;
        /* . Calculate the center. */
        if ( work == NULL ) status = Status_OutOfMemory ;
        else Coordinates3_CenterRaw ( self, selection, weights, Vector3_Data ( work ), &status ) ;
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate a center to a raw real array.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* Could add stride here? */
void Coordinates3_CenterRaw ( const Coordinates3 *self, const Selection *selection, const RealArray1D *weights, Real *data, Status *status )
{
    if ( ( self != NULL ) && ( data != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer  i, n = self->extent0, s ;
        auto Real     w, wTotal = 0.0e+00, x, y, z ;
        for ( i = 0 ; i < 3 ; i++ ) data[i] = 0.0e+00 ;
        /* . Treat the four cases. */
        if ( selection == NULL )
        {
            if ( weights == NULL )
            {
                wTotal = ( Real ) n ;
                for ( i = 0 ; i < n ; i++ )
                {
                    Coordinates3_GetRow ( self, i, x, y, z ) ;
                    data[0] += x ;
                    data[1] += y ;
                    data[2] += z ;
                }
            }
            else
            {
                wTotal = RealArray1D_Sum ( weights ) ;
                for ( i = 0 ; i < n ; i++ )
                {
                    w = Array1D_Item ( weights, i ) ;
                    Coordinates3_GetRow ( self, i, x, y, z ) ;
                    data[0] += w * x ;
                    data[1] += w * y ;
                    data[2] += w * z ;
                }
            }
        }
        else
        {
            if ( weights == NULL )
            {
                wTotal = ( Real ) selection->capacity ;
                for ( s = 0 ; s < selection->capacity ; s++ )
                {
                    i = selection->indices[s] ;
                    Coordinates3_GetRow ( self, i, x, y, z ) ;
                    data[0] += x ;
                    data[1] += y ;
                    data[2] += z ;
                }
            }
            else
            {
                for ( s = 0 ; s < selection->capacity ; s++ )
                {
                    i = selection->indices[s] ;
                    w = Array1D_Item ( weights, i ) ;
                    Coordinates3_GetRow ( self, i, x, y, z ) ;
                    data[0] += w * x ;
                    data[1] += w * y ;
                    data[2] += w * z ;
                    wTotal  += w ;
                }
            }
        }
        /* . Scaling. */
        if ( wTotal != 0.0e+00 ) { for ( i = 0 ; i < 3 ; i++ ) data[i] /= wTotal ; }
        else Status_Set ( status, Status_AlgorithmError ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate a dihedral angle between four points.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Coordinates3_Dihedral ( const Coordinates3 *self, const Integer i, const Integer j, const Integer k, const Integer l )
{
    Real phi = 0.0e+00 ;
    if ( self != NULL )
    {
        auto Real cosphi, m2, mn, mx, my, mz, n2, norm2, nx, ny, nz, rkj, rkj2, sinphi, xij, yij, zij, xkj, ykj, zkj, xlk, ylk, zlk ;
        /* . Coordinate displacements. */
        Coordinates3_DifferenceRow ( self, i, j, xij, yij, zij ) ;
        Coordinates3_DifferenceRow ( self, k, j, xkj, ykj, zkj ) ;
        Coordinates3_DifferenceRow ( self, l, k, xlk, ylk, zlk ) ;
        rkj2 = xkj * xkj + ykj * ykj + zkj * zkj ;
        rkj  = sqrt ( rkj2 ) ;
        /* . m and n. */
        mx = yij * zkj - zij * ykj ;
        my = zij * xkj - xij * zkj ;
        mz = xij * ykj - yij * xkj ;
        nx = ylk * zkj - zlk * ykj ;
        ny = zlk * xkj - xlk * zkj ;
        nz = xlk * ykj - ylk * xkj ;
        m2 = mx * mx + my * my + mz * mz ;
        n2 = nx * nx + ny * ny + nz * nz ;
        mn = sqrt ( m2 * n2 ) ;
        /* . Cosine and sine of the dihedral. */
        cosphi  =       (  mx * nx +  my * ny +  mz * nz ) / mn ;
        sinphi  = rkj * ( xij * nx + yij * ny + zij * nz ) / mn ;
        /* . Normalize. */
        norm2   = cosphi * cosphi + sinphi * sinphi ;
        cosphi /= norm2 ;
        /* . Calculate the angle. */
        phi = UNITS_ANGLE_RADIANS_TO_DEGREES * acos ( cosphi ) * Sign ( sinphi ) ;
    }
    return phi ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate a distance between two points.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Coordinates3_Distance ( const Coordinates3 *self, const Integer i, const Integer j )
{
    Real rij = 0.0e+00, xij, yij, zij ;
    if ( self != NULL )
    {
        Coordinates3_DifferenceRow ( self, i, j, xij, yij, zij ) ;
        rij = sqrt ( xij * xij + yij * yij + zij * zij ) ;
    }
    return rij ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Finding the enclosing orthorhombic box for the matrix.
! . Origin is the bottom left hand corner of the box and extents are the extents of the box sides.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Coordinates3_EnclosingOrthorhombicBox ( const Coordinates3 *self, const Selection *selection, const RealArray1D *radii, Vector3 *origin, Vector3 *extents )
{
    if ( ( self != NULL ) && ( origin != NULL ) && ( extents != NULL ) )
    {
        if ( self->extent0 > 0 )
        {
            auto Real x, y, z, xmax, ymax, zmax, xmin, ymin, zmin ;
            auto Integer    i, s ;
            /* . The four cases. */
            if ( selection == NULL )
            {
                /* . Initialization. */
                Coordinates3_GetRow ( self, 0, x, y, z ) ;
                xmax = x ; xmin = x ;
                ymax = y ; ymin = y ;
                zmax = z ; zmin = z ;
                /* . Without radii. */
                if ( radii == NULL )
                {
                    for ( i = 1 ; i < self->extent0 ; i++ )
                    {
                        Coordinates3_GetRow ( self, i, x, y, z ) ;
                        if ( x > xmax ) xmax = x ;
                        if ( x < xmin ) xmin = x ;
                        if ( y > ymax ) ymax = y ;
                        if ( y < ymin ) ymin = y ;
                        if ( z > zmax ) zmax = z ;
                        if ( z < zmin ) zmin = z ;
                    }
                }
                /* . With radii. */
                else
                {
                    auto Real r ;
                    r = Array1D_Item ( radii, 0 ) ;
                    xmax += r ; xmin -= r ;
                    ymax += r ; ymin -= r ;
                    zmax += r ; zmin -= r ;
                    for ( i = 1 ; i < self->extent0 ; i++ )
                    {
                        r = Array1D_Item ( radii, i ) ;
                        Coordinates3_GetRow ( self, i, x, y, z ) ;
                        if ( x + r > xmax ) xmax = x + r ;
                        if ( x - r < xmin ) xmin = x - r ;
                        if ( y + r > ymax ) ymax = y + r ;
                        if ( y - r < ymin ) ymin = y - r ;
                        if ( z + r > zmax ) zmax = z + r ;
                        if ( z - r < zmin ) zmin = z - r ;
                    }
                }
            }
            else
            {
                /* . Initialization. */
                i = selection->indices[0] ;
                Coordinates3_GetRow ( self, i, x, y, z ) ;
                xmax = x ; xmin = x ;
                ymax = y ; ymin = y ;
                zmax = z ; zmin = z ;
                /* . Without radii. */
                if ( radii == NULL )
                {
                    for ( s = 1 ; s < selection->capacity ; s++ )
                    {
                        i = selection->indices[s] ;
                        Coordinates3_GetRow ( self, i, x, y, z ) ;
                        if ( x > xmax ) xmax = x ;
                        if ( x < xmin ) xmin = x ;
                        if ( y > ymax ) ymax = y ;
                        if ( y < ymin ) ymin = y ;
                        if ( z > zmax ) zmax = z ;
                        if ( z < zmin ) zmin = z ;
                    }
                }
                /* . With radii. */
                else
                {
                    auto Real r ;
                    r = Array1D_Item ( radii, 0 ) ;
                    xmax += r ; xmin -= r ;
                    ymax += r ; ymin -= r ;
                    zmax += r ; zmin -= r ;
                    for ( s = 1 ; s < selection->capacity ; s++ )
                    {
                        i = selection->indices[s] ;
                        r = Array1D_Item ( radii, i ) ;
                        Coordinates3_GetRow ( self, i, x, y, z ) ;
                        if ( x + r > xmax ) xmax = x + r ;
                        if ( x - r < xmin ) xmin = x - r ;
                        if ( y + r > ymax ) ymax = y + r ;
                        if ( y - r < ymin ) ymin = y - r ;
                        if ( z + r > zmax ) zmax = z + r ;
                        if ( z - r < zmin ) zmin = z - r ;
                    }
                }
            }
            /* . Set the values. */
            Vector3_Item ( origin , 0 ) = xmin ;
            Vector3_Item ( origin , 1 ) = ymin ;
            Vector3_Item ( origin , 2 ) = zmin ;
            Vector3_Item ( extents, 0 ) = xmax - xmin ;
            Vector3_Item ( extents, 1 ) = ymax - ymin ;
            Vector3_Item ( extents, 2 ) = zmax - zmin ;
        }
        else
        {
            Vector3_Set ( extents, 0.0e+00 ) ;
            Vector3_Set ( origin , 0.0e+00 ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Generate coordinates from a regular grid.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Coordinates3_FromRegularGrid ( Coordinates3 *self, const RegularGrid *grid, Selection *selection, Status *status )
{
    if ( ( self != NULL ) && ( grid != NULL ) && Status_IsOK ( status ) )
    {
        auto Integer n, nGrid ;
        nGrid = RegularGrid_NumberOfGridPoints ( grid ) ;
        if ( selection == NULL ) n = nGrid ;
        else                     n = selection->capacity ;
        if ( ( selection != NULL ) && ( Selection_UpperBound ( selection ) > nGrid ) ) Status_Set ( status, Status_IndexOutOfRange      ) ;
        else if ( ( grid->ndimensions != 3 ) || ( Coordinates3_Rows ( self ) != n )  ) Status_Set ( status, Status_NonConformableArrays ) ; /* . < n? */
        else
        {
            auto IntegerArray1D *indices = IntegerArray1D_AllocateWithExtent ( grid->ndimensions, status ) ; IntegerArray1D_Set ( indices, 0 ) ;
            if ( indices != NULL )
            {
                /* . Do all points. */
                if ( selection == NULL )
                {
                    auto Integer d, g, i, n ;
                    Array1D_Item ( indices, grid->ndimensions-1 ) = -1 ;
                    for ( g = n = 0 ; g < nGrid ; g++ )
                    {
                        for ( d = grid->ndimensions-1 ; d >= 0 ; d-- )
                        {
                            i = Array1D_Item ( indices, d ) + 1 ;
                            if ( i >= grid->dimensions[d].bins ) Array1D_Item ( indices, d ) = 0 ;
                            else { Array1D_Item ( indices, d ) = i ; break ; }
                        }
                        for ( d = 0 ; d < grid->ndimensions ; d++, n++ ) self->data[n] = ( ( Real ) Array1D_Item ( indices, d ) ) * grid->dimensions[d].binSize + grid->dimensions[d].midPointLower ;
                    }
                }
                /* . Do selected points. */
                else
                {
                    auto Integer d, i, n, s ;
                    for ( n = s = 0 ; s < selection->capacity ; s++ )
                    {
                        i = selection->indices[s] ;
                        for ( d = 0 ; d < grid->ndimensions ; d++ )
                        {
                            Array1D_Item ( indices, d ) = i / grid->dimensions[d].stride ;
                            i %= grid->dimensions[d].stride ;
                        }
                        for ( d = 0 ; d < grid->ndimensions ; d++, n++ ) self->data[n] = ( ( Real ) Array1D_Item ( indices, d ) ) * grid->dimensions[d].binSize + grid->dimensions[d].midPointLower ;
                    }
                }
                /* . Finish up. */
                IntegerArray1D_Deallocate ( &indices ) ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Gathering.
! . The operation is from a sparse other, indexed by selection, to a compact self.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Coordinates3_Gather ( Coordinates3 *self, const Coordinates3 *other, const Selection *selection )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        if ( selection == NULL ) Coordinates3_CopyTo ( other, self, NULL ) ;
        else
        {
            auto       Real    *newData ;
            auto const Real    *oldData ;
            auto       Integer  i, j, n ;

            for ( i = 0, n = 0, newData = Coordinates3_RowPointer ( self, 0 ) ; i < selection->capacity ; i++ )
            {
                for ( j = 0, oldData = Coordinates3_RowPointer ( other, selection->indices[i] ) ; j < self->extent1 ; j++, n++ ) newData[n] = oldData[j] ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Gathering with adding and scaling.
! . The operation is from a sparse other, indexed by selection, to a compact self.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Coordinates3_GatherAdd ( Coordinates3 *self, const Real alpha, const Coordinates3 *other, const Selection *selection )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        if ( selection == NULL ) Coordinates3_Add ( self, alpha, other, NULL ) ;
        else
        {
            auto       Real    *newData ;
            auto const Real    *oldData ;
            auto       Integer  i, j, n ;
/*
printf ( "\nOther:\n" ) ;
Coordinates3_Print ( other ) ;
printf ( "Selection:\n" ) ;
for ( i = 0 ; i < selection->capacity ; i++ ) printf ( "%5d", selection->indices[i] ) ;
printf ( "\nBefore:\n" ) ;
Coordinates3_Print ( self ) ;
*/
            for ( i = 0, n = 0, newData = Coordinates3_RowPointer ( self, 0 ) ; i < selection->capacity ; i++ )
            {
                for ( j = 0, oldData = Coordinates3_RowPointer ( other, selection->indices[i] ) ; j < self->extent1 ; j++, n++ ) newData[n] += ( alpha * oldData[j] ) ;
            }
/*
printf ( "\nAfter:\n" ) ;
Coordinates3_Print ( self ) ;
*/
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Identify points on a grid that are covered by the coordinates each of which has a specific radius.
! . An occupied grid point is defined in two ways:
!   1. A coordinate sphere overlaps the midpoint of a grid box (QMIDPOINTOVERLAP).
!   2. A coordinate sphere overlaps any part of a grid box (! QMIDPOINTOVERLAP).
! . Points with zero radii or less are skipped (assumed to have no overlap).
!---------------------------------------------------------------------------------------------------------------------------------*/
Status Coordinates3_IdentifyOccupiedGridPoints ( Coordinates3 *self, const RegularGrid *grid, const RealArray1D *radii, const Boolean QMIDPOINTOVERLAP, Selection **occupied )
{
    Status status = Status_OK ;
    if ( ( self != NULL ) && ( grid != NULL ) && ( radii != NULL ) && ( occupied != NULL ) )
    {
        (*occupied) = NULL ;
        if ( ( grid->ndimensions == 3 ) && ( Coordinates3_Rows ( self ) == View1D_Extent ( radii ) ) )
        {
            auto BooleanBlock *flags = NULL ;
            auto Integer nGrid ;
            nGrid = RegularGrid_NumberOfGridPoints ( grid ) ;
            flags = BooleanBlock_Allocate ( nGrid, NULL ) ;
            BooleanBlock_Set ( flags, False ) ;
            if ( flags != NULL )
            {
                IntegerArray1D *indices = NULL, *lower = NULL, *range = NULL ;
                /* . Allocate other arrays. */
                indices = IntegerArray1D_AllocateWithExtent ( grid->ndimensions, NULL ) ;
                lower   = IntegerArray1D_AllocateWithExtent ( grid->ndimensions, NULL ) ;
                range   = IntegerArray1D_AllocateWithExtent ( grid->ndimensions, NULL ) ;

                /* . Mid-point overlap. */
                if ( QMIDPOINTOVERLAP )
                {
                    auto Integer d, i, index = 0, g, l, nboxes, p, u ;
                    auto Real    c, radius, radius2, r2, *row ;
                    for ( p = 0 ; p < Coordinates3_Rows ( self ) ; p++ )
                    {
                        radius = Array1D_Item ( radii, p ) ;
                        if ( radius >= 0.0e+00 )
                        {
                            radius2 = radius * radius ;
                            row     = Coordinates3_RowPointer ( self, p ) ;
                            /* . Get the bounds for the search. */
                            for ( d = 0, nboxes = 1 ; d < grid->ndimensions ; d++ )
                            {
                                /* . Lower. */
                                c = row[d] - radius - grid->dimensions[d].lower ;
                                l = ( Integer ) floor ( c / grid->dimensions[d].binSize ) ;
                                if ( c > ( ( Real ) l + 0.5e+00 ) * grid->dimensions[d].binSize ) l += 1 ;
                                if ( l < 0 ) l = 0 ;
                                else if ( l >= grid->dimensions[d].bins ) l = grid->dimensions[d].bins - 1 ;
                                /* . Upper. */
                                c = row[d] + radius - grid->dimensions[d].lower ;
                                u = ( Integer ) floor ( c / grid->dimensions[d].binSize ) ;
                                if ( c < ( ( Real ) u + 0.5e+00 ) * grid->dimensions[d].binSize ) u -= 1 ;
                                if ( u < 0 ) u = 0 ;
                                else if ( u >= grid->dimensions[d].bins ) u = grid->dimensions[d].bins - 1 ;
                                /* . Set range. */
                                Array1D_Item ( lower, d ) = l ;
                                Array1D_Item ( range, d ) = Maximum ( u - l + 1, 0 ) ;
                                nboxes *= Array1D_Item ( range, d ) ;
                            }
                            /* . Loop over boxes. */
                            IntegerArray1D_Set ( indices, 0 ) ;
                            Array1D_Item ( indices, grid->ndimensions-1 ) = -1 ;
                            for ( g = 0 ; g < nboxes ; g++ )
                            {
                                for ( d = grid->ndimensions-1 ; d >= 0 ; d-- )
                                {
                                    i = Array1D_Item ( indices, d ) + 1 ;
                                    if ( i >= Array1D_Item ( range, d ) ) Array1D_Item ( indices, d ) = 0 ;
                                    else { Array1D_Item ( indices, d ) = i ; break ; }
                                }
                                for ( d = index = 0, r2 = 0.0e+00 ; d < grid->ndimensions ; d++ )
                                {
                                    index += ( Array1D_Item ( indices, d ) + Array1D_Item ( lower, d ) ) * grid->dimensions[d].stride ;
                                    r2    += pow ( ( ( ( Real ) ( Array1D_Item ( indices, d ) + Array1D_Item ( lower, d ) ) ) * grid->dimensions[d].binSize + grid->dimensions[d].midPointLower - row[d] ), 2 ) ;
                                }
# ifdef DEBUG
if ( index >= nGrid ) { printf ( "ERROR> Grid point index out of range %d %d\n", index, nGrid ) ; exit ( 9999 ) ; }
# endif
                                if ( r2 <= radius2 ) Block_Item ( flags, index ) = True ;
                            }
                        }
                    }
                }
                /* . Overlap of any part of a box. */
                else
                {
                    auto Integer d, i, index, g, l, nboxes, p, u ;
                    auto Real    c, cl, cu, radius, radius2, r2, *row ;
                    for ( p = 0 ; p < Coordinates3_Rows ( self ) ; p++ )
                    {
                        radius = Array1D_Item ( radii, p ) ;
                        if ( radius >= 0.0e+00 )
                        {
                            radius2 = radius * radius ;
                            row     = Coordinates3_RowPointer ( self, p ) ;
                            /* . Get the bounds for the search. */
                            for ( d = 0, nboxes = 1 ; d < grid->ndimensions ; d++ )
                            {
                                /* . Lower. */
                                c = row[d] - radius - grid->dimensions[d].lower ;
                                l = ( Integer ) floor ( c / grid->dimensions[d].binSize ) ;
                                if ( l < 0 ) l = 0 ;
                                else if ( l >= grid->dimensions[d].bins ) l = grid->dimensions[d].bins - 1 ;
                                /* . Upper. */
                                c = row[d] + radius - grid->dimensions[d].lower ;
                                u = ( Integer ) floor ( c / grid->dimensions[d].binSize ) ;
                                if ( u < 0 ) u = 0 ;
                                else if ( u >= grid->dimensions[d].bins ) u = grid->dimensions[d].bins - 1 ;
                                /* . Set range. */
                                Array1D_Item ( lower, d ) = l ;
                                Array1D_Item ( range, d ) = Maximum ( u - l + 1, 0 ) ;
                                nboxes *= Array1D_Item ( range, d ) ;
                            }
                            /* . Loop over boxes. */
                            IntegerArray1D_Set ( indices, 0 ) ;
                            Array1D_Item ( indices, grid->ndimensions-1 ) = -1 ;
                            for ( g = 0 ; g < nboxes ; g++ )
                            {
                                for ( d = grid->ndimensions-1 ; d >= 0 ; d-- )
                                {
                                    i = Array1D_Item ( indices, d ) + 1 ;
                                    if ( i >= Array1D_Item ( range, d ) ) Array1D_Item ( indices, d ) = 0 ;
                                    else { Array1D_Item ( indices, d ) = i ; break ; }
                                }
                                for ( d = index = 0, r2 = 0.0e+00 ; d < grid->ndimensions ; d++ )
                                {
                                    index += ( Array1D_Item ( indices, d ) + Array1D_Item ( lower, d ) ) * grid->dimensions[d].stride ;
                                    cl     = ( ( Real ) ( Array1D_Item ( indices, d ) + Array1D_Item ( lower, d ) ) ) * grid->dimensions[d].binSize + grid->dimensions[d].lower ;
                                    cu     = cl + grid->dimensions[d].binSize ;
                                    if      ( row[d] < cl ) r2 += pow ( cl - row[d], 2 ) ;
                                    else if ( row[d] > cu ) r2 += pow ( row[d] - cu, 2 ) ;
                                }
# ifdef DEBUG
if ( index >= nGrid ) { printf ( "ERROR> Grid point index out of range %d %d\n", index, nGrid ) ; exit ( 9999 ) ; }
# endif
                                if ( r2 <= radius2 ) Block_Item ( flags, index ) = True ;
                            }
                        }
                    }
                }
                /* . Finish up. */
                IntegerArray1D_Deallocate ( &indices ) ;
                IntegerArray1D_Deallocate ( &lower   ) ;
                IntegerArray1D_Deallocate ( &range   ) ;
                (*occupied) = Selection_FromBooleans ( flags->capacity, flags->items, &status ) ;
            }
            BooleanBlock_Deallocate ( &flags ) ;
        }
        else status = Status_NonConformableArrays ;
    }
    return status ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the inertia matrix.
! . The matrix are moved to their center for calculation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Coordinates3_InertiaMatrix ( const Coordinates3 *self, const Selection *selection, const RealArray1D *weights, SymmetricMatrix *inertia )
{
    if ( ( self != NULL ) && ( inertia != NULL ) )
    {
        auto Integer  i, s ;
        auto Real     cx, cy, cz, w, x, y, z, xx = 0.0e+00, xy = 0.0e+00, xz = 0.0e+00, yy = 0.0e+00, yz = 0.0e+00, zz = 0.0e+00 ;
        auto Vector3 *center = NULL ;
        /* . Center. */
        Coordinates3_Center ( self, selection, weights, &center ) ;
        cx = Vector3_Item ( center, 0 ) ;
        cy = Vector3_Item ( center, 1 ) ;
        cz = Vector3_Item ( center, 2 ) ;
        /* . The inertia matrix. */
        /* . The four cases. */
        if ( selection == NULL )
        {
            if ( weights == NULL )
            {
                for ( i = 0 ; i < self->extent0 ; i++ )
                {
	            Coordinates3_GetRow ( self, i, x, y, z ) ;
	            x  -= cx ;
	            y  -= cy ;
	            z  -= cz ;
                    xx += x * x ;
                    xy += x * y ;
                    xz += x * z ;
                    yy += y * y ;
                    yz += y * z ;
                    zz += z * z ;
                }
            }
            else
            {
                for ( i = 0 ; i < self->extent0 ; i++ )
                {
                    w = Array1D_Item ( weights, i ) ;
	            Coordinates3_GetRow ( self, i, x, y, z ) ;
	            x  -= cx ;
	            y  -= cy ;
	            z  -= cz ;
                    xx += w * x * x ;
                    xy += w * x * y ;
                    xz += w * x * z ;
                    yy += w * y * y ;
                    yz += w * y * z ;
                    zz += w * z * z ;
                }
            }
        }
        else
        {
            if ( weights == NULL )
            {
                for ( s = 0 ; s < selection->capacity ; s++ )
                {
                    i = selection->indices[s] ;
	            Coordinates3_GetRow ( self, i, x, y, z ) ;
	            x  -= cx ;
	            y  -= cy ;
	            z  -= cz ;
                    xx += x * x ;
                    xy += x * y ;
                    xz += x * z ;
                    yy += y * y ;
                    yz += y * z ;
                    zz += z * z ;
                }
            }
            else
            {
                for ( s = 0 ; s < selection->capacity ; s++ )
                {
                    i = selection->indices[s] ;
                    w = Array1D_Item ( weights, i ) ;
	            Coordinates3_GetRow ( self, i, x, y, z ) ;
	            x  -= cx ;
	            y  -= cy ;
	            z  -= cz ;
                    xx += w * x * x ;
                    xy += w * x * y ;
                    xz += w * x * z ;
                    yy += w * y * y ;
                    yz += w * y * z ;
                    zz += w * z * z ;
                }
            }
        }
        SymmetricMatrix_Item ( inertia, 0, 0 ) = yy + zz ; /* xx */
        SymmetricMatrix_Item ( inertia, 1, 0 ) =    - xy ; /* xy */
        SymmetricMatrix_Item ( inertia, 1, 1 ) = xx + zz ; /* yy */
        SymmetricMatrix_Item ( inertia, 2, 0 ) =    - xz ; /* xz */
        SymmetricMatrix_Item ( inertia, 2, 1 ) =    - yz ; /* yz */
        SymmetricMatrix_Item ( inertia, 2, 2 ) = xx + yy ; /* zz */
        /* . Finish up. */
        Vector3_Deallocate ( &center ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set up a conforming grid given a set of coordinates and another grid.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Coordinates3_MakeConformingGrid ( const Coordinates3 *self, const Selection    *andSelection   ,
                                                                 const RegularGrid  *grid           ,
                                                                 RegularGrid       **conformingGrid ,
                                                                 IntegerArray1D    **offSet         ,
                                                                 Status             *status         )
{
# ifdef DEBUGPRINTING
printf ( "\nConforming Grid Pointers: %p %p %p %p %p\n", self, andSelection, grid, conformingGrid, offSet ) ;
# endif
    if ( ( self != NULL ) && ( grid != NULL ) && ( conformingGrid != NULL ) && ( offSet != NULL ) )
    {
        auto Integer         d, lowerN, stride, upperN[3] ;
        auto Real            gridSize ;
        auto IntegerArray1D *newOrigin = NULL ;
        auto RegularGrid    *newGrid   = NULL ;
        auto Vector3        *lower     = NULL, *upper = NULL ;
        /* . Initialization. */
        (*conformingGrid) = NULL ;
        (*offSet)         = NULL ;
# ifdef DEBUGPRINTING
printf ( "\nConforming Grid Data: %d %d %d %d\n", self->extent0, self->extent1, Selection_Capacity ( andSelection ), RegularGrid_NumberOfGridPoints ( grid ) ) ;
# endif
        /* . Allocation. */
        newGrid   = RegularGrid_Allocate    ( 3, status ) ;
        newOrigin = IntegerArray1D_AllocateWithExtent ( 3, status ) ;
        lower     = Vector3_Allocate ( ) ;
        upper     = Vector3_Allocate ( ) ;
        if ( ( newGrid != NULL ) && ( newOrigin != NULL ) && ( lower != NULL ) && ( upper != NULL ) )
        {
            /* . Find the upper and lower limits of the coordinates. */
            Coordinates3_EnclosingOrthorhombicBox ( self, andSelection, NULL, lower, upper ) ;
            Vector3_Add ( upper, 1.0e+00, lower, NULL ) ;

            /* . Find the cell indices of lower which gives the offset. */
            IntegerArray1D_Set ( newOrigin, 0 ) ;
            RegularGrid_FindCellIndicesOfPoint ( grid, Vector3_Data ( lower ), False, Array1D_Data ( newOrigin ), NULL ) ;

            /* . Find the cell indices of upper. */
            RegularGrid_FindCellIndicesOfPoint ( grid, Vector3_Data ( upper ), False, upperN, NULL ) ;

            /* . Basic data. */
# ifdef DEBUGPRINTING
printf ( "\nConforming Grid Data:\n" ) ;
# endif
            for ( d = 0 ; d < 3 ; d++ )
            {
                gridSize = grid->dimensions[d].binSize ;
                lowerN   = Array1D_Item ( newOrigin, d ) ;
                newGrid->dimensions[d].bins    = upperN[d] - lowerN + 1 ;
                newGrid->dimensions[d].binSize = gridSize ;
                newGrid->dimensions[d].lower   = grid->dimensions[d].lower + ( Real ) ( lowerN ) * gridSize ;
# ifdef DEBUGPRINTING
printf ( "%6d %6d %6d %6d %10.3f %10.3f %10.3f %10.3f %10.3f\n", d, lowerN, upperN[d], newGrid->dimensions[d].bins, gridSize,
grid->dimensions[d].lower, newGrid->dimensions[d].lower, Vector3_Item ( lower, d ), Vector3_Item ( upper, d ) ) ;
# endif
            }

            /* . Remaining data (do in reverse order for the strides). */
            stride = 1 ;
            for ( d = 2 ; d >= 0 ; d-- )
            {
                newGrid->dimensions[d].midPointLower = newGrid->dimensions[d].lower + 0.5e+00 * newGrid->dimensions[d].binSize ;
                newGrid->dimensions[d].upper         = newGrid->dimensions[d].lower + ( Real ) ( newGrid->dimensions[d].bins ) * newGrid->dimensions[d].binSize ;
                newGrid->dimensions[d].stride        = stride ;
                stride *= newGrid->dimensions[d].bins ;
            }

            /* . Finish up. */
            Vector3_Deallocate ( &lower ) ;
            Vector3_Deallocate ( &upper ) ;
            (*conformingGrid) = newGrid   ;
            (*offSet)         = newOrigin ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a conforming grid and occupancy suitable for pairlist generation.
!---------------------------------------------------------------------------------------------------------------------------------*/
/*
! . This should be improved by allowing mixed non-periodic/periodic grids and the option of having parts of a periodic grid (rather than the entire thing).
*/
void Coordinates3_MakeConformingGridAndOccupancy ( const Coordinates3    *self           ,
                                                   const Selection       *andSelection   ,
                                                   const RegularGrid     *grid           ,
                                                   RegularGrid          **conformingGrid ,
                                                   RegularGridOccupancy **occupancy      ,
                                                   IntegerArray1D       **offSet         ,
                                                   Status                *status         )
{
    if ( ( self != NULL ) && ( self->extent0 > 0 ) && ( grid != NULL ) && ( conformingGrid != NULL ) && ( occupancy != NULL ) && ( offSet != NULL ) )
    {
        if ( RegularGrid_IsPeriodic ( grid ) )
        {
            (*conformingGrid) = RegularGrid_Clone       ( grid, status ) ;
            (*offSet        ) = IntegerArray1D_AllocateWithExtent (    3, status ) ;
            IntegerArray1D_Set ( (*offSet), 0 ) ;
        }
        else Coordinates3_MakeConformingGrid ( self, andSelection, grid, conformingGrid, offSet, status ) ;
        (*occupancy) = RegularGridOccupancy_FromGridAndPoints ( (*conformingGrid), self, status ) ;
    }
    else
    {
        if ( conformingGrid != NULL ) (*conformingGrid) = NULL ;
        if ( occupancy      != NULL ) (*occupancy)      = NULL ;
        if ( offSet         != NULL ) (*offSet)         = NULL ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set up grid from a set of coordinates.
!---------------------------------------------------------------------------------------------------------------------------------*/
RegularGrid *Coordinates3_MakeGrid ( const Coordinates3 *self, const Selection *andSelection, const Real gridSize, Status *status )
{
    RegularGrid *grid = NULL ;
    if ( self != NULL )
    {
        auto Integer  d, stride ;
        auto Real     cells, extent, lower, midPoint, upper ;
        auto Vector3 *extents = NULL, *origin = NULL ;

        /* . Allocation. */
        extents = Vector3_Allocate ( ) ;
        origin  = Vector3_Allocate ( ) ;

        /* . Find the extents of the coordinates. */
        Coordinates3_EnclosingOrthorhombicBox ( self, andSelection, NULL, origin, extents ) ;

        /* . Allocation. */
        grid = RegularGrid_Allocate ( 3, status ) ;

        /* . Basic data. */
        for ( d = 0 ; d < 3 ; d++ )
        {
            extent   = Vector3_Item ( extents, d ) ;
            lower    = Vector3_Item ( origin , d ) ;
            upper    = lower + extent ;
            midPoint = 0.5e+00 * ( lower + upper ) ;
            cells    = Maximum ( ceil ( extent / gridSize ), 1 ) ;
            grid->dimensions[d].bins    = ( Integer ) cells ;
            grid->dimensions[d].binSize = gridSize ;
            grid->dimensions[d].lower   = midPoint - 0.5e+00 * ( cells * gridSize ) ;
        }

        /* . Remaining data (do in reverse order for the strides). */
        stride = 1 ;
        for ( d = 2 ; d >= 0 ; d-- )
        {
            grid->dimensions[d].midPointLower = grid->dimensions[d].lower + 0.5e+00 * grid->dimensions[d].binSize ;
            grid->dimensions[d].upper         = grid->dimensions[d].lower + ( Real ) ( grid->dimensions[d].bins ) * grid->dimensions[d].binSize ;
            grid->dimensions[d].stride        = stride ;
            stride *= grid->dimensions[d].bins ;
        }

        /* . Finish up. */
        Vector3_Deallocate ( &extents ) ;
        Vector3_Deallocate ( &origin  ) ;
    }
    return grid ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a grid and occupancy suitable for pairlist generation.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Note that the andSelection here can be a superset of the andSelection passed to the generation procedures. */
void Coordinates3_MakeGridAndOccupancy ( const Coordinates3    *self         ,
                                         const Selection       *andSelection ,
                                         const Real             gridSize     ,
                                         RegularGrid          **grid         ,
                                         RegularGridOccupancy **occupancy    ,
                                         Status                *status       )
{
    if ( ( self != NULL ) && ( self->extent0 > 0 ) )
    {
        (*grid)      = Coordinates3_MakeGrid ( self, andSelection, gridSize, status ) ;
        (*occupancy) = RegularGridOccupancy_FromGridAndPoints ( (*grid), self, status ) ;
    }
    else { (*grid) = NULL ; (*occupancy) = NULL ; }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Make a periodic grid and occupancy suitable for pairlist generation.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . Note that the andSelection here can be a superset of the andSelection passed to the generation procedures. */
void Coordinates3_MakePeriodicGridAndOccupancy ( const Coordinates3    *self      ,
                                                 const Vector3         *boxSize   ,
                                                 const Real             gridSize  ,
                                                 RegularGrid          **grid      ,
                                                 RegularGridOccupancy **occupancy ,
                                                 Status                *status    )
{
    if ( ( self != NULL ) && ( self->extent0 > 0 ) )
    {
        (*grid)      = RegularGrid_MakePeriodicGrid3 ( boxSize, gridSize, status ) ;
        (*occupancy) = RegularGridOccupancy_FromGridAndPoints ( (*grid), self, status ) ;
    }
    else { (*grid) = NULL ; (*occupancy) = NULL ; }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the moments of inertia.
! . The axes are not unique (there is a sign dependence).
!---------------------------------------------------------------------------------------------------------------------------------*/
void Coordinates3_MomentsOfInertia ( const Coordinates3 *self, const Selection *selection, const RealArray1D *weights, Vector3 *moments, Matrix33 *axes )
{
    /* . Initialization. */
    Vector3_Set  ( moments, 0.0e+00  ) ;
    Matrix33_Set ( axes,    BadValue ) ;
    if ( ( self != NULL ) && ( moments != NULL ) )
    {
        auto Real det = 0.0e+00 ;
        auto SymmetricMatrix *inertia = NULL ;
        /* . Get the moments and axes. */
        inertia = SymmetricMatrix_AllocateWithExtent ( 3, NULL ) ;
        Coordinates3_InertiaMatrix ( self, selection, weights, inertia ) ;
        SymmetricMatrix_EigenvaluesSolve     ( inertia, False, 0, 3, moments, axes, False, NULL ) ;
        SymmetricMatrix_Deallocate           ( &inertia ) ;
        /* . Ensure that the determinant of "axes" is positive. */
        if ( axes != NULL )
        {
            det = Matrix33_Determinant ( axes ) ;
            if ( det < 0.0e+00 )
            {
                auto Real v ;
                auto Integer    i ;
                for ( i = 0 ; i < 3 ; i++ )
                {
                    v = Matrix33_Item ( axes, i, 0 ) ;
                    Matrix33_Item ( axes, i, 0 ) = -1.0e+00 * v ;
                }
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the radius of gyration.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Coordinates3_RadiusOfGyration ( const Coordinates3 *self, const Selection *selection, const RealArray1D *weights )
{
    Real rgyr = 0.0e+00 ;
    if ( self != NULL )
    {
        auto Real cx, cy, cz, wTotal = 0.0e+00, w, x, y, z ;
        auto Integer    i, s ;
        auto Vector3 *center = NULL ;
        /* . Center. */
        Coordinates3_Center ( self, selection, weights, &center ) ;
        cx = Vector3_Item ( center, 0 ) ;
        cy = Vector3_Item ( center, 1 ) ;
        cz = Vector3_Item ( center, 2 ) ;
        /* . Radius of gyration. */
        /* . The four cases. */
        if ( selection == NULL )
        {
            if ( weights == NULL )
            {
                for ( i = 0 ; i < self->extent0 ; i++ )
                {
                    Coordinates3_GetRow ( self, i, x, y, z ) ;
                    rgyr  += pow ( x - cx, 2 ) + pow ( y - cy, 2 ) + pow ( z - cz, 2 ) ;
                    wTotal += 1.0e+00 ;
                }
            }
            else
            {
                for ( i = 0 ; i < self->extent0 ; i++ )
                {
                    Coordinates3_GetRow ( self, i, x, y, z ) ;
                    w      = Array1D_Item ( weights, i ) ;
                    rgyr  += w * ( pow ( x - cx, 2 ) + pow ( y - cy, 2 ) + pow ( z - cz, 2 ) ) ;
                    wTotal += w ;
                }
            }
        }
        else
        {
            if ( weights == NULL )
            {
                for ( s = 0 ; s < selection->capacity ; s++ )
                {
                    i = selection->indices[s] ;
                    Coordinates3_GetRow ( self, i, x, y, z ) ;
                    rgyr  += pow ( x - cx, 2 ) + pow ( y - cy, 2 ) + pow ( z - cz, 2 ) ;
                    wTotal += 1.0e+00 ;
                }
            }
            else
            {
                for ( s = 0 ; s < selection->capacity ; s++ )
                {
                    i = selection->indices[s] ;
                    Coordinates3_GetRow ( self, i, x, y, z ) ;
                    w      = Array1D_Item ( weights, i ) ;
                    rgyr  += w * ( pow ( x - cx, 2 ) + pow ( y - cy, 2 ) + pow ( z - cz, 2 ) ) ;
                    wTotal += w ;
                }
            }
        }
        /* . Scaling. */
        if ( wTotal == 0.0e+00 ) rgyr  = 0.0e+00 ;
        else                    rgyr /= wTotal   ;
        if ( rgyr   > 0.0e+00 ) rgyr  = sqrt ( rgyr ) ;
        else                    rgyr  = 0.0e+00 ;
        /* . Finish up. */
        Vector3_Deallocate ( &center ) ;
    }
    return rgyr ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate the root mean square deviation deviation between two data sets.
!---------------------------------------------------------------------------------------------------------------------------------*/
Real Coordinates3_RootMeanSquareDeviation ( const Coordinates3 *self, const Coordinates3 *other, const Selection *selection, const RealArray1D *weights )
{
    Real rmsd = 0.0e+00 ;
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Real wTotal = 0.0e+00, w, x1, x2, y1, y2, z1, z2 ;
        auto Integer    i, s ;
        /* . The RMS deviation. */
        /* . The four cases. */
        if ( selection == NULL )
        {
            if ( weights == NULL )
            {
                for ( i = 0 ; i < self->extent0 ; i++ )
                {
                    Coordinates3_GetRow ( self,  i, x1, y1, z1 ) ;
                    Coordinates3_GetRow ( other, i, x2, y2, z2 ) ;
                    rmsd   += pow ( x1 - x2, 2 ) + pow ( y1 - y2, 2 ) + pow ( z1 - z2, 2 ) ;
                    wTotal += 1.0e+00 ;
                }
            }
            else
            {
                for ( i = 0 ; i < self->extent0 ; i++ )
                {
                    Coordinates3_GetRow ( self,  i, x1, y1, z1 ) ;
                    Coordinates3_GetRow ( other, i, x2, y2, z2 ) ;
                    w       = Array1D_Item ( weights, i ) ;
                    rmsd   += w * ( pow ( x1 - x2, 2 ) + pow ( y1 - y2, 2 ) + pow ( z1 - z2, 2 ) ) ;
                    wTotal += w ;
                }
            }
        }
        else
        {
            if ( weights == NULL )
            {
                for ( s = 0 ; s < selection->capacity ; s++ )
                {
                    i = selection->indices[s] ;
                    Coordinates3_GetRow ( self,  i, x1, y1, z1 ) ;
                    Coordinates3_GetRow ( other, i, x2, y2, z2 ) ;
                    rmsd   += pow ( x1 - x2, 2 ) + pow ( y1 - y2, 2 ) + pow ( z1 - z2, 2 ) ;
                    wTotal += 1.0e+00 ;
                }
            }
            else
            {
                for ( s = 0 ; s < selection->capacity ; s++ )
                {
                    i = selection->indices[s] ;
                    Coordinates3_GetRow ( self,  i, x1, y1, z1 ) ;
                    Coordinates3_GetRow ( other, i, x2, y2, z2 ) ;
                    w       = Array1D_Item ( weights, i ) ;
                    rmsd   += w * ( pow ( x1 - x2, 2 ) + pow ( y1 - y2, 2 ) + pow ( z1 - z2, 2 ) ) ;
                    wTotal += w ;
                }
            }
        }
        /* . Scaling. */
        if ( wTotal == 0.0e+00 ) rmsd  = 0.0e+00 ;
        else                     rmsd /= wTotal  ;
        if ( rmsd   > 0.0e+00  ) rmsd  = sqrt ( rmsd ) ;
        else                     rmsd  = 0.0e+00 ;
    }
    return rmsd ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Rotation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Coordinates3_Rotate ( Coordinates3 *self, const Matrix33 *rotation, const Selection *selection )
{
    if ( ( self != NULL ) && ( rotation != NULL ) )
    {
        auto Real r00, r01, r02, r10, r11, r12, r20, r21, r22, x0, x1, y0, y1, z0, z1 ;
        auto Integer    i, s ;
        r00 = Matrix33_Item ( rotation, 0, 0 ) ;
        r01 = Matrix33_Item ( rotation, 0, 1 ) ;
        r02 = Matrix33_Item ( rotation, 0, 2 ) ;
        r10 = Matrix33_Item ( rotation, 1, 0 ) ;
        r11 = Matrix33_Item ( rotation, 1, 1 ) ;
        r12 = Matrix33_Item ( rotation, 1, 2 ) ;
        r20 = Matrix33_Item ( rotation, 2, 0 ) ;
        r21 = Matrix33_Item ( rotation, 2, 1 ) ;
        r22 = Matrix33_Item ( rotation, 2, 2 ) ;
        if ( selection == NULL )
        {
            for ( i = 0 ; i < self->extent0 ; i++ )
            {
                Coordinates3_GetRow ( self, i, x0, y0, z0 ) ;
                x1 = r00 * x0 + r01 * y0 + r02 * z0 ;
                y1 = r10 * x0 + r11 * y0 + r12 * z0 ;
                z1 = r20 * x0 + r21 * y0 + r22 * z0 ;
                Coordinates3_SetRow ( self, i, x1, y1, z1 ) ;
            }
        }
        else
        {
            for ( s = 0 ; s < selection->capacity ; s++ )
            {
                i = selection->indices[s] ;
                Coordinates3_GetRow ( self, i, x0, y0, z0 ) ;
                x1 = r00 * x0 + r01 * y0 + r02 * z0 ;
                y1 = r10 * x0 + r11 * y0 + r12 * z0 ;
                z1 = r20 * x0 + r21 * y0 + r22 * z0 ;
                Coordinates3_SetRow ( self, i, x1, y1, z1 ) ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Calculate a set of orthonormal rotation and translation vectors.
! . The input flags indicate which vectors should be calculated.
! . If the weights array is present, the vectors are determined in terms
! . of the weighted matrix ( i.e. sqrt ( weights ) ).
! . |dimension| is the requested extent of each vector. It is only used if
! . it is bigger than self->extent0 * self->extent1.
!---------------------------------------------------------------------------------------------------------------------------------*/
Integer Coordinates3_RotationTranslationVectors ( const Coordinates3 *self    ,
                                                  const RealArray1D  *weights ,
                                                  const Boolean       QRx     ,
                                                  const Boolean       QRy     ,
                                                  const Boolean       QRz     ,
                                                  const Boolean       QTx     ,
                                                  const Boolean       QTy     ,
                                                  const Boolean       QTz     ,
                                                        RealArray2D  *vectors ,
                                                        Status       *status  )
{
    Integer nVectors = 0 ;
    /* . Determine the number of vectors. */
    if ( QRx ) nVectors++ ; if ( QRy ) nVectors++ ; if ( QRz ) nVectors++ ;
    if ( QTx ) nVectors++ ; if ( QTy ) nVectors++ ; if ( QTz ) nVectors++ ;
    /* . Initialization. */
    RealArray2D_Set ( vectors, 0.0e+00 ) ;
    if ( ( self != NULL ) && ( nVectors > 0 ) && ( vectors != NULL ) && Status_IsOK ( status ) )
    {
        if ( ( View2D_Rows ( vectors ) >= Coordinates3_Size ( self ) ) && ( View2D_Columns ( vectors ) >= nVectors ) )
        {
            auto Integer      i, iatom, inc ;
            auto Real         cx = 0.0e+00, cy = 0.0e+00, cz = 0.0e+00, w, x, y, z ;
            auto RealArray1D *wts    = NULL ;
            auto Vector3     *center = NULL ;
            /* . Get the weights. */
            wts = RealArray1D_AllocateWithExtent ( self->extent0, NULL ) ;
            if ( weights == NULL ) RealArray1D_Set    ( wts, 1.0e+00 ) ;
            else                   RealArray1D_CopyTo ( weights, wts, NULL ) ;
            /* . Get the center. */
            if ( QRx || QRy || QRz )
            {
                Coordinates3_Center ( self, NULL, wts, &center ) ;
                cx = Vector3_Item ( center, 0 ) ;
                cy = Vector3_Item ( center, 1 ) ;
                cz = Vector3_Item ( center, 2 ) ;
            }
            /* . Calculate the square root of the weights. */
            if ( weights != NULL )
            {
                for ( i = 0 ; i < wts->extent ; i++ )
                {
                    w = sqrt ( Array1D_Item ( wts, i ) ) ;
                    Array1D_Item ( wts, i ) = w ;
                }
            }
            /* . Calculate the vectors. */
            for ( iatom = 0 ; iatom < self->extent0 ; iatom++ )
            {
                inc = 0 ;
                w   = Array1D_Item ( wts, iatom ) ;
                Coordinates3_GetRow ( self, iatom, x, y, z ) ;
                if ( QTx ) { RealArray2D_SetItem ( vectors, 3*iatom    , inc, w, NULL ) ; inc++ ; }
                if ( QTy ) { RealArray2D_SetItem ( vectors, 3*iatom + 1, inc, w, NULL ) ; inc++ ; }
                if ( QTz ) { RealArray2D_SetItem ( vectors, 3*iatom + 2, inc, w, NULL ) ; inc++ ; }
                if ( QRx ) { RealArray2D_SetItem ( vectors, 3*iatom + 2, inc,   w * ( y - cy ), NULL ) ;
                             RealArray2D_SetItem ( vectors, 3*iatom + 1, inc, - w * ( z - cz ), NULL ) ; inc++ ; }
                if ( QRy ) { RealArray2D_SetItem ( vectors, 3*iatom,     inc,   w * ( z - cz ), NULL ) ;
                             RealArray2D_SetItem ( vectors, 3*iatom + 2, inc, - w * ( x - cx ), NULL ) ; inc++ ; }
                if ( QRz ) { RealArray2D_SetItem ( vectors, 3*iatom + 1, inc,   w * ( x - cx ), NULL ) ;
                             RealArray2D_SetItem ( vectors, 3*iatom,     inc, - w * ( y - cy ), NULL ) ; }
            }
            /* . Create a linearly independent set of orthonormal vectors. */
            nVectors = RealArray2D_GramSchmidtOrthogonalize ( vectors, NULL, NULL, NULL, NULL ) ;
            /* . Finish up. */
            RealArray1D_Deallocate ( &wts    ) ;
            Vector3_Deallocate     ( &center ) ;
        }
        else Status_Set ( status, Status_NonConformableArrays ) ;
    }
    else nVectors = 0 ;
    return nVectors ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Scale each row by a unique scaling factor.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Coordinates3_ScaleRows ( Coordinates3 *self, const RealArray1D *rowScalingFactors, Status *status )
{
    if ( self != NULL )
    {
        auto Integer i, n ;
	auto Real    f, x, y, z ;
	if ( self->extent0 != View1D_Extent ( rowScalingFactors ) )
	{
	    n = Minimum ( self->extent0, View1D_Extent ( rowScalingFactors ) ) ;
	    Status_Set ( status, Status_NonConformableArrays ) ;
	}
	else n = self->extent0 ;
	for ( i = 0 ; i < n ; i++ )
	{
	    f = Array1D_Item ( rowScalingFactors, i ) ;
	    Coordinates3_GetRow ( self, i,   x,   y,   z ) ;
	    Coordinates3_SetRow ( self, i, f*x, f*y, f*z ) ;
	}
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Set the values of certain rows.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Coordinates3_SetByRow (       Coordinates3 *self      ,
                             const Selection    *selection ,
                             const Real          alpha     ,
                                   Status       *status    )
{
    if ( ( self != NULL ) && Status_IsOK ( status ) )
    {
        if ( selection == NULL ) Coordinates3_Set ( self, alpha ) ;
        else
        {
            if ( Selection_UpperBound ( selection ) <= Coordinates3_Rows ( self ) )
            {
                auto Integer i, j, s ;
                for ( i = 0 ; i < selection->capacity ; i++ )
                {
                    s = selection->indices[i] ;
                    for ( j = 0 ; j < self->extent1 ; j++ ) Array2D_Item ( self, s, j ) = alpha ;
                }
            }
            else Status_Set ( status, Status_IndexOutOfRange ) ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Scattering.
! . The operation is from a compact self to a sparse other, indexed by selection.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Coordinates3_Scatter ( const Coordinates3 *self, Coordinates3 *other, const Selection *selection )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        if ( selection == NULL ) Coordinates3_CopyTo ( self, other, NULL ) ;
        else
        {
            auto       Real    *newData ;
            auto const Real    *oldData ;
            auto       Integer  i, j, n ;
            for ( i = 0, n = 0, oldData = Coordinates3_RowPointer ( self, 0 ) ; i < selection->capacity ; i++ )
            {
                for ( j = 0, newData = Coordinates3_RowPointer ( other, selection->indices[i] ) ; j < self->extent1 ; j++, n++ ) newData[j] = oldData[n] ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Scattering with adding and scaling.
! . The operation is from a compact self to a sparse other, indexed by selection.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Coordinates3_ScatterAdd ( const Coordinates3 *self, const Real alpha, Coordinates3 *other, const Selection *selection )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        if ( selection == NULL ) Coordinates3_Add ( other, alpha, self, NULL ) ;
        else
        {
            auto       Real    *newData ;
            auto const Real    *oldData ;
            auto       Integer  i, j, n ;
            for ( i = 0, n = 0, oldData = Coordinates3_RowPointer ( self, 0 ) ; i < selection->capacity ; i++ )
            {
                for ( j = 0, newData = Coordinates3_RowPointer ( other, selection->indices[i] ) ; j < self->extent1 ; j++, n++ ) newData[j] += ( alpha * oldData[n] ) ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Superimpose one coordinate set (self) upon another (other).
! . |selection| determines the calculation of the transformation but not
! . which rows are transformed.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Coordinates3_Superimpose ( Coordinates3 *self, const Coordinates3 *other, const Selection *selection, const RealArray1D *weights, Matrix33 *rotation, Vector3 *translation )
{
    if ( ( self != NULL ) && ( other != NULL ) )
    {
        auto Integer          i, s ;
        auto Real             cx, cy, cz, diagonal, w, x1, x2, y1, y2, z1, z2, xx = 0.0e+00, xy = 0.0e+00, xz = 0.0e+00,
                                                                               yx = 0.0e+00, yy = 0.0e+00, yz = 0.0e+00,
                                                                               zx = 0.0e+00, zy = 0.0e+00, zz = 0.0e+00 ;
        auto Matrix33        *localRotation = NULL ;
        auto RealArray1D     *eigenValues   ;
        auto RealArray2D     *eigenVectors  ;
        auto Vector3         *center = NULL ;
        auto SymmetricMatrix *m             ;
        /* . Get the center of self. */
        center = Vector3_Allocate ( ) ;
        Coordinates3_Center ( self, selection, weights, &center ) ;
        /* . Translate self to its center. */
        Vector3_Scale ( center, -1.0e+00 ) ;
        Coordinates3_Translate ( self, center, NULL ) ;
        if ( translation != NULL ) Vector3_CopyTo ( center, translation, NULL ) ;
        /* . Determine the center of other. */
        Coordinates3_Center ( other, selection, weights, &center ) ;
        cx = Vector3_Item ( center, 0 ) ;
        cy = Vector3_Item ( center, 1 ) ;
        cz = Vector3_Item ( center, 2 ) ;
        /* . Determine the M matrix components. */
        /* . The four cases. */
        if ( selection == NULL )
        {
            if ( weights == NULL )
            {
                for ( i = 0, diagonal = 0.0e+00 ; i < self->extent0 ; i++ )
                {
                    Coordinates3_GetRow ( self,  i, x1, y1, z1 ) ;
                    Coordinates3_GetRow ( other, i, x2, y2, z2 ) ; x2 -= cx ; y2 -= cy ; z2 -= cz ;
                    diagonal += x1 * x1 + y1 * y1 + z1 * z1 + x2 * x2 + y2 * y2 + z2 * z2 ;
                    xx       += 2.0e+00 * x1 * x2 ;
                    xy       += 2.0e+00 * x1 * y2 ;
                    xz       += 2.0e+00 * x1 * z2 ;
                    yx       += 2.0e+00 * y1 * x2 ;
                    yy       += 2.0e+00 * y1 * y2 ;
                    yz       += 2.0e+00 * y1 * z2 ;
                    zx       += 2.0e+00 * z1 * x2 ;
                    zy       += 2.0e+00 * z1 * y2 ;
                    zz       += 2.0e+00 * z1 * z2 ;
                }
            }
            else
            {
                for ( i = 0, diagonal = 0.0e+00 ; i < self->extent0 ; i++ )
                {
                    Coordinates3_GetRow ( self,  i, x1, y1, z1 ) ;
                    Coordinates3_GetRow ( other, i, x2, y2, z2 ) ; x2 -= cx ; y2 -= cy ; z2 -= cz ;
                    w         = Array1D_Item ( weights, i ) ;
                    diagonal += w * ( x1 * x1 + y1 * y1 + z1 * z1 + x2 * x2 + y2 * y2 + z2 * z2 ) ;
                    xx       += 2.0e+00 * w * x1 * x2 ;
                    xy       += 2.0e+00 * w * x1 * y2 ;
                    xz       += 2.0e+00 * w * x1 * z2 ;
                    yx       += 2.0e+00 * w * y1 * x2 ;
                    yy       += 2.0e+00 * w * y1 * y2 ;
                    yz       += 2.0e+00 * w * y1 * z2 ;
                    zx       += 2.0e+00 * w * z1 * x2 ;
                    zy       += 2.0e+00 * w * z1 * y2 ;
                    zz       += 2.0e+00 * w * z1 * z2 ;
                }
            }
        }
        else
        {
            if ( weights == NULL )
            {
                for ( s = 0, diagonal = 0.0e+00 ; s < selection->capacity ; s++ )
                {
                    i = selection->indices[s] ;
                    Coordinates3_GetRow ( self,  i, x1, y1, z1 ) ;
                    Coordinates3_GetRow ( other, i, x2, y2, z2 ) ; x2 -= cx ; y2 -= cy ; z2 -= cz ;
                    diagonal += x1 * x1 + y1 * y1 + z1 * z1 + x2 * x2 + y2 * y2 + z2 * z2 ;
                    xx       += 2.0e+00 * x1 * x2 ;
                    xy       += 2.0e+00 * x1 * y2 ;
                    xz       += 2.0e+00 * x1 * z2 ;
                    yx       += 2.0e+00 * y1 * x2 ;
                    yy       += 2.0e+00 * y1 * y2 ;
                    yz       += 2.0e+00 * y1 * z2 ;
                    zx       += 2.0e+00 * z1 * x2 ;
                    zy       += 2.0e+00 * z1 * y2 ;
                    zz       += 2.0e+00 * z1 * z2 ;
                }
            }
            else
            {
                for ( s = 0, diagonal = 0.0e+00 ; s < selection->capacity ; s++ )
                {
                    i = selection->indices[s] ;
                    Coordinates3_GetRow ( self,  i, x1, y1, z1 ) ;
                    Coordinates3_GetRow ( other, i, x2, y2, z2 ) ; x2 -= cx ; y2 -= cy ; z2 -= cz ;
                    w         = Array1D_Item ( weights, i ) ;
                    diagonal += w * ( x1 * x1 + y1 * y1 + z1 * z1 + x2 * x2 + y2 * y2 + z2 * z2 ) ;
                    xx       += 2.0e+00 * w * x1 * x2 ;
                    xy       += 2.0e+00 * w * x1 * y2 ;
                    xz       += 2.0e+00 * w * x1 * z2 ;
                    yx       += 2.0e+00 * w * y1 * x2 ;
                    yy       += 2.0e+00 * w * y1 * y2 ;
                    yz       += 2.0e+00 * w * y1 * z2 ;
                    zx       += 2.0e+00 * w * z1 * x2 ;
                    zy       += 2.0e+00 * w * z1 * y2 ;
                    zz       += 2.0e+00 * w * z1 * z2 ;
                }
            }
        }
        /* . Fill the M matrix. */
        m = SymmetricMatrix_AllocateWithExtent ( 4, NULL ) ;
        SymmetricMatrix_Item ( m, 0, 0 ) = - xx - yy - zz + diagonal ;
        SymmetricMatrix_Item ( m, 1, 0 ) =   yz - zy ;
        SymmetricMatrix_Item ( m, 1, 1 ) = - xx + yy + zz + diagonal ;
        SymmetricMatrix_Item ( m, 2, 0 ) =   zx - xz ;
        SymmetricMatrix_Item ( m, 2, 1 ) = - xy - yx ;
        SymmetricMatrix_Item ( m, 2, 2 ) =   xx - yy + zz + diagonal ;
        SymmetricMatrix_Item ( m, 3, 0 ) =   xy - yx ;
        SymmetricMatrix_Item ( m, 3, 1 ) = - xz - zx ;
        SymmetricMatrix_Item ( m, 3, 2 ) = - yz - zy ;
        SymmetricMatrix_Item ( m, 3, 3 ) =   xx + yy - zz + diagonal ;
        /* . Find the eigenvalues and eigenvectors. */
        eigenValues  = RealArray1D_AllocateWithExtent ( 4   , NULL ) ;
        eigenVectors = RealArray2D_AllocateWithExtents ( 4, 4, NULL ) ;
        SymmetricMatrix_EigenvaluesSolve    ( m, False, 0, 4, eigenValues, eigenVectors, False, NULL ) ;
        /* . Determine the rotation from the eigenvector of smallest eigenvalue. */
        Matrix33_RotationFromQuaternion ( &localRotation, Array2D_Item ( eigenVectors, 0, 0 ),
                                                          Array2D_Item ( eigenVectors, 1, 0 ),
                                                          Array2D_Item ( eigenVectors, 2, 0 ),
                                                          Array2D_Item ( eigenVectors, 3, 0 ) ) ;
        /* . Transform the matrix (localRotation before translation). */
        Coordinates3_Rotate    ( self, localRotation, NULL ) ;
        Coordinates3_Translate ( self, center,        NULL ) ;
        if ( rotation    != NULL ) Matrix33_CopyTo ( localRotation, rotation, NULL ) ;
        if ( translation != NULL ) Vector3_Add ( translation, 1.0e+00, center, NULL ) ;
        /* . Finish up. */
        Matrix33_Deallocate        ( &localRotation ) ;
        RealArray1D_Deallocate     ( &eigenValues   ) ;
        RealArray2D_Deallocate     ( &eigenVectors  ) ;
        SymmetricMatrix_Deallocate ( &m             ) ;
        Vector3_Deallocate         ( &center        ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Perform a principal axis transformation.
! . |selection| determines the calculation of the transformation but not which rows are transformed.
! . Note that the transformation is not unique as it depends on the sign of the axes of inertia (the eigenvectors).
!---------------------------------------------------------------------------------------------------------------------------------*/
void Coordinates3_ToPrincipalAxes ( Coordinates3 *self, const Selection *selection, const RealArray1D *weights )
{
    if ( self != NULL )
    {
        auto Matrix33 *axes = NULL ;
        auto Vector3  *center = NULL, *moments = NULL ;
        /* . Translate. */
        Coordinates3_TranslateToCenter ( self, selection, weights ) ;
        /* . Get the rotation. */
        moments = Vector3_Allocate  ( ) ;
        axes    = Matrix33_Allocate ( ) ;
        Coordinates3_MomentsOfInertia ( self, selection, weights, moments, axes ) ;
        Matrix33_Transpose ( axes, NULL ) ;
        Coordinates3_Rotate ( self, axes, NULL ) ;
        /* . Finish up. */
        Matrix33_Deallocate ( &axes    ) ;
        Vector3_Deallocate  ( &center  ) ;
        Vector3_Deallocate  ( &moments ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Transformation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Coordinates3_Transform ( Coordinates3 *self, const Transformation3 *transformation3, const Selection *selection )
{
    if ( ( self != NULL ) && ( transformation3 != NULL ) )
    {
        Coordinates3_Rotate    ( self, transformation3->rotation,    selection ) ;
        Coordinates3_Translate ( self, transformation3->translation, selection ) ;
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Translation.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Coordinates3_Translate ( Coordinates3 *self, const Vector3 *translation, const Selection *selection )
{
    if ( ( self != NULL ) && ( translation != NULL ) )
    {
        auto Integer i, s ;
        auto Real    tx, ty, tz, x, y, z ;
        tx = Vector3_Item ( translation, 0 ) ;
        ty = Vector3_Item ( translation, 1 ) ;
        tz = Vector3_Item ( translation, 2 ) ;
        if ( selection == NULL )
        {
            for ( i = 0 ; i < self->extent0 ; i++ )
            {
                Coordinates3_GetRow ( self, i, x, y, z ) ;
                x += tx ; y += ty ; z += tz ;
                Coordinates3_SetRow ( self, i, x, y, z ) ;
            }
        }
        else
        {
            for ( s = 0 ; s < selection->capacity ; s++ )
            {
                i = selection->indices[s] ;
                Coordinates3_GetRow ( self, i, x, y, z ) ;
                x += tx ; y += ty ; z += tz ;
                Coordinates3_SetRow ( self, i, x, y, z ) ;
            }
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Translation to center.
! . |selection| determines the calculation of the center but not which
! . rows are translated.
!---------------------------------------------------------------------------------------------------------------------------------*/
void Coordinates3_TranslateToCenter ( Coordinates3 *self, const Selection *selection, const RealArray1D *weights )
{
    if ( self != NULL )
    {
        auto Vector3 *center = NULL ;
        /* . Center. */
        Coordinates3_Center ( self, selection, weights, &center ) ;
        Vector3_Scale ( center, -1.0e+00 ) ;
        /* . Translate. */
        Coordinates3_Translate ( self, center, NULL ) ;
        /* . Finish up. */
        Vector3_Deallocate ( &center ) ;
    }
}
