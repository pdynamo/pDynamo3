/*==================================================================================================================================
! . Marching cubes algorithm.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
!
! . Code for marching cubes adapted from the following:
!
! * @file    MarchingCubes.cpp
! * @author  Thomas Lewiner <thomas.lewiner@polytechnique.org>
! * @author  Math Dept, PUC-Rio
! * @version 0.2
! * @date    12/08/2002
! *
! * @brief   MarchingCubes Algorithm
!
! . Counterclockwise vertex order for the triangles.
!
!---------------------------------------------------------------------------------------------------------------------------------*/

# include <math.h>
# include <stdlib.h>
# include <stdio.h>

# include "Array2D_N3Macros.h"
# include "Boolean.h"
# include "Integer.h"
# include "IntegerArrayND.h"
# include "MarchingCubes.h"
# include "NumericalMacros.h"

/*----------------------------------------------------------------------------------------------------------------------------------
! . Parameters.
!---------------------------------------------------------------------------------------------------------------------------------*/
/* . A small value for face and interior testing. */
# define _Epsilon 1.0e-10

/* . Factor for estimating the starting number of polygons. */
# define _PolygonFactor0 4

/* . Polygon number increments. */
# define _PolygonIncrement 5000

/* . A small value to reduce numerical problems in linear interpolation. */
# define _SafeMinimum 1.0e-10

/* . Starting vertex count. */
# define _VertexCount0 10000

/* . Vertex number increments. */
# define _VertexIncrement1 10000
# define _VertexIncrement2  5000

/*----------------------------------------------------------------------------------------------------------------------------------
! . Local procedures.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void    AddTriangle         (       Integer        *nTriangles       ,
                                           IntegerArray2D *triangles        ,
                                     const IntegerArrayND *intersections    ,
                                     const Integer         i                ,
                                     const Integer         j                ,
                                     const Integer         k                ,
                                     const Integer8       *trig             ,
                                     const Integer8        n                ,
                                     const Integer         v12              ) ;
static Integer AddInteriorVertex   ( const Integer         i                ,
                                     const Integer         j                ,
                                     const Integer         k                ,
                                     const IntegerArrayND *intersections    ,
                                           Integer        *nVertices        ,
                                           RealArray2D    *vertices         ,
                                           RealArray2D    *normals          ) ;
static void    GetGradient         ( const RealArrayND    *data             ,
                                     const Integer         i                ,
                                     const Integer         j                ,
                                     const Integer         k                ,
                                     const Integer         nI               ,
                                     const Integer         nJ               ,
                                     const Integer         nK               ,
                                           Real           *gX               ,
                                           Real           *gY               ,
                                           Real           *gZ               ) ;
static Real    LinearlyInterpolate ( const Real            f0               ,
                                     const Real            f1               ,
                                     const Real            u0               ) ;
static void    PrintCube           ( const Real           *cube             ) ;
static void    ProcessCube         ( const Integer         i                ,
                                     const Integer         j                ,
                                     const Integer         k                ,
                                     const Real           *cube             ,
                                     const Integer8        cubeCase         ,
                                     const Integer8        configuration    ,
                                     const IntegerArrayND *intersections    ,
                                           Integer        *nVertices        ,
                                           RealArray2D    *vertices         ,
                                           RealArray2D    *normals          ,
                                           Integer        *nTriangles       ,
                                           IntegerArray2D *triangles        ) ;
static Boolean TestFace            ( const Real           *cube             ,
                                     const Integer8        face             ) ;
static Boolean TestInterior        ( const Real           *cube             ,
                                     const Integer8        cubeCase         ,
                                     const Integer8        configuration    ,
                                     const Integer8        subConfiguration ,
                                     const Integer8        s                ) ;

/*----------------------------------------------------------------------------------------------------------------------------------
! . Generate an isosurface given data on a regular 3-D grid.
! . Valid vertex normals, polygons and vertices arrays must be supplied on entry.
!---------------------------------------------------------------------------------------------------------------------------------*/
void MarchingCubes_Isosurface3D ( const RegularGrid    *grid           ,
                                  const RealArrayND    *data           ,
                                  const Real            isoValue       ,
                                        RealArray2D    *polygonNormals ,
                                        IntegerArray2D *polygons       ,
                                        RealArray2D    *vertexNormals  ,
                                        RealArray2D    *vertices       ,
                                        Status         *status         )
{
    if ( ( grid           != NULL ) &&
         ( data           != NULL ) &&
         ( polygonNormals != NULL ) &&
         ( polygons       != NULL ) &&
         ( vertexNormals  != NULL ) &&
         ( vertices       != NULL ) &&
         Status_IsOK ( status ) )
    {
        /* . Check that the grid has three dimensions and the data is compatible with the grid. */
        if ( ( grid->ndimensions == 3 ) && RegularGrid_IsConformingRealArrayND ( grid, data ) )
        {
            auto Cardinal8       tableEntry ;
            auto Integer         d, i, index, indices[4], j, k, nCubes, nI, nJ, nK, nold, nPolygons, nPolygons0, nVertices, nVertices0, p ;
            auto Integer8        cubeCase ;
            auto IntegerArrayND *intersections = NULL ;
            auto Real            cube[8], f0, faxis[3], gX, gY, gZ, hX, hY, hZ, u ;
            auto RealArray1D     view1D ;
            /* . Get the grid extents. */
            nI = data->view->extents[0] ;
            nJ = data->view->extents[1] ;
            nK = data->view->extents[2] ;
            /* . Allocate the intersections. */
            indices[0]    = nI ;
            indices[1]    = nJ ;
            indices[2]    = nK ;
            indices[3]    =  3 ;
            intersections = IntegerArrayND_AllocateWithShape ( 4, indices, status ) ;
            IntegerArrayND_Set ( intersections, -1, status ) ;
            /* . Ensure that there is enough space for the vertices - be conservative. */
            nCubes     = RegularGrid_NumberOfGridPoints ( grid ) ;
            nVertices0 = Minimum ( 3 * nCubes, _VertexCount0 ) ;
            RealArray2D_ResizeWithInitializer ( vertexNormals, nVertices0, 0.0e+00, status ) ;
            RealArray2D_ResizeWithInitializer ( vertices     , nVertices0, 0.0e+00, status ) ;
            if ( ! Status_IsOK ( status ) ) goto FinishUp ;
            /* . Compute intersections for each cube along the cube edges (almost all of them). */
            for ( i = nVertices = 0 ; i < nI ; i++ )
            {
                for ( j = 0 ; j < nJ ; j++ )
                {
                    for ( k = 0 ; k < nK ; k++ )
                    {
                        /* . Make sure that there is enough space for the maximum three vertices that can be added here. */
                        if ( nVertices + 3 > nVertices0 )
                        {
                            nVertices0 += _VertexIncrement1 ;
                            RealArray2D_ResizeWithInitializer ( vertexNormals, nVertices0, 0.0e+00, status ) ;
                            RealArray2D_ResizeWithInitializer ( vertices     , nVertices0, 0.0e+00, status ) ;
                            if ( ! Status_IsOK ( status ) ) goto FinishUp ;
                        }
                        /* . Get the function values along the lower corner of the cube. */
                        f0 = ArrayND_Item3D ( data, i, j, k ) - isoValue ;
                        if ( i < nI - 1 ) faxis[0] = ArrayND_Item3D ( data, i+1, j, k ) - isoValue ;
                        else              faxis[0] = f0 ;
                        if ( j < nJ - 1 ) faxis[1] = ArrayND_Item3D ( data, i, j+1, k ) - isoValue ;
                        else              faxis[1] = f0 ;
                        if ( k < nK - 1 ) faxis[2] = ArrayND_Item3D ( data, i, j ,k+1 ) - isoValue ;
                        else              faxis[2] = f0 ;
                        /* . Determine if there are any vertices. */
                        nold = nVertices ;
                        if ( f0 < 0 )
                        {
                            for ( d = 0 ; d < 3 ; d++ )
                            {
                                if ( faxis[d] >= 0 )
                                {
                                    Array2D_Item   ( vertices, nVertices, d )    = LinearlyInterpolate ( f0, faxis[d], 1.0e+00 ) ;
                                    ArrayND_Item4D ( intersections, i, j, k, d ) = nVertices ;
                                    nVertices ++ ;
                                }
                            }
                        }
                        else
                        {
                            for ( d = 0 ; d < 3 ; d++ )
                            {
                                if ( faxis[d] < 0 )
                                {
                                    Array2D_Item   ( vertices, nVertices, d )    = LinearlyInterpolate ( f0, faxis[d], 0.0e+00 ) ;
                                    ArrayND_Item4D ( intersections, i, j, k, d ) = nVertices ;
                                    nVertices ++ ;
                                }
                            }
                        }
                        /* . Process the vertex information. */
                        if ( nVertices > nold )
                        {
                            /* . Get the gradient at i, j, k. */
                            GetGradient ( data, i, j, k, nI, nJ, nK, &gX, &gY, &gZ ) ;
                            /* . Loop over the dimensions. */
                            for ( d = 0 ; d < 3 ; d++ )
                            {
                                index = ArrayND_Item4D ( intersections, i, j, k, d ) ;
                                if ( index >= 0 )
                                {
                                    u = Array2D_Item ( vertices, index, d ) ;
                                    Array2D_IncrementRowN3 ( vertices, index, ( Real ) i, ( Real ) j, ( Real ) k ) ;
                                    /* . Get the gradient along the appropriate axis. */
                                    GetGradient ( data, i + IJKTERMS[d][0], j + IJKTERMS[d][1], k + IJKTERMS[d][2], nI, nJ, nK, &hX, &hY, &hZ ) ;
                                    /* . Determine the normal. */
                                    Array2D_SetRowN3 ( vertexNormals, index, ( 1.0e+00 - u ) * gX + u * hX ,
                                                                             ( 1.0e+00 - u ) * gY + u * hY ,
                                                                             ( 1.0e+00 - u ) * gZ + u * hZ ) ;
                                }
                            }
                        }
                    }
                }
            }
            /* . Ensure that there is enough space for the polygons. */
            nPolygons0 = _PolygonFactor0 * nVertices ;
            IntegerArray2D_ResizeWithInitializer ( polygons, nPolygons0, -1, status ) ;
            if ( ! Status_IsOK ( status ) ) goto FinishUp ;
            /* . Process each cube. */
            for ( i = nPolygons = 0 ; i < nI - 1 ; i++ )
            {
                for ( j = 0 ; j < nJ - 1 ; j++ )
                {
                    for ( k = 0 ; k < nK - 1 ; k++ )
                    {
                        /* . Make sure that there is enough space for the maximum one vertex and twelve polygons that can be added here. */
                        if ( nPolygons + 12 > nPolygons0 )
                        {
                            nPolygons0 += _PolygonIncrement ;
                            IntegerArray2D_ResizeWithInitializer ( polygons, nPolygons0, -1, status ) ;
                            if ( ! Status_IsOK ( status ) ) goto FinishUp ;
                        }
                        if ( nVertices +  1 > nVertices0 )
                        {
                            nVertices0 += _VertexIncrement2 ;
                            RealArray2D_ResizeWithInitializer ( vertexNormals, nVertices0, 0.0e+00, status ) ;
                            RealArray2D_ResizeWithInitializer ( vertices     , nVertices0, 0.0e+00, status ) ;
                            if ( ! Status_IsOK ( status ) ) goto FinishUp ;
                        }
                        /* . Determine to which case the cube belongs (values between 0 and 255). */
                        tableEntry = 0 ;
                        for ( p = 0 ; p < 8 ; ++p )
                        {
                            cube[p] = ArrayND_Item3D ( data,  i+((p^(p>>1))&1), j+((p>>1)&1), k+((p>>2)&1) ) - isoValue ;
                            if ( cube[p] >= 0 ) tableEntry += 1 << p ;
                        }
                        /* . Process the cube if necessary. */
                        cubeCase = CUBECASES[tableEntry][0] ;
                        if ( cubeCase > 0 ) ProcessCube ( i, j, k                  ,
                                                          cube                     ,
                                                          cubeCase                 ,
                                                          CUBECASES[tableEntry][1] ,
                                                          intersections            ,
                                                          &nVertices               ,
                                                          vertices                 ,
                                                          vertexNormals            ,
                                                          &nPolygons               ,
                                                          polygons                 ) ;
                    }
                }
            }
            /* . Resize the surface data structure. */
            RealArray2D_Resize    ( polygonNormals, nPolygons, status ) ;
            IntegerArray2D_Resize ( polygons      , nPolygons, status ) ;
            RealArray2D_Resize    ( vertexNormals , nVertices, status ) ;
            RealArray2D_Resize    ( vertices      , nVertices, status ) ;
            if ( ! Status_IsOK ( status ) ) goto FinishUp ;
            /* . Scale the vertices and vertex normals to the grid coordinates. */
            for ( d = 0 ; d < 3 ; d++ )
            {
                RealArray2D_ColumnView ( vertexNormals, d, False, &view1D, NULL ) ;
                RealArray1D_Scale      ( &view1D      , grid->dimensions[d].binSize ) ;
                RealArray2D_ColumnView ( vertices     , d, False, &view1D, NULL ) ;
                RealArray1D_Scale      ( &view1D      , grid->dimensions[d].binSize ) ;
            }
            /* . Translate the vertex coordinates to the correct origin. */
            for ( d = 0 ; d < 3 ; d++ )
            {
                RealArray2D_ColumnView ( vertices, d, False, &view1D, NULL ) ;
                RealArray1D_Increment  ( &view1D, grid->dimensions[d].midPointLower ) ;
            }
            /* . Normalize the vertex normals. */
            for ( i = 0 ; i < View2D_Rows ( vertexNormals ) ; i++ )
            {
                RealArray2D_RowView   ( vertexNormals, i, False, &view1D, NULL ) ;
                RealArray1D_Normalize ( &view1D, NULL, NULL ) ;
            }
            /* . Initialize the polygon normals. */
            RealArray2D_Set ( polygonNormals, 0.0e+00 ) ;
            /* . Finish up. */
        FinishUp:
            IntegerArrayND_Deallocate ( &intersections ) ;
        }
        /* . Argument error. */
        else Status_Set ( status, Status_InvalidArgument ) ;
    }
}

/*==================================================================================================================================
! . Marching cubes auxiliary procedures.
!=================================================================================================================================*/
/*----------------------------------------------------------------------------------------------------------------------------------
! . Add triangles.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void AddTriangle (       Integer        *nTriangles    ,
                                IntegerArray2D *triangles     ,
                          const IntegerArrayND *intersections ,
                          const Integer         i             ,
                          const Integer         j             ,
                          const Integer         k             ,
                          const Integer8       *trig          ,
                          const Integer8        n             ,
                          const Integer         v12           )
{
    Integer t, tv[3] ;
    for( t = 0 ; t < 3*n ; t++ )
    {
        switch ( trig[t] )
        {
            case  0 : tv[ t % 3 ] = ArrayND_Item4D ( intersections, i  , j  , k  , 0 ) ; break ;
            case  1 : tv[ t % 3 ] = ArrayND_Item4D ( intersections, i+1, j  , k  , 1 ) ; break ;
            case  2 : tv[ t % 3 ] = ArrayND_Item4D ( intersections, i  , j+1, k  , 0 ) ; break ;
            case  3 : tv[ t % 3 ] = ArrayND_Item4D ( intersections, i  , j  , k  , 1 ) ; break ;
            case  4 : tv[ t % 3 ] = ArrayND_Item4D ( intersections, i  , j  , k+1, 0 ) ; break ;
            case  5 : tv[ t % 3 ] = ArrayND_Item4D ( intersections, i+1, j  , k+1, 1 ) ; break ;
            case  6 : tv[ t % 3 ] = ArrayND_Item4D ( intersections, i  , j+1, k+1, 0 ) ; break ;
            case  7 : tv[ t % 3 ] = ArrayND_Item4D ( intersections, i  , j  , k+1, 1 ) ; break ;
            case  8 : tv[ t % 3 ] = ArrayND_Item4D ( intersections, i  , j  , k  , 2 ) ; break ;
            case  9 : tv[ t % 3 ] = ArrayND_Item4D ( intersections, i+1, j  , k  , 2 ) ; break ;
            case 10 : tv[ t % 3 ] = ArrayND_Item4D ( intersections, i+1, j+1, k  , 2 ) ; break ;
            case 11 : tv[ t % 3 ] = ArrayND_Item4D ( intersections, i  , j+1, k  , 2 ) ; break ;
            case 12 : tv[ t % 3 ] = v12 ; break ;
            default : break ;
        }
        if ( tv[t%3] == -1 ) { printf ( "Marching Cubes: invalid triangle %d %d %d %d\n", (*nTriangles) + 1, i, j, k ) ; }
        if ( t%3 == 2 )
        {
            Array2D_SetRowN3 ( triangles, (*nTriangles), tv[0], tv[1], tv[2] ) ;
            (*nTriangles)++ ;
        }
    }
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Add an interior vertex using the average of the intersection points of a cube.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Integer AddInteriorVertex ( const Integer         i             ,
                                   const Integer         j             ,
                                   const Integer         k             ,
                                   const IntegerArrayND *intersections ,
                                         Integer        *nVertices     ,
                                         RealArray2D    *vertices      ,
                                         RealArray2D    *normals       )
{
    Integer c, current, di, dj, dk, n, t, vid ;
    Real    scale, x, y, z ;
    current = (*nVertices) ;
    for ( c = n = 0 ; c < 3 ; c++ )
    {
        for ( t = 0 ; t < 4 ; t++ )
        {
            di  = i + IVERTEXTERMS[c][t][0] ;
            dj  = j + IVERTEXTERMS[c][t][1] ;
            dk  = k + IVERTEXTERMS[c][t][2] ;
            vid = ArrayND_Item4D ( intersections, di, dj, dk, c ) ;
            if ( vid != -1 )
            {
                Array2D_GetRowN3       ( vertices, vid,     x, y, z ) ;
                Array2D_IncrementRowN3 ( vertices, current, x, y, z ) ;
                Array2D_GetRowN3       ( normals,  vid,     x, y, z ) ;
                Array2D_IncrementRowN3 ( normals,  current, x, y, z ) ;
                n++ ;
            }
        }
    }
    scale = 1.0e+00 / ( Real ) n ;
    Array2D_ScaleRowN3 ( vertices, current, scale ) ;
    Array2D_ScaleRowN3 ( vertices, current, scale ) ;
    (*nVertices)++ ;
    return current ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Gradient calculation at a grid point using finite differences.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void GetGradient ( const RealArrayND *data ,
                          const Integer      i    ,
                          const Integer      j    ,
                          const Integer      k    ,
                          const Integer      nI   ,
                          const Integer      nJ   ,
                          const Integer      nK   ,
                                Real        *gX   ,
                                Real        *gY   ,
                                Real        *gZ   )
{
    if ( i > 0 )
    {
        if ( i < nI - 1 ) (*gX) = 0.5e+00 * ( ArrayND_Item3D ( data, i+1, j, k ) - ArrayND_Item3D ( data, i-1, j, k ) ) ;
        else              (*gX) =             ArrayND_Item3D ( data, i  , j, k ) - ArrayND_Item3D ( data, i-1, j, k )   ;
    }
    else                  (*gX) =             ArrayND_Item3D ( data, i+1, j, k ) - ArrayND_Item3D ( data, i  , j, k )   ;
    if ( j > 0 )
    {
        if ( j < nJ - 1 ) (*gY) = 0.5e+00 * ( ArrayND_Item3D ( data, i, j+1, k ) - ArrayND_Item3D ( data, i, j-1, k ) ) ;
        else              (*gY) =             ArrayND_Item3D ( data, i, j  , k ) - ArrayND_Item3D ( data, i, j-1, k )   ;
    }
    else                  (*gY) =             ArrayND_Item3D ( data, i, j+1, k ) - ArrayND_Item3D ( data, i, j  , k )   ;
    if ( k > 0 )
    {
        if ( k < nK - 1 ) (*gZ) = 0.5e+00 * ( ArrayND_Item3D ( data, i, j, k+1 ) - ArrayND_Item3D ( data, i, j, k-1 ) ) ;
        else              (*gZ) =             ArrayND_Item3D ( data, i, j, k   ) - ArrayND_Item3D ( data, i, j, k-1 )   ;
    }
    else                  (*gZ) =             ArrayND_Item3D ( data, i, j, k+1 ) - ArrayND_Item3D ( data, i, j, k   )   ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Linear interpolation.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Real LinearlyInterpolate ( const Real f0, const Real f1, const Real u0 )
{
    Real delta, u = u0 ;
    delta = f0 - f1 ;
    if ( fabs ( delta ) > _SafeMinimum ) u = f0 / delta ;
    return u ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Print a cube for debugging.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void PrintCube ( const Real *cube ) { printf ( "\t%f %f %f %f %f %f %f %f\n", cube[0], cube[1], cube[2], cube[3], cube[4], cube[5], cube[6], cube[7] ) ; }

/*----------------------------------------------------------------------------------------------------------------------------------
! . Tesselate a cube.
! . This procedure adds at most one vertex and at most 12 triangles. As no checks are made on storage here, it should be ensured
! . that the vertex and triangle arrays have at least this amount of space available on entry.
!---------------------------------------------------------------------------------------------------------------------------------*/
static void ProcessCube ( const Integer         i             ,
                          const Integer         j             ,
                          const Integer         k             ,
                          const Real           *cube          ,
                          const Integer8        cubeCase      ,
                          const Integer8        configuration ,
                          const IntegerArrayND *intersections ,
                                Integer        *nVertices     ,
                                RealArray2D    *vertices      ,
                                RealArray2D    *normals       ,
                                Integer        *nTriangles    ,
                                IntegerArray2D *triangles     )
{
    Integer8 subConfiguration = 0 ;
    Integer  v12 = -1 ;
    switch ( cubeCase )
    {
/* . Done earlier.
    case  0 :
        break ;
*/
    case  1 :
        AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING1[configuration], 1, -1 ) ; break ;

    case  2 :
        AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING2[configuration], 2, -1 ) ; break ;

    case  3 :
        if ( TestFace ( cube, TEST3[configuration] ) ) AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING3_2[configuration], 4, -1 ) ; /* .  3.2 */
        else AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING3_1[configuration], 2, -1 ) ; /* .  3.1 */
        break ;

    case  4 :
        if ( TestInterior ( cube, cubeCase, configuration, subConfiguration, TEST4[configuration] ) ) AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING4_1[configuration], 2, -1 ) ; /* . 4.1.1 */
        else AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING4_2[configuration], 6, -1 ) ; /* .  4.1.2 */
        break ;

    case  5 :
        AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING5[configuration], 3, -1 ) ;
        break ;

    case  6 :
        if ( TestFace ( cube, TEST6[configuration][0] ) )
            AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING6_2[configuration], 5, -1 ) ; /* .  6.2 */
        else
        {
            if ( TestInterior ( cube, cubeCase, configuration, subConfiguration, TEST6[configuration][1] ) ) AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING6_1_1[configuration], 3, -1 ) ; /* .  6.1.1 */
            else AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING6_1_2[configuration], 7, -1 ) ; /* .  6.1.2 */
        }
        break ;

    case  7 :
        if ( TestFace ( cube, TEST7[configuration][0] ) ) subConfiguration +=  1 ;
        if ( TestFace ( cube, TEST7[configuration][1] ) ) subConfiguration +=  2 ;
        if ( TestFace ( cube, TEST7[configuration][2] ) ) subConfiguration +=  4 ;
        switch ( subConfiguration )
        {
            case 0 :
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING7_1[configuration], 3, -1 ) ; break ;
            case 1 :
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING7_2[configuration][0], 5, -1 ) ; break ;
            case 2 :
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING7_2[configuration][1], 5, -1 ) ; break ;
            case 3 :
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING7_3[configuration][0], 9, v12 ) ; break ;
            case 4 :
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING7_2[configuration][2], 5, -1 ) ; break ;
            case 5 :
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING7_3[configuration][1], 9, v12 ) ; break ;
            case 6 :
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING7_3[configuration][2], 9, v12 ) ; break ;
            case 7 :
                if ( TestInterior ( cube, cubeCase, configuration, subConfiguration, TEST7[configuration][3] ) ) AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING7_4_2[configuration], 9, -1 ) ;
                else AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING7_4_1[configuration], 5, -1 ) ;
              break ;
        } ;
        break ;

    case  8 :
        AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING8[configuration], 2, -1 ) ;
        break ;

    case  9 :
        AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING9[configuration], 4, -1 ) ;
        break ;

    case 10 :
        if ( TestFace ( cube, TEST10[configuration][0] ) )
        {
            if ( TestFace ( cube, TEST10[configuration][1] ) ) AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING10_1_1_[configuration], 4, -1 ) ; /* .  10.1.1 */
            else
            {
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING10_2[configuration], 8, v12 ) ; /* .  10.2 */
            }
        }
        else
        {
            if ( TestFace ( cube, TEST10[configuration][1] ) )
            {
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING10_2_[configuration], 8, v12 ) ; /* .  10.2 */
            }
            else
            {
                if ( TestInterior ( cube, cubeCase, configuration, subConfiguration, TEST10[configuration][2] ) ) AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING10_1_1[configuration], 4, -1 ) ; /* .  10.1.1 */
                else AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING10_1_2[configuration], 8, -1 ) ; /* .  10.1.2 */
            }
        }
        break ;

    case 11 :
        AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING11[configuration], 4, -1 ) ;
        break ;

    case 12 :
        if ( TestFace ( cube, TEST12[configuration][0] ) )
        {
            if ( TestFace ( cube, TEST12[configuration][1] ) ) AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING12_1_1_[configuration], 4, -1 ) ; /* .  12.1.1 */
            else
            {
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING12_2[configuration], 8, v12 ) ; /* .  12.2 */
            }
        }
        else
        {
            if ( TestFace ( cube, TEST12[configuration][1] ) )
            {
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING12_2_[configuration], 8, v12 ) ; /* .  12.2 */
            }
            else
            {
                if ( TestInterior ( cube, cubeCase, configuration, subConfiguration, TEST12[configuration][2] ) ) AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING12_1_1[configuration], 4, -1 ) ; /* .  12.1.1 */
                else AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING12_1_2[configuration], 8, -1 ) ; /* .  12.1.2 */
            }
        }
        break ;

    case 13 :
        if ( TestFace ( cube, TEST13[configuration][0] ) ) subConfiguration +=  1 ;
        if ( TestFace ( cube, TEST13[configuration][1] ) ) subConfiguration +=  2 ;
        if ( TestFace ( cube, TEST13[configuration][2] ) ) subConfiguration +=  4 ;
        if ( TestFace ( cube, TEST13[configuration][3] ) ) subConfiguration +=  8 ;
        if ( TestFace ( cube, TEST13[configuration][4] ) ) subConfiguration += 16 ;
        if ( TestFace ( cube, TEST13[configuration][5] ) ) subConfiguration += 32 ;
        switch ( SUBCONFIGURATION13[subConfiguration] )
        {
            case 0 :/* 13.1 */
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_1[configuration], 4, -1 ) ; break ;

            case 1 :/* 13.2 */
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_2[configuration][0], 6, -1 ) ; break ;
            case 2 :/* 13.2 */
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_2[configuration][1], 6, -1 ) ; break ;
            case 3 :/* 13.2 */
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_2[configuration][2], 6, -1 ) ; break ;
            case 4 :/* 13.2 */
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_2[configuration][3], 6, -1 ) ; break ;
            case 5 :/* 13.2 */
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_2[configuration][4], 6, -1 ) ; break ;
            case 6 :/* 13.2 */
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_2[configuration][5], 6, -1 ) ; break ;

            case 7 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_3[configuration][0], 10, v12 ) ; break ;
            case 8 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_3[configuration][1], 10, v12 ) ; break ;
            case 9 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_3[configuration][2], 10, v12 ) ; break ;
            case 10 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_3[configuration][3], 10, v12 ) ; break ;
            case 11 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_3[configuration][4], 10, v12 ) ; break ;
            case 12 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_3[configuration][5], 10, v12 ) ; break ;
            case 13 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_3[configuration][6], 10, v12 ) ; break ;
            case 14 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_3[configuration][7], 10, v12 ) ; break ;
            case 15 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_3[configuration][8], 10, v12 ) ; break ;
            case 16 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_3[configuration][9], 10, v12 ) ; break ;
            case 17 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_3[configuration][10], 10, v12 ) ; break ;
            case 18 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_3[configuration][11], 10, v12 ) ; break ;

            case 19 :/* 13.4 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_4[configuration][0], 12, v12 ) ; break ;
            case 20 :/* 13.4 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_4[configuration][1], 12, v12 ) ; break ;
            case 21 :/* 13.4 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_4[configuration][2], 12, v12 ) ; break ;
            case 22 :/* 13.4 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_4[configuration][3], 12, v12 ) ; break ;

            case 23 :/* 13.5 */
                subConfiguration = 0 ;
                if ( TestInterior ( cube, cubeCase, configuration, subConfiguration, TEST13[configuration][6] ) ) AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_5_1[configuration][0], 6, -1 ) ;
                else AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_5_2[configuration][0], 10, -1 ) ;
                break ;
            case 24 :/* 13.5 */
                subConfiguration = 1 ;
                if ( TestInterior ( cube, cubeCase, configuration, subConfiguration, TEST13[configuration][6] ) ) AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_5_1[configuration][1], 6, -1 ) ;
                else AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_5_2[configuration][1], 10, -1 ) ;
                break ;
            case 25 :/* 13.5 */
                subConfiguration = 2 ;
                if ( TestInterior ( cube, cubeCase, configuration, subConfiguration, TEST13[configuration][6] ) ) AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_5_1[configuration][2], 6, -1 ) ;
                else AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_5_2[configuration][2], 10, -1 ) ;
                break ;
            case 26 :/* 13.5 */
                subConfiguration = 3 ;
                if ( TestInterior ( cube, cubeCase, configuration, subConfiguration, TEST13[configuration][6] ) ) AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_5_1[configuration][3], 6, -1 ) ;
                else AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_5_2[configuration][3], 10, -1 ) ;
                break ;

            case 27 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_3_[configuration][0], 10, v12 ) ; break ;
            case 28 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_3_[configuration][1], 10, v12 ) ; break ;
            case 29 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_3_[configuration][2], 10, v12 ) ; break ;
            case 30 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_3_[configuration][3], 10, v12 ) ; break ;
            case 31 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_3_[configuration][4], 10, v12 ) ; break ;
            case 32 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_3_[configuration][5], 10, v12 ) ; break ;
            case 33 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_3_[configuration][6], 10, v12 ) ; break ;
            case 34 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_3_[configuration][7], 10, v12 ) ; break ;
            case 35 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_3_[configuration][8], 10, v12 ) ; break ;
            case 36 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_3_[configuration][9], 10, v12 ) ; break ;
            case 37 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_3_[configuration][10], 10, v12 ) ; break ;
            case 38 :/* 13.3 */
                v12 = AddInteriorVertex ( i, j, k, intersections, nVertices, vertices, normals ) ;
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_3_[configuration][11], 10, v12 ) ; break ;

            case 39 :/* 13.2 */
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_2_[configuration][0], 6, -1 ) ; break ;
            case 40 :/* 13.2 */
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_2_[configuration][1], 6, -1 ) ; break ;
            case 41 :/* 13.2 */
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_2_[configuration][2], 6, -1 ) ; break ;
            case 42 :/* 13.2 */
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_2_[configuration][3], 6, -1 ) ; break ;
            case 43 :/* 13.2 */
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_2_[configuration][4], 6, -1 ) ; break ;
            case 44 :/* 13.2 */
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_2_[configuration][5], 6, -1 ) ; break ;

            case 45 :/* 13.1 */
                AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING13_1_[configuration], 4, -1 ) ; break ;

            default :
                printf ( "Marching Cubes: Impossible case 13?\n" ) ;  PrintCube ( cube ) ;
        }
        break ;
    case 14 :
        AddTriangle ( nTriangles, triangles, intersections, i, j, k, TILING14[configuration], 4, -1 ) ;
        break ;
    } ;
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Tests if the s of the tesselation of the cube should be connected by the interior of an ambiguous face.
! . Returns True if the face contains a part of the surface.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Boolean TestFace ( const Real *cube, const Integer8 face )
{
    Real A,B,C,D ;
    switch ( face )
    {
        case -1 : case 1 :  A = cube[0] ;  B = cube[4] ;  C = cube[5] ;  D = cube[1] ;  break ;
        case -2 : case 2 :  A = cube[1] ;  B = cube[5] ;  C = cube[6] ;  D = cube[2] ;  break ;
        case -3 : case 3 :  A = cube[2] ;  B = cube[6] ;  C = cube[7] ;  D = cube[3] ;  break ;
        case -4 : case 4 :  A = cube[3] ;  B = cube[7] ;  C = cube[4] ;  D = cube[0] ;  break ;
        case -5 : case 5 :  A = cube[0] ;  B = cube[3] ;  C = cube[2] ;  D = cube[1] ;  break ;
        case -6 : case 6 :  A = cube[4] ;  B = cube[7] ;  C = cube[6] ;  D = cube[5] ;  break ;
        default : printf ( "Invalid face code %d\n", face ) ;  PrintCube ( cube ) ;  A = B = C = D = 0 ;
    } ;
    if ( fabs ( A*C - B*D ) < _Epsilon ) return face >= 0 ;
    return face * A * ( A*C - B*D ) >= 0  ;  /* . face and A invert signs. */
}

/*----------------------------------------------------------------------------------------------------------------------------------
! . Tests if the components of the tesselation of the cube should be connected through the interior of the cube.
! . If the interior is empty returns True for s = 7 and False for s = -7.
!---------------------------------------------------------------------------------------------------------------------------------*/
static Boolean TestInterior ( const Real     *cube             ,
                              const Integer8  cubeCase         ,
                              const Integer8  configuration    ,
                              const Integer8  subConfiguration ,
                              const Integer8  s                )
{
    Integer8 edge = -1, test = 0 ; /* . edge is the reference edge of the triangulation. */
    Real     At = 0, Bt = 0, Ct = 0, Dt = 0 ;
    Real     a, b, t ;
    switch ( cubeCase )
    {
        case  4 :
        case 10 :
            a = ( cube[4] - cube[0] ) * ( cube[6] - cube[2] ) - ( cube[7] - cube[3] ) * ( cube[5] - cube[1] ) ;
            b =  cube[2] * ( cube[4] - cube[0] ) + cube[0] * ( cube[6] - cube[2] ) - cube[1] * ( cube[7] - cube[3] ) - cube[3] * ( cube[5] - cube[1] ) ;
            t = - b / ( 2 * a ) ;
            if ( ( t < 0 ) || ( t > 1 ) ) return s > 0 ;
            At = cube[0] + ( cube[4] - cube[0] ) * t ;
            Bt = cube[3] + ( cube[7] - cube[3] ) * t ;
            Ct = cube[2] + ( cube[6] - cube[2] ) * t ;
            Dt = cube[1] + ( cube[5] - cube[1] ) * t ;
            break ;
        case  6 :
        case  7 :
        case 12 :
        case 13 :
            switch ( cubeCase )
            {
                case  6 : edge = TEST6 [configuration][2] ; break ;
                case  7 : edge = TEST7 [configuration][4] ; break ;
                case 12 : edge = TEST12[configuration][3] ; break ;
                case 13 : edge = TILING13_5_1[configuration][subConfiguration][0] ; break ;
            }
            switch ( edge )
            {
            case  0 :
                t  = cube[0] / ( cube[0] - cube[1] ) ;
                At = 0 ;
                Bt = cube[3] + ( cube[2] - cube[3] ) * t ;
                Ct = cube[7] + ( cube[6] - cube[7] ) * t ;
                Dt = cube[4] + ( cube[5] - cube[4] ) * t ;
                break ;
            case  1 :
                t  = cube[1] / ( cube[1] - cube[2] ) ;
                At = 0 ;
                Bt = cube[0] + ( cube[3] - cube[0] ) * t ;
                Ct = cube[4] + ( cube[7] - cube[4] ) * t ;
                Dt = cube[5] + ( cube[6] - cube[5] ) * t ;
                break ;
            case  2 :
                t  = cube[2] / ( cube[2] - cube[3] ) ;
                At = 0 ;
                Bt = cube[1] + ( cube[0] - cube[1] ) * t ;
                Ct = cube[5] + ( cube[4] - cube[5] ) * t ;
                Dt = cube[6] + ( cube[7] - cube[6] ) * t ;
                break ;
            case  3 :
                t  = cube[3] / ( cube[3] - cube[0] ) ;
                At = 0 ;
                Bt = cube[2] + ( cube[1] - cube[2] ) * t ;
                Ct = cube[6] + ( cube[5] - cube[6] ) * t ;
                Dt = cube[7] + ( cube[4] - cube[7] ) * t ;
                break ;
            case  4 :
                t  = cube[4] / ( cube[4] - cube[5] ) ;
                At = 0 ;
                Bt = cube[7] + ( cube[6] - cube[7] ) * t ;
                Ct = cube[3] + ( cube[2] - cube[3] ) * t ;
                Dt = cube[0] + ( cube[1] - cube[0] ) * t ;
                break ;
            case  5 :
                t  = cube[5] / ( cube[5] - cube[6] ) ;
                At = 0 ;
                Bt = cube[4] + ( cube[7] - cube[4] ) * t ;
                Ct = cube[0] + ( cube[3] - cube[0] ) * t ;
                Dt = cube[1] + ( cube[2] - cube[1] ) * t ;
                break ;
            case  6 :
                t  = cube[6] / ( cube[6] - cube[7] ) ;
                At = 0 ;
                Bt = cube[5] + ( cube[4] - cube[5] ) * t ;
                Ct = cube[1] + ( cube[0] - cube[1] ) * t ;
                Dt = cube[2] + ( cube[3] - cube[2] ) * t ;
                break ;
            case  7 :
                t  = cube[7] / ( cube[7] - cube[4] ) ;
                At = 0 ;
                Bt = cube[6] + ( cube[5] - cube[6] ) * t ;
                Ct = cube[2] + ( cube[1] - cube[2] ) * t ;
                Dt = cube[3] + ( cube[0] - cube[3] ) * t ;
                break ;
            case  8 :
                t  = cube[0] / ( cube[0] - cube[4] ) ;
                At = 0 ;
                Bt = cube[3] + ( cube[7] - cube[3] ) * t ;
                Ct = cube[2] + ( cube[6] - cube[2] ) * t ;
                Dt = cube[1] + ( cube[5] - cube[1] ) * t ;
                break ;
            case  9 :
                t  = cube[1] / ( cube[1] - cube[5] ) ;
                At = 0 ;
                Bt = cube[0] + ( cube[4] - cube[0] ) * t ;
                Ct = cube[3] + ( cube[7] - cube[3] ) * t ;
                Dt = cube[2] + ( cube[6] - cube[2] ) * t ;
                break ;
            case 10 :
                t  = cube[2] / ( cube[2] - cube[6] ) ;
                At = 0 ;
                Bt = cube[1] + ( cube[5] - cube[1] ) * t ;
                Ct = cube[0] + ( cube[4] - cube[0] ) * t ;
                Dt = cube[3] + ( cube[7] - cube[3] ) * t ;
                break ;
            case 11 :
                t  = cube[3] / ( cube[3] - cube[7] ) ;
                At = 0 ;
                Bt = cube[2] + ( cube[6] - cube[2] ) * t ;
                Ct = cube[1] + ( cube[5] - cube[1] ) * t ;
                Dt = cube[0] + ( cube[4] - cube[0] ) * t ;
                break ;
            default : printf ( "Invalid edge %d\n", edge ) ;  PrintCube ( cube ) ;  break ;
        }
        break ;
        default : printf ( "Invalid ambiguous case %d\n", cubeCase ) ;  PrintCube ( cube ) ;  break ;
    }
    if ( At >= 0 ) test ++   ;
    if ( Bt >= 0 ) test += 2 ;
    if ( Ct >= 0 ) test += 4 ;
    if ( Dt >= 0 ) test += 8 ;
    switch ( test )
    {
        case  0 : return s > 0 ;
        case  1 : return s > 0 ;
        case  2 : return s > 0 ;
        case  3 : return s > 0 ;
        case  4 : return s > 0 ;
        case  5 : if ( ( At * Ct - Bt * Dt ) <  _Epsilon ) return s > 0 ; break ;
        case  6 : return s > 0 ;
        case  7 : return s < 0 ;
        case  8 : return s > 0 ;
        case  9 : return s > 0 ;
        case 10 : if ( ( At * Ct - Bt * Dt ) >= _Epsilon ) return s > 0 ; break ;
        case 11 : return s < 0 ;
        case 12 : return s > 0 ;
        case 13 : return s < 0 ;
        case 14 : return s < 0 ;
        case 15 : return s < 0 ;
    }
    return ( s < 0 ) ;
}
