"""Marching cubes algorithm for generating isosurfaces in 3-D."""

from .PolygonalSurface import PolygonalSurface
from .SurfacesError    import SurfacesError

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def MarchingCubes_Isosurface3D ( RegularGrid grid not None ,
                                 RealArrayND data not None ,
                                 CReal       isoValue      ):
    """Construct an isosurface from regular 3-D grid data using the marching cubes algorithm."""
    cdef IntegerArray2D polygons
    cdef RealArray2D    polygonNormals
    cdef RealArray2D    vertexNormals
    cdef RealArray2D    vertices
    cdef CStatus        cStatus = CStatus_OK
    # . Guess sizes for initial allocation - resized as necessary within the marching cubes function.
    nVertices = min ( 3 * grid.size, 10000 )
    nPolygons = 4 * nVertices
    ( polygonNormals, polygons, vertexNormals, vertices ) = PolygonalSurface.AllocateArrays ( 3, nPolygons, nVertices )
    CMarchingCubes_Isosurface3D ( grid.cObject           ,
                                  data.cObject           ,
                                  isoValue               ,
                                  polygonNormals.cObject ,
                                  polygons.cObject       ,
                                  vertexNormals.cObject  ,
                                  vertices.cObject       ,
                                  &cStatus               )
    if cStatus != CStatus_OK: raise SurfacesError ( "Error generating isosurface." )
    surface = PolygonalSurface.FromArrays ( polygonNormals, polygons, vertexNormals, vertices )
    surface.MakePolygonNormals ( )
    return surface

