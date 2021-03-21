"""An object for polygonal surfaces."""

from   pCore         import DataType           , \
                            logFile            , \
                            LogFileActive      , \
                            SummarizableObject
from ..Arrays        import Array
from ..Geometry3     import Vector3
from  .SurfacesError import SurfacesError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PolygonalSurface ( SummarizableObject ):
    """An object for polygonal surfaces."""

    # . No extra summarizable.
    _attributable = dict ( SummarizableObject._attributable )
    _classLabel   = "Polygonal Surface"
    _attributable.update ( { "polygonNormals" : None ,
                             "polygons"       : None ,
                             "vertexNormals"  : None ,
                             "vertices"       : None } )

    @staticmethod
    def AllocateArrays ( rank, polygons, vertices ):
        """Allocate arrays with given sizes."""
        polygonNormals = Array.WithShape ( [ polygons, rank ]                              ) ; polygonNormals.Set  ( 0.0 )
        polygons       = Array.WithShape ( [ polygons, rank ], dataType = DataType.Integer ) ; polygons.Set        ( -1  )
        vertexNormals  = Array.WithShape ( [ vertices, rank ]                              ) ; vertexNormals.Set   ( 0.0 )
        vertices       = Array.WithShape ( [ vertices, rank ]                              ) ; vertices.Set        ( 0.0 )
        return ( polygonNormals, polygons, vertexNormals, vertices )

    def BuildFromArrays ( self, polygonNormals, polygons, vertexNormals, vertices ):
        """Fill a raw object from a set of arrays."""
        # . The normals can be None.
        isOK = ( ( polygons is not None ) and ( vertices is not None ) and ( polygons.shape[1] == vertices.shape[1] ) )
        if polygonNormals is not None: isOK = isOK and ( polygonNormals.shape == polygons.shape )
        if vertexNormals  is not None: isOK = isOK and ( vertexNormals.shape  == vertices.shape )
        if not isOK: raise SurfacesError ( "Constructing surface from incompatible arrays." )
        self.polygonNormals = polygonNormals
        self.polygons       = polygons
        self.vertexNormals  = vertexNormals
        self.vertices       = vertices

    @classmethod
    def FromArrays ( selfClass, polygonNormals, polygons, vertexNormals, vertices ):
        """Constructor given a set of arrays."""
        self = selfClass ( )
        self.BuildFromArrays ( polygonNormals, polygons, vertexNormals, vertices )
        return self

    def MakePolygonNormals ( self ):
        """Make polygon normals."""
        if self.rank == 3:
            a = Vector3.Null ( )
            b = Vector3.Null ( )
#            from ..Arrays import ArrayPrint
#            ArrayPrint ( self.polygonNormals[-1,:] )
#            print ( self.polygons.shape, self.polygonNormals.shape )
            for p in range ( self.polygons.rows ):
                self.vertices[self.polygons[p,1],:].CopyTo ( a ) ; a.Add ( self.vertices[self.polygons[p,0],:], scale = -1.0 )
                self.vertices[self.polygons[p,2],:].CopyTo ( b ) ; b.Add ( self.vertices[self.polygons[p,0],:], scale = -1.0 )
                a.Cross ( b )
                a.Normalize ( )
 #               ArrayPrint ( a )
#                print ( p, a.Norm2 ( ) )
#                try:
                a.CopyTo ( self.polygonNormals[p,:] )
#                except:
#                    print ( "\nCopying error.\n" )
        else:
            raise SurfacesError ( "Unable to make polygon normals for surfaces with rank other than 3." )

    def MakeVertexNormalsFromPolygonalNormals ( self ):
        """Make vertex normals from polygon normals."""
        self.vertexNormals.Set ( 0.0 )
        for p in range ( self.polygons.rows ):
            n = self.polygonNormals[p,:]
            for r in range ( self.rank ):
                self.vertexNormals[self.polygons[p,r],:].Add ( n )
        for v in range ( self.vertexNormals.rows ):
            self.vertexNormals[v,:].Normalize ( )

    def OriginAndExtents ( self ):
        """Get the origin and extents of the surface."""
        extents = Array.WithExtent ( self.rank )
        origin  = Array.WithExtent ( self.rank )
        for d in range ( self.rank ):
            column = self.vertices[:,d]
            lower  = column.Minimum ( )
            upper  = column.Maximum ( )
            extents[d] = upper - lower
            origin [d] = lower
        return ( origin, extents )

    def SummaryItems ( self ):
        """Summary items."""
        items = super ( PolygonalSurface, self ).SummaryItems ( )
        items.extend ( [ ( "Polygons", "{:d}".format ( self.polygons.rows ) ) ,
                         ( "Rank"    , "{:d}".format ( self.rank          ) ) ,
                         ( "Vertices", "{:d}".format ( self.vertices.rows ) ) ] )
        return items

    @classmethod
    def WithSizes ( selfClass, rank, polygons, vertices ):
        """Constructor with sizes."""
        ( polygonNormals, polygons, vertexNormals, vertices ) = selfClass.AllocateArrays ( rank, polygons, vertices )
        return selfClass.FromArrays ( polygonNormals, polygons, vertexNormals, vertices )

    # . Properties.
    @property
    def  rank ( self ): return ( self.vertices.shape[1] )

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass
