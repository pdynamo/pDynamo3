"""Classes and functions for writing ASCII STL files."""

import math

from  pCore            import logFile          , \
                              LogFileActive    , \
                              TextFileWriter
from .PolygonalSurface import PolygonalSurface

# . Need to add option for writing multiple surfaces to the same file (as in STLFileReader).

# . The STL file requires that the vertices of each facet are output in counterclockwise order as viewed from the exterior of the
# . surface. Likewise, the facet normals should be pointing outwards. This is not automatically done in pDynamo - it depende on
# . how the surface is generated. However, there are several tools available to correct the surface (e.g. admesh).

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
_DefaultName = "Surface"

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class STLFileWriter ( TextFileWriter ):
    """Class for writing ASCII STL files."""

    _classLabel = "STL File Writer"

    @classmethod
    def PathFromPolygonalSurface ( selfClass, path, surface, log = logFile, name = None ):
        """Write a polygonal surface to an STL file."""
        outFile = selfClass.FromPath  ( path )
        outFile.WritePolygonalSurface ( surface, name = name )

    def WritePolygonalSurface ( self, surface, format = "{:.6e} ", name = None ):
        """Write a surface."""
        if name is None: name = _DefaultName
        normals  = surface.polygonNormals
        polygons = surface.polygons
        vertices = surface.vertices
        self.Open ( )
        self.file.write ( "solid {:s}\n".format ( name ) )
        for p in range ( polygons.rows ):
            self.file.write ( "facet normal " )
            for c in range ( 3 ): self.file.write ( format.format ( normals[p,c] ) )
            self.file.write ( "\n    outer loop" )
            for v in polygons[p,:]:
                self.file.write ( "\n        vertex " )
                for c in range ( 3 ): self.file.write ( format.format ( vertices[v,c] ) )
            self.file.write ( "\n    endloop\nendfacet\n" )
        self.file.write ( "endsolid {:s}\n".format ( name ) )
        self.Close ( )

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass
