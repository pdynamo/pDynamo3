"""Classes and functions for writing OOGL Off files."""

import math

from  pCore            import logFile          , \
                              LogFileActive    , \
                              TextFileWriter
from .PolygonalSurface import PolygonalSurface

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class OOGLOffFileWriter ( TextFileWriter ):
    """Class for writing OOGL Off files."""

    _classLabel = "OOGL Off File Writer"

    @classmethod
    def PathFromPolygonalSurface ( selfClass, path, surface, log = logFile ):
        """Write a polygonal surface to an OOGL Off file."""
        outFile = selfClass.FromPath  ( path    )
        outFile.WritePolygonalSurface ( surface )

    def WritePolygonalSurface ( self, surface ):
        """Write a surface."""
        # . Initialization.
        self.Open ( )
        polygons = surface.polygons
        rank     = surface.rank
        vertices = surface.vertices
        # . Header.
        if rank == 3: self.file.write ( "NOFF\n" )
        else:         self.file.write ( "NnOFF\n{:d}\n".format ( rank ) )
        # . Number of vertices, faces and edges (latter not used).
        self.file.write ( "{:d} {:d} 0\n\n".format ( vertices.rows, polygons.rows ) )
        # . Vertices.
        format = 2 * rank * " {:.10f}" + "\n"
        for i in range ( vertices.rows ):
            self.file.write ( format.format ( vertices[i,:] ) )
        # . Polygons.
        self.file.write ( "\n" )
        format = "{:d}".format ( rank ) + rank * " {:d}" + "\n"
        for i in range ( polygons.rows ):
            self.file.write ( format.format ( polygons[i,:] ) )
        # . Finish up.
        self.Close ( )

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass
