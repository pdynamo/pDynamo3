"""Classes and functions for reading OOGL Off files."""

from  pCore            import logFile, LogFileActive, TextFileReader
from .PolygonalSurface import PolygonalSurface

#===================================================================================================================================
# . OOGLOff file reader class.
#===================================================================================================================================
class OOGLOffFileReader ( TextFileReader ):
    """OOGLOffFileReader is the class for OOGL Off files that are to be read."""

    _classLabel = "OOGL Off File Reader"

    def Parse ( self, log = logFile ):
        """Parse the data on the file."""
        if not self.isParsed:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
            # . Open the file.
            self.Open ( )
            try:
                # . Header line.
                line       = self.GetLine ( )
                hasNormals = ( line.find ( "N" ) > -1 )
                if line.find ( "n" ) > -1:
                    items = self.GetTokens ( converters = [ int ] )
                    rank  = items[0]
                else:
                    rank  = 3
                # . Counters line.
                items     = self.GetTokens ( converters = [ int, int, int ] )
                nVertices = items[0]
                nPolygons = items[1]
                # . Allocate the object.
                self.surface = PolygonalSurface.WithSizes ( rank, nPolygons, nVertices )
                normals      = self.surface.vertexNormals
                polygons     = self.surface.polygons
                vertices     = self.surface.vertices
                # . Empty line.
                self.GetLine ( )
                # . Vertex lines.
                nFields = rank
                if hasNormals: nFields *= 2
                converters = nFields * [ float ]
                for i in range ( nVertices ):
                    items = self.GetTokens ( converters = converters )
                    for c in range ( rank ): vertices[i,c] = items[c]
                    if hasNormals:
                        for c in range ( rank ): normals[i,c] = items[c+rank]
                # . Empty line.
                self.GetLine ( )
                # . Face lines.
                nFields    = rank + 1
                converters = nFields * [ int ]
                for i in range ( nPolygons ):
                    items = self.GetTokens ( converters = converters )
                    for c in range ( rank ): polygons[i,c] = items[c+1]
            except EOFError:
                pass
            # . Close the file.
            self.WarningStop ( )
            self.Close ( )
            # . Set the parsed flag and some other options.
            self.log     = None
            self.isParsed = True

    @classmethod
    def PathToPolygonalSurface ( selfClass, path, log = logFile ):
        """Read a surface from an OOGL Off file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse ( )
        return inFile.ToPolygonalSurface ( )

    def ToPolygonalSurface ( self ):
        """Return a polygonal surface."""
        surface = None
        if self.isParsed: surface = self.surface
        return surface

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
