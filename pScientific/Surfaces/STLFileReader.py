"""Classes and functions for reading ASCII STL files."""

from  pCore            import logFile, LogFileActive, TextFileReader
from .PolygonalSurface import PolygonalSurface
from .SurfacesError    import SurfacesError

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
 # . The vertex format is a reasonable way of getting unique vertices without calculating distances.
_DefaultVertexFormat = "{:12.6e} {:12.6e} {:12.6e}"
_DefaultName         = "Surface"

#===================================================================================================================================
# . STL file reader class.
#===================================================================================================================================
class STLFileReader ( TextFileReader ):
    """STLFileReader is the class for STL files that are to be read."""

    _classLabel = "STL File Reader"

    def Parse ( self, log = logFile ):
        """Parse the data on the file."""
        if not self.isParsed:
            if LogFileActive ( log ): self.log = log
            self.Open ( )
            facets       = None
            isTerminated = False
            name         = None
            surfaces     = []
            while True:
                try:
                    # . Header line.
                    tokens = self.GetTokens ( signalWarnings = False )
                    if tokens[0] != "solid" : raise SurfacesError ( "Missing \"solid\" keyword in STL file." )
                    if len ( tokens ) > 1: name = tokens[1]
                    else:                  name = _DefaultName
                    isTerminated = False
                    # . Facets.
                    facets = []
                    while True:
                        tokens = self.GetTokens ( )
                        if ( len ( tokens ) >= 5 ) and ( tokens[0] == "facet" ) and ( tokens[1] == "normal" ):
                            facet  = [ [ float ( t ) for t in tokens[2:5] ] ]
                            line   = self.GetLine ( )
                            if line == "outer loop":
                                for v in range ( 3 ):
                                    tokens = self.GetTokens ( converters = [ None, float, float, float ] )
                                    if tokens.pop ( 0 ) == "vertex": facet.append ( tokens )
                                    else: self.Warning ( "Invalid \"vertex\" keyword line.", True )
                                facets.append ( facet )
                                if self.GetLine ( ) != "endloop" : self.Warning ( "Invalid \"endloop\" keyword line." , True )
                                if self.GetLine ( ) != "endfacet": self.Warning ( "Invalid \"endfacet\" keyword line.", True )
                            else:
                                self.Warning ( "Invalid \"outer loop\" keyword line.", True )
                        elif ( len ( tokens ) >= 1 ) and ( tokens[0] == "endsolid" ):
                            surfaces.append ( ( name, facets ) )
                            facets       = []
                            name         = None
                            isTerminated = True
                            break
                        else:
                            self.Warning ( "Missing \"endsolid\" or \"facet normal\" keyword line.", True )
                except EOFError:
                    break
            if not isTerminated:
                self.Warning ( "Missing \"endsolid\" keyword in STL file.", True )
                if ( name   is not None ) and \
                   ( facets is not None ) and \
                   ( len ( facets ) > 0 ): surfaces.append ( ( name, facets ) )
            self.WarningStop     ( )
            self.Close           ( )
            self.isParsed = True
            self.log      = None
            self.ProcessSurfaces ( surfaces )

    @classmethod
    def PathToPolygonalSurface ( selfClass, path, log = logFile, name = None ):
        """Read a surface from an STL file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse ( )
        return inFile.ToPolygonalSurface ( name = name )

    def ProcessSurfaces ( self, surfaces ):
        """Process the raw input surfaces."""
        names = set ( [ name for ( name, facets ) in surfaces ] )
        if len ( names ) != len ( surfaces ): raise SurfacesError ( "Duplicate \"solid\" names in STL file." )
        self.surfaces = {}
        for ( name, facets ) in surfaces:
            polygons = []
            strings  = {}
            vertices = []
            for facet in facets:
                polygon = []
                for ( x, y, z ) in facet[1:4]:
                    string = _DefaultVertexFormat.format ( x, y, z )
                    vIndex = strings.get ( string, -1 )
                    if vIndex < 0:
                        vIndex          = len ( vertices )
                        strings[string] = vIndex
                        vertices.append ( ( x, y, z ) )
                    polygon.append ( vIndex )
                polygons.append ( polygon )
            surface       = PolygonalSurface.WithSizes ( 3, len ( polygons ), len ( vertices ) )
            surface.label = name
            for ( p, ( i, j, k ) ) in enumerate ( polygons ):
                surface.polygons[p,0] = i
                surface.polygons[p,1] = j
                surface.polygons[p,2] = k
            for ( v, ( x, y, z ) ) in enumerate ( vertices ):
                surface.vertices[v,0] = x
                surface.vertices[v,1] = y
                surface.vertices[v,2] = z
            surface.MakePolygonNormals ( )
            surface.MakeVertexNormalsFromPolygonalNormals ( )
            self.surfaces[name] = surface

    def ToPolygonalSurface ( self, name = _DefaultName ):
        """Return a polygonal surface."""
        surface = None
        if self.isParsed:
            if len ( self.surfaces ) == 1: surface = list ( self.surfaces.values ( ) )[0]
            else:                          surface = self.surfaces.get ( name, None )
        return surface

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
