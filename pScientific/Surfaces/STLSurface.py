"""A class for STL surfaces."""

# . A basic (and slow!) class for STL surfaces.
# . Mainly designed for checking surfaces for consistency.
# . More advanced functionality, such as repairing, can be found in admesh and other programs.

import random

from   collections      import defaultdict
from   pCore            import Align            , \
                               logFile          , \
                               LogFileActive
from ..Geometry3        import Coordinates3     , \
                               Vector3
from  .PolygonalSurface import PolygonalSurface
from  .SurfacesError    import SurfacesError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class STLSurface:
    """A STL surface."""

    def _CheckEdgeOrientations ( self ):
        """Check edge orientations."""
        # . Returns the number of paired edges with parallel, as opposed to antiparallel, directions.
        n        = 0
        polygons = self.surface.polygons
        for ( ( v1, v2 ), ( p, q ) ) in self.pairedEdges.items ( ):
            pV = list ( polygons[p,:] ) ; p1 = pV.index ( v1 ) ; p2 = pV.index ( v2 ) ; p12 = ( p1 - p2 ) % 3
            qV = list ( polygons[q,:] ) ; q1 = qV.index ( v1 ) ; q2 = qV.index ( v2 ) ; q12 = ( q1 - q2 ) % 3
            if p12 == q12: n += 1
        return n

    def _FindOpenComponents ( self ):
        """Find which components are open."""
        neighbors = self.polygonEdgeNeighbors
        unpaired  = { p for pq in self.unpairedEdges.values ( ) for p in pq }
        flags     = []
        for component in self.connectedComponents:
            common     = unpaired & set ( component )
            nNeighbors = [ len ( neighbors[p] ) for p in component ]
            flags.append ( ( len ( common ) > 0 ) or any ( n != 3 for n in nNeighbors ) )
        return flags

    def _MakeAreas ( self ):
        """The surface area."""
        a        = Vector3.Null ( )
        b        = Vector3.Null ( )
        polygons = self.surface.polygons
        vertices = self.surface.vertices
        areas = []
        for component in self.connectedComponents:
            area = 0.0
            for p in component:
                v1 = polygons[p,0]
                v2 = polygons[p,1]
                v3 = polygons[p,2]
                for c in range ( 3 ):
                    a[c] = vertices[v2,c] - vertices[v1,c]
                    b[c] = vertices[v3,c] - vertices[v1,c]
                a.Cross ( b ) # . The unnormalized normal.
                area += a.Norm2 ( )
            areas.append ( 0.5 * area )
        return areas

    def _MakeComponents ( self ):
        """Find the connected components."""
        nP         = self.surface.polygons.rows
        neighbors  = self.polygonEdgeNeighbors
        isChecked  = [ False for p in range ( nP ) ]
        components = []
        for p in range ( nP ):
            if not isChecked[p]:
                isChecked[p] = True
                component    = [ p ]
                toCheck      = [ p ]
                components.append ( component )
                while len ( toCheck ) > 0:
                    q = toCheck.pop ( 0 )
                    for n in neighbors[q]:
                        if not isChecked[n]:
                            isChecked[n] = True
                            component.append ( n )
                            toCheck.append   ( n )
                component.sort ( )
        return components

    def _MakeEdges ( self ):
        """Find dictionaries of paired and unpaired edges."""
        polygons = self.surface.polygons
        edges    = defaultdict ( list )
        for p in range ( polygons.rows ):
            vertices       = list ( polygons[p,:] )
            ( v1, v2, v3 ) = sorted ( vertices ) # . Assumes no vertices identical! Maybe should check.
            edges[(v1,v2)].append ( p )
            edges[(v1,v3)].append ( p )
            edges[(v2,v3)].append ( p )
        paired   = { key : value for ( key, value ) in edges.items ( ) if len ( value ) == 2 }
        unpaired = { key : value for ( key, value ) in edges.items ( ) if len ( value ) != 2 }
        return ( paired, unpaired )

    def _MakePolygonEdgeNeighbors ( self ):
        """Find the polygon edge neighbors."""
        neighbors = [ [] for p in range ( self.surface.polygons.rows ) ]
        if len ( self.unpairedEdges ) > 0:
            for pq in self.unpairedEdges.values ( ):
                for n in pq: neighbors[n].extend ( pq )
            neighbors = [ sorted ( set ( n ) ) for n in neighbors ]
        for ( p, q ) in self.pairedEdges.values ( ):
            neighbors[p].append ( q )
            neighbors[q].append ( p )
        return neighbors

    def _MakeSubsurfaces ( self ):
        """Split the surface into separate subsurfaces (= connected components)."""
        components = self.connectedComponents
        if len ( components ) == 1:
            subsurfaces = [ self.surface ]
        else:
            surface  = self.surface
            polygons = self.surface.polygons
            if surface.label is None: label = "component"
            else:                     label = surface.label
            subsurfaces = []
            for ( c, component ) in enumerate ( components ):
                vertices = set ( )
                for p in component:
                    vertices.update ( polygons[p,:] )
                vertices = sorted ( vertices )
                vIndices = { p : i for ( i, p ) in enumerate ( vertices ) }
                subsurface       = PolygonalSurface.WithSizes ( 3, len ( component ), len ( vertices ) )
                subsurface.label = "{:s}_{:d}".format ( label, c )
                for ( i, p ) in enumerate ( component ):
                    for c in range ( 3 ):
                        subsurface.polygonNormals[i,c] = surface.polygonNormals[p,c]
                        subsurface.polygons      [i,c] = vIndices[surface.polygons[p,c]]
                for ( i, v ) in enumerate ( vertices ):
                    for c in range ( 3 ):
                        subsurface.vertexNormals[i,c] = surface.vertexNormals[v,c]
                        subsurface.vertices     [i,c] = surface.vertices     [v,c]
                subsurfaces.append ( subsurface )
        return subsurfaces

    def _MakeVolumes ( self ):
        """The volume."""
        # . The volumes of open components are always zero.
        a        = Vector3.Null ( )
        b        = Vector3.Null ( )
        h        = Vector3.Null ( )
        o        = Vector3.Null ( )
        polygons = self.surface.polygons
        vertices = self.surface.vertices
        volumes  = []
        for ( component, isOpen ) in zip ( self.connectedComponents, self.openComponents ):
            volume = 0.0
            if not isOpen:
                vO     = polygons[component[0],0] # . Arbitrary choice of origin.
                for c in range ( 3 ): o[c] = vertices[vO,c]
                for p in component:
                    v1 = polygons[p,0]
                    v2 = polygons[p,1]
                    v3 = polygons[p,2]
                    for c in range ( 3 ):
                        a[c] = vertices[v2,c] - vertices[v1,c]
                        b[c] = vertices[v3,c] - vertices[v1,c]
                        h[c] = vertices[v1,c] - o[c]
                    a.Cross ( b ) # . Unnormalized normal.
                    volume += a.Dot ( h ) # . Sign is important.
                volume /= 6.0
            volumes.append ( volume )
        return volumes

    def _ReorientPolygons ( self ):
        """Reorient the polygons to have the same orientation."""
        polygons     = self.surface.polygons
        edges        = self.pairedEdges.copy ( )
        nP           = polygons.rows
        # . Make polygon edges.
        polygonEdges = []
        for p in range ( nP ):
            vertices       = list ( polygons[p,:] )
            ( v1, v2, v3 ) = sorted ( vertices )
            polygonEdges.append ( ( (v1,v2), (v1,v3), (v2,v3) ) )
        # . Loop over edges and check that they have the same orientation in their host polygons.
        # . Treat each connected component separately.
        # . A polygon can only be swapped once - when it's first edge is checked.
        nSwaps = 0
        for component in self.connectedComponents:
            toCheck = [ component[0] ]
            while len ( toCheck ) > 0:
                p = toCheck.pop ( 0 )
                for key in polygonEdges[p]:
                    pq = edges.pop ( key, None )
                    if pq is not None:
                        if pq[0] == p: q = pq[1]
                        else:          q = pq[0]
                        toCheck.append ( q )
                        ( v1, v2 ) = key
                        pV = list ( polygons[p,:] ) ; p1 = pV.index ( v1 ) ; p2 = pV.index ( v2 ) ; p12 = ( p1 - p2 ) % 3
                        qV = list ( polygons[q,:] ) ; q1 = qV.index ( v1 ) ; q2 = qV.index ( v2 ) ; q12 = ( q1 - q2 ) % 3
                        if p12 == q12:
                            t       = polygons[q,q1] ; polygons[q,q1] = polygons[q,q2] ; polygons[q,q2] = t
                            nSwaps += 1
        return nSwaps

    @classmethod
    def FromPolygonalSurface ( selfClass, surface ):
        """Constructor from a polygonal surface."""
        if surface.rank != 3: raise SurfacesError ( "Invalid surface for STL surface construction." )
        self         = selfClass ( )
        self.surface = surface
        return self

    def MakePolygonCentroids ( self, normalStep = None ):
        """Return the polygon centroids with an optional step along the normals."""
        polygons = self.surface.polygons
        vertices = self.surface.vertices
        nP       = polygons.rows
        xyz      = Coordinates3.WithExtent ( nP ) ; xyz.Set ( 0.0 )
        for p in range ( nP ):
            for r in range ( 3 ):
                xyz[p,:].Add ( vertices[polygons[p,r],:] )
            xyz[p,:].Scale ( 1.0 / 3.0 )
        if normalStep is not None:
            xyz.Add ( self.surface.polygonNormals, scale = normalStep )
        return xyz

    def OrientPolygons ( self, log = logFile ):
        """Orient the polygons."""
        # . Antiparallel edges.
        parallelEdges0 = self._CheckEdgeOrientations ( )
        isOK           = ( parallelEdges0 == 0 )
        messages       = [ ( "Wrongly oriented edges", "{:d}".format ( parallelEdges0 ) ) ]
        if not isOK:
            reorientations = self._ReorientPolygons      ( )
            parallelEdges  = self._CheckEdgeOrientations ( )
            isOK           = ( parallelEdges == 0 )
            messages.extend ( [ ( "Number of reorientations"                  , "{:d}".format ( reorientations ) ) ,
                                ( "Wrongly oriented edges after reorientation", "{:d}".format ( parallelEdges  ) ) ] )
        # . Negative volumes - only for closed surfaces.
        if isOK:
            self.surface.MakePolygonNormals ( )
            self.__dict__["_volumes"] = None
            negatives0 = 0
            polygons   = self.surface.polygons
            for ( component, volume ) in zip ( self.connectedComponents, self.volumes ):
                if volume < 0.0:
                    negatives0 += 1
                    for p in component: # . Swap last two vertices.
                        t             = polygons[p,1]
                        polygons[p,1] = polygons[p,2]
                        polygons[p,2] = t
            isOK = ( negatives0 == 0 )
            messages.append ( ( "Number of negative volumes", "{:d}".format ( negatives0 ) ) )
            if not isOK:
                self.surface.MakePolygonNormals ( )
                self.__dict__["_volumes"] = None
                negatives = sum ( [ 1 for v in self.volumes if v < 0.0 ] )
                isOK      = ( negatives == 0 )
                messages.append ( ( "Number of negative volumes after correction", "{:d}".format ( negatives ) ) )
        # . Reporting.
        if LogFileActive ( log ):
            if ( parallelEdges0 == 0 ) and ( negatives0 == 0 ):
                if self.surface.label is None:
                    log.Paragraph ( "STL surface correctly oriented." )
                else:
                    log.Paragraph ( "STL surface \"{:s}\" correctly oriented.".format ( self.surface.label ) )
            else:
                title = "STL Surface Polygon Orientation"
                if self.surface.label is not None: title = "{:s} For \"{:s}\"".format ( title, self.surface.label )
                messages.append ( ( "Surface correctly oriented", "{:s}".format ( repr ( isOK ) ) ) )
                n = 0
                for ( key, value ) in messages: n = max ( len ( key ), n )
                table = logFile.GetTable ( columns = [ 3, 3, n+5, 10 ] )
                table.Start ( )
                table.Title ( title )
                for ( i, ( key, value ) ) in enumerate ( messages ):
                    table.Entry ( repr ( i+1 ) )
                    table.Entry ( "" )
                    table.Entry ( key   , align = Align.Left )
                    table.Entry ( value )
                table.Stop ( )
        return isOK

    def RandomizePolygonVertices ( self ):
        """Randomize the vertices of each polygon."""
        polygons = self.surface.polygons
        for p in range ( polygons.rows ):
            vertices = list ( polygons[p,:] )
            random.shuffle ( vertices )
            polygons[p,0] = vertices[0]
            polygons[p,1] = vertices[1]
            polygons[p,2] = vertices[2]

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( logFile ):
            surface = self.surface
            if surface.label is None: title = "STL Surface Summary"
            else:                     title = "Summary for STL Surface \"{:s}\"".format ( surface.label )
            log.SummaryOfItems ( self.SummaryItems ( ), title = title )

    def SummaryItems ( self ):
        """Summary items."""
        return [ ( "Components"     , "{:d}".format   ( len ( self.connectedComponents   ) ) ) ,
                 ( "Open Components", "{:d}".format   ( self.openComponents.count ( True ) ) ) ,
                 ( "Polygons"       , "{:d}".format   ( self.surface.polygons.rows         ) ) ,
                 ( "Vertices"       , "{:d}".format   ( self.surface.vertices.rows         ) ) ,
                 ( "Area"           , "{:.1f}".format ( self.area                          ) ) ,
                 ( "Volume"         , "{:.1f}".format ( self.volume                        ) ) ]

    @property
    def area ( self ):
        return sum ( self.areas )

    @property
    def areas ( self ):
        items = self.__dict__.get ( "_areas", None )
        if items is None:
            items = self._MakeAreas ( )
            self.__dict__["_areas"] = items
        return items

    @property
    def connectedComponents ( self ):
        items = self.__dict__.get ( "_components", None )
        if items is None:
            items = self._MakeComponents ( )
            self.__dict__["_components"] = items
        return items

    @property
    def isClosed ( self ):
        """Is the surface closed?"""
        return not self.isOpen

    @property
    def isOpen ( self ):
        """Is the surface open?"""
        return any ( self.openComponents )

    @property
    def openComponents ( self ):
        """Flags to indicate whether a component is open or closed."""
        items = self.__dict__.get ( "_openComponents", None )
        if items is None:
            items = self._FindOpenComponents ( )
            self.__dict__["_openComponents"] = items
        return items

    @property
    def pairedEdges ( self ):
        paired = self.__dict__.get ( "_pairedEdges", None )
        if paired is None:
            ( paired, unpaired ) = self._MakeEdges ( )
            self.__dict__["_pairedEdges"  ] =   paired
            self.__dict__["_unpairedEdges"] = unpaired
        return paired

    @property
    def polygonEdgeNeighbors ( self ):
        items = self.__dict__.get ( "_polygonEdgeNeighbors", None )
        if items is None:
            items = self._MakePolygonEdgeNeighbors ( )
            self.__dict__["_polygonEdgeNeighbors"] = items
        return items

    @property
    def subsurfaces ( self ):
        items = self.__dict__.get ( "_subsurfaces", None )
        if items is None:
            items = self._MakeSubsurfaces ( )
            self.__dict__["_subsurfaces"] = items
        return items

    @property
    def unpairedEdges ( self ):
        unpaired = self.__dict__.get ( "_unpairedEdges", None )
        if unpaired is None:
            ( paired, unpaired ) = self._MakeEdges ( )
            self.__dict__["_pairedEdges"  ] =   paired
            self.__dict__["_unpairedEdges"] = unpaired
        return unpaired

    @property
    def volume ( self ):
        return sum ( self.volumes )

    @property
    def volumes ( self ):
        items = self.__dict__.get ( "_volumes", None )
        if items is None:
            items = self._MakeVolumes ( )
            self.__dict__["_volumes"] = items
        return items

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass
