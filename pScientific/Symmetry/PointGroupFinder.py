"""Point group tests."""

import glob, math, os

from   collections       import defaultdict
from   pCore             import Align                 , \
                                AttributableObject    , \
                                Clone                 , \
                                logFile               , \
                                LogFileActive         , \
                                Selection
from  .PointGroup        import PointGroups_FromYAML
from  .SymmetryError     import SymmetryError
from  .SymmetryOperation import Identity              , \
                                ImproperRotation      , \
                                Inversion             , \
                                ProperRotation        , \
                                Reflection
from ..Arrays            import Array                 , \
                                ArrayPrint2D
from ..Geometry3         import Vector3
from ..Graph             import BiconnectedComponents , \
                                Edge                  , \
                                EdgeVectorToPath      , \
                                Graph                 , \
                                Node                  , \
                                PathToEdgeVector      , \
                                VismaraRelevantCycles

#
# . This algorithm is quite slow for very high symmetry 3D-graphs (e.g. icosahedral groups) due to the ring finding. It may be
# . possible to optimize it by not including every group of connections in PostulateHigherOrderRotationAxes. However, simple
# . restrictions (like checking the length of the number of connections of a given distance) do not work.
#
# . Possible additions include:
#
# - Better way of handling indeterminate possible C2 axes.
# - Flexible point group assignment (if keys don't match).
# - Tolerances higher for larger Cn.
#
# - Include a pointGroupFinderState to avoid putting results into Finder itself.
#

#
# . Irreducible representation determination is quite sensitive to the degeneracy tolerances. Too high a tolerance bunches
# . vectors or states together and then illegal combinations of IRs can occur. This might be improved by a better algorithm.
# . As high an accuracy as possible is best in the CI, SCF or normal mode calculation.
#
# . Accidental degeneracies are not yet handled.
#
# . Also problems will arise if the wavefunction does not have the symmetry of the nuclear framework (possible in some cases,
# . particularly when there are distorted geometries).
#

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Options.
_IncludeCompositeImproperRotations = False

# . Tolerances.
# . The results can be quite sensitive to these tolerances.
_AngleTolerance          = 0.1 # . In radians.
_DistanceTolerance       = 0.1 # . In Angstroms.

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class PointGroupFinder ( AttributableObject ):
    """Class for finding the point group of an object."""

    # . Attributes.
    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "angleTolerance"    : _AngleTolerance    ,
                             "distanceTolerance" : _DistanceTolerance ,
                             "pointGroups"       : None               } )

    def _CheckOptions ( self ):
        """Check options."""
        super ( PointGroupFinder, self )._CheckOptions ( )
        # . Set additional tolerances.
        self.cosineParallelTolerance      = math.fabs ( 1.0 - math.cos ( self.angleTolerance ) )
        self.cosinePerpendicularTolerance = math.fabs (       math.sin ( self.angleTolerance ) )
        # . Point groups.
        self.pointGroups = PointGroups_FromYAML ( )

    #===============================================================================================================================
    def AddAdditionalCycles ( self, graph, cycles ):
        """Add cycles by fusing existing cycles."""

        # . Do something only if there are sufficient cycles.
        if len ( cycles ) > 1:

            # . Initialization.
            newCycles = []

            # . Find the cycles each node is involved in.
            maximumCycleLength = 0
            nodeCycles         = {}
            for cycle in cycles:
                maximumCycleLength = max ( maximumCycleLength, len ( cycle ) )
                for node in cycle:
                    existing = nodeCycles.get ( node, [] )
                    existing.append ( cycle )
                    nodeCycles[node] = existing

            # . Loop over nodes.
            # . Only do something if there are sufficient cycles to fuse.
            for ( node, local ) in nodeCycles.items ( ):
                if len ( local ) > maximumCycleLength:

                    # . Find the indices of the edges that this node is involved in.
                    edgeIndices = []
                    for ( e, edge ) in enumerate ( graph.edges ):
                        if ( edge.node1 is node ) or ( edge.node2 is node ): edgeIndices.append ( e )

                    # . Double loop over cycles.
                    for ( c, cycle ) in enumerate ( local[:-1] ):
                        fusedCycle = PathToEdgeVector ( graph, cycle, closePath = True )
                        for other in local[c+1:]:
                            edgeVector = PathToEdgeVector ( graph, other, closePath = True )

                            # . Check for common edges at node.
                            overlap = False
                            for e in edgeIndices:
                                if edgeVector[e]:
                                    overlap = True
                                    break

                            # . Fuse cycles and then check for edges at node.
                            if overlap:
                                for i in range ( len ( edgeVector ) ):
                                    fusedCycle[i] ^= edgeVector[i]
                                found = False
                                for e in edgeIndices:
                                    if fusedCycle[e]:
                                        found = True
                                        break
                                if not found: break

                        # . Check the fused cycle.
                        # . Only include a path if it is a cycle with length > maximumCycleLength and that does not contain the node.
                        found = False
                        for e in edgeIndices:
                            if fusedCycle[e]:
                                found = True
                                break
                        if not found:
                            try:
                                newCycle = EdgeVectorToPath ( graph, fusedCycle )[:-1]
                                if len ( newCycle ) > maximumCycleLength: newCycles.append ( newCycle )
                            except:
                                pass

            # . Finish up.
            cycles.extend ( newCycles )

    #===============================================================================================================================
    def CheckForAlignment ( self, order, iRing, iCenter, jRing, jCenter, coordinates3 ):
        """Check two rings for an alignment which is appropriate for an improper rotation."""
        # . Initialization.
        areAligned    = False
        improperOrder = 0
        # . Get minimum absolute angle between a center-node vector of one ring and all those of the other ring.
        iVector = Vector3.Null ( )
        jVector = Vector3.Null ( )
        iNode   = iRing[0]
        for i in range ( 3 ): iVector[i] = coordinates3[iNode,i] - iCenter[i]
        iVector.Normalize ( )
        angles = []
        for jNode in jRing:
            for i in range ( 3 ): jVector[i] = coordinates3[jNode,i] - jCenter[i]
            jVector.Normalize ( )
            dot = iVector.Dot ( jVector )
            if   dot >  1.0: dot =  1.0
            elif dot < -1.0: dot = -1.0
            angles.append ( math.fabs ( math.acos ( dot ) ) )
        angles.sort ( )
        minimumAngle = angles[0]
        # . Check for alignment.
        areEclipsed  = ( minimumAngle <= self.angleTolerance )
        areStaggered = ( math.fabs ( minimumAngle - math.pi / float ( order ) ) <= self.angleTolerance )
        assert not ( areEclipsed and areStaggered ), "Logic error."
        # . Staggered rings imply an S2n that cannot be generated from the ring Cn and sigmah.
        if areStaggered:
            areAligned    = True
            improperOrder = 2 * order
        # . Eclipsed rings imply existing Cn and sigmah operations.
        elif areEclipsed and _IncludeCompositeImproperRotations:
            areAligned    = True
            improperOrder = order
        return ( areAligned, improperOrder )

    #===============================================================================================================================
    def CheckForPlanarity ( self, nodes, coordinates3 ):
        """Check for planarity.

        Calculate the axis parameters if the ring is planar.
        """
        # . Initialization.
        axis     = None
        center   = None
        distance = None
        isPlanar = True
        # . Determine planarity.
        n        = len ( nodes )
        target   = math.pi * ( 1.0 - 2.0 / float ( n ) )
        for j in range ( n ):
            jIndex = nodes[j]
            if j == 0:
                iIndex = nodes[-1]
                kIndex = nodes[ 1]
            elif j == ( n - 1 ):
                iIndex = nodes[j-1]
                kIndex = nodes[0]
            else:
                iIndex = nodes[j-1]
                kIndex = nodes[j+1]
            angle = math.radians ( coordinates3.Angle ( iIndex, jIndex, kIndex ) )
            if math.fabs ( angle - target ) > self.angleTolerance:
                isPlanar = False
                break
        # . Determine axis parameters.
        if isPlanar:
            # . Calculate the coordinates of the ring center.
            axis   = Vector3.Null ( )
            center = Vector3.Null ( )
            for iNode in nodes:
                for i in range ( 3 ): center[i] += coordinates3[iNode,i]
            for i in range ( 3 ): center[i] /= float ( n )
            distance = center.Norm2 ( )
            # . The center is at the origin so determine the ring normal.
            # . Arbitrarily use the coordinates of the first two nodes in the ring.
            # . Direction is also arbitrary.
            if distance <= self.distanceTolerance:
                iNode    = nodes[0]
                jNode    = nodes[1]
                axis[0]  = coordinates3[iNode,1]*coordinates3[jNode,2] - coordinates3[jNode,1]*coordinates3[iNode,2]
                axis[1]  = coordinates3[iNode,2]*coordinates3[jNode,0] - coordinates3[jNode,2]*coordinates3[iNode,0]
                axis[2]  = coordinates3[iNode,0]*coordinates3[jNode,1] - coordinates3[jNode,0]*coordinates3[iNode,1]
                axis.Normalize ( )
                distance = 0.0
                for i in range ( 3 ): center[i] = 0.0
            # . Normalize.
            else:
                for i in range ( 3 ): axis[i] = center[i]
                axis.Scale ( 1.0 / distance )
        return ( isPlanar, distance, axis, center )

    #===============================================================================================================================
    def FindPlanarRings ( self, connections, coordinates3 ):
        """Find all planar rings given a set of connections between nodes."""
        # . Initialization.
        planarRings = {}
        # . Create the graph.
        graph = Graph ( )
        nodes = {}
        for connection in connections:
            for i in connection:
                if i not in nodes:
                    node = Node ( )
                    node.index = i
                    nodes[i] = node
                    graph.AddNode ( node )
        for ( i, j ) in connections:
            graph.AddEdge ( Edge.WithNodes ( nodes[i], nodes[j] ) )
        # . Get the cycles.
        cycles = VismaraRelevantCycles ( graph )
        # . Add additional cycles by fusing the relevant cycles for each node.
        self.AddAdditionalCycles ( graph, cycles )
        # . Loop over cycles.
        for cycle in cycles:
            ring = []
            for node in cycle: ring.append ( node.index )
            ( isPlanar, distance, axis, center ) = self.CheckForPlanarity ( ring, coordinates3 )
            if isPlanar:
                n     = len ( ring )
                rings = planarRings.get ( n, [] )
                rings.append ( ( ring, distance, axis, center ) )
                planarRings[n] = rings
        # . Finish up.
        return planarRings

    #===============================================================================================================================
    def Find3DGraphPointGroup ( self, nodeTypes, coordinates3, doCharacterSymmetryOperations = True, log = logFile, weights = None ):
        """Find the point group of a 3D graph."""
        # . Initialization.
        symmetryOperations = {}
        if weights is None:
            weights = Array.WithExtent ( len ( nodeTypes ) )
            for ( i, t ) in enumerate ( nodeTypes ): weights[i] = float ( t )

        # . Translate the coordinates to their weighted center.
        center = coordinates3.Center ( weights = weights )
        center.Scale ( -1.0 )
        coordinates3.Translate ( center )

        # . Classify nodes into possible symmetry classes.
        nodeGroups = self.PostulateSymmetryEquivalentNodes ( nodeTypes, coordinates3, log = log )

        # . Check the moments of inertia for basic information about the graph.
        self.ProcessMomentsOfInertia ( nodeGroups, coordinates3, weights, symmetryOperations )

        # . Find higher-order operations.
        ( properRotations, improperRotations ) = self.PostulateHigherOrderRotationAxes ( nodeGroups, coordinates3, log = log )
        self.ProcessHigherOrderRotationAxes ( nodeGroups, coordinates3, properRotations, improperRotations, symmetryOperations )

        # . Find order-2 operations (C2 rotation, inversion, reflection).
        hasIndeterminateC2Axes = self.ProcessOrderTwoOperations ( nodeGroups, coordinates3, symmetryOperations )
        if hasIndeterminateC2Axes: self.ProcessIndeterminateC2Axes ( nodeGroups, coordinates3, symmetryOperations )

        # . Remove redundant S4s.
        self.RemoveRedundantS4s ( symmetryOperations )

        # . Append the identity.
        symmetryOperations["E"] = [ Identity.WithDefaults ( ) ]

        # . Printing.
        keys = list ( symmetryOperations.keys ( ) )
        keys.sort ( )
        if LogFileActive ( log ):
            if len ( symmetryOperations ) > 0:
                table = log.GetTable ( columns = [ 20, 20 ] )
                table.Start  ( )
                table.Title  ( "Symmetry Operations" )
                for key in keys:
                    table.Entry ( key, align = Align.Left )
                    table.Entry ( "{:d}".format ( len ( symmetryOperations[key] ) ) )
                table.Stop ( )
            else: log.Paragraph ( "No symmetry operations identified." )

        # . Generate the operation key.
        items = []
        for key in keys:
            n = len ( symmetryOperations[key] )
            if n == 1: items.append ( key )
            else:      items.append ( repr ( n ) + "*" + key )
        operationKey = " ".join ( items )

        # . Identify point group from a point group dictionary.
        # . Maybe need to be cleverer here if there is no match by searching for best compromise group.
        pointGroup                  = self.pointGroups.get ( operationKey, None )
        characterSymmetryOperations = None
        try:
            characterSymmetryOperations = self.SetUpIrreducibleRepresentationCalculation ( nodeGroups, coordinates3, symmetryOperations, pointGroup )
        except:
            if doCharacterSymmetryOperations: raise

        # . Finish up.
        return { "Character Symmetry Operations" : characterSymmetryOperations ,
                 "Coordinates3"                  : coordinates3                ,
                 "Node Groups"                   : nodeGroups                  ,
                 "Point Group"                   : pointGroup                  ,
                 "Symmetry Operations"           : symmetryOperations          }

    #===============================================================================================================================
    def PostulateHigherOrderRotationAxes ( self, nodeGroups, coordinates3, log = None ):
        """Analyze group connectivity by distance to determine possible proper and improper rotation axes.

        Only axes of order >= 3 are searched for. In addition, composite symmetry operations are omitted.
        Thus, for example all S(2n+1) can be skipped as they automatically imply a C(2n+1) axis with a
        perpendicular mirror plane.
        """

        # . Initialization.
        # . Possible symmetry operations.
        improperRotations = {}
        properRotations   = {}

        # . Process each set of groups in turn.
        for ( node, groups ) in nodeGroups.items ( ):
            for group in groups:

                # . Skip small groups.
                if len ( group ) > 2:

                    # . Get the all internode distances within the group.
                    distances = []
                    for i in range ( 1, len ( group ) ):
                        iNode = group[i]
                        for j in range ( i ):
                            jNode = group[j]
                            distances.append ( ( coordinates3.Distance ( iNode, jNode ), ( iNode, jNode ) ) )

                    # . Separate distances by value.
                    distances.sort ( )
                    ( oldr, nodes ) = distances[0]
                    connectionSets = [ [ oldr, [ nodes ] ] ]
                    if len ( distances ) > 1:
                        for ( r, nodes ) in distances[1:]:
                            if math.fabs ( r - oldr ) <= self.distanceTolerance:
                                connectionSets[-1][1].append ( nodes )
                            else:
                                if len ( connectionSets[-1][1] ) < 2: connectionSets.pop ( -1 )
                                connectionSets.append ( [ r, [ nodes ] ] )
                                oldr = r
                        if len ( connectionSets[-1][1] ) < 2: connectionSets.pop ( -1 )

                    # . Process connections.
                    for ( r, connections ) in connectionSets:

                        # . Look for possible S4s.
                        # . This search involves looking for perpendicular connections equidistant from the origin but on opposite sides.
                        if len ( connections ) > 1:

                            # . Get data for the connections.
                            connectionData = []
                            for ( iNode, jNode ) in connections:
                                iAxis = Vector3.Null ( )
                                for i in range ( 3 ): iAxis[i] = coordinates3[iNode,i] + coordinates3[jNode,i]
                                distance = iAxis.Norm2 ( )
                                # . Centers at the origin do not contribute to S4 operations.
                                if distance > self.distanceTolerance:
                                    iAxis.Scale ( 1.0 / distance )
                                    iVector = Vector3.Null ( )
                                    for i in range ( 3 ): iVector[i] = coordinates3[iNode,i] - coordinates3[jNode,i]
                                    iVector.Normalize ( )
                                    connectionData.append ( ( distance, iAxis, iVector ) )

                            # . Check for appropriate pairs.
                            for i in range ( 1, len ( connectionData ) ):
                                ( iDistance, iAxis, iVector ) = connectionData[i]
                                for j  in range ( i ):
                                    ( jDistance, jAxis, jVector ) = connectionData[j]
                                    if ( math.fabs ( iDistance - jDistance     ) <= self.distanceTolerance            ) and \
                                       ( math.fabs ( iAxis.Dot ( jAxis ) + 1.0 ) <= self.cosineParallelTolerance      ) and \
                                       ( math.fabs ( iVector.Dot ( jVector )   ) <= self.cosinePerpendicularTolerance ):
                                        axes     = improperRotations.get ( 4, [] )
                                        isUnique = True
                                        for axis in axes:
                                            if ( math.fabs ( math.fabs ( iAxis.Dot ( axis ) ) - 1.0 ) <= self.cosineParallelTolerance ):
                                                isUnique = False
                                                break
                                        if isUnique:
                                            axes.append ( iAxis )
                                            improperRotations[4] = axes

                        # . Higher rotations (n > 2).
                        if len ( connections ) > 2:

                            # . Find all planar rings.
                            planarRings = self.FindPlanarRings ( connections, coordinates3 )

                            # . Proper rotations.
                            for ( order, rings ) in planarRings.items ( ):
                                for ( iRing, iDistance, iAxis, iCenter ) in rings:
                                    axes     = properRotations.get ( order, [] )
                                    isUnique = True
                                    for axis in axes:
                                        if ( math.fabs ( math.fabs ( iAxis.Dot ( axis ) ) - 1.0 ) <= self.cosineParallelTolerance ):
                                            isUnique = False
                                            break
                                    if isUnique:
                                        axes.append ( iAxis )
                                        properRotations[order] = axes

                            # . Improper rotations.
                            # . Find aligned pairs of parallel rings of the same size equidistant from the origin but on opposite sides.
                            for ( order, rings ) in planarRings.items ( ):
                                for i in range ( 1, len ( rings ) ):
                                    ( iRing, iDistance, iAxis, iCenter ) = rings[i]
                                    if iDistance > self.distanceTolerance:
                                        for j in range ( i ):
                                            ( jRing, jDistance, jAxis, jCenter ) = rings[j]
                                            if ( math.fabs ( iDistance - jDistance     ) <= self.distanceTolerance ) and \
                                               ( math.fabs ( iAxis.Dot ( jAxis ) + 1.0 ) <= self.cosineParallelTolerance   ):
                                               ( areAligned, improperOrder ) = self.CheckForAlignment ( order, iRing, iCenter, jRing, jCenter, coordinates3 )
                                               if areAligned:
                                                   axes     = improperRotations.get ( improperOrder, [] )
                                                   isUnique = True
                                                   for axis in axes:
                                                       if ( math.fabs ( math.fabs ( iAxis.Dot ( axis ) ) - 1.0 ) <= self.cosineParallelTolerance ):
                                                           isUnique = False
                                                           break
                                                   if isUnique:
                                                       axes.append ( iAxis )
                                                       improperRotations[improperOrder] = axes

        # . Print results.
        if LogFileActive ( log ):
            for ( tag, rotations ) in ( ( "Proper", properRotations ), ( "Improper", improperRotations ) ):
                if len ( rotations ) <= 0: log.Paragraph ( "No candidate " + tag.lower ( ) + " rotation axes found." )
                else:
                    table = log.GetTable ( columns = [ 20, 20 ] )
                    table.Start  ( )
                    table.Title  ( "Candidate {:s} Rotation Axes".format ( tag ) )
                    for order in sorted ( rotations.keys ( ) ):
                        table.Entry ( "{:d}".format (                 order    ) )
                        table.Entry ( "{:d}".format ( len ( rotations[order] ) ) )
                    table.Stop ( )

        # . Finish up.
        return ( properRotations, improperRotations )

    #===============================================================================================================================
    def PostulateSymmetryEquivalentNodes ( self, nodeTypes, coordinates3, log = None ):
        """Group nodes by type and distance to origin."""
        nodeData   = defaultdict ( list )
        nodeGroups = {}
        for ( i, nodeType ) in enumerate ( nodeTypes ):
            r = 0.0
            for j in range ( 3 ): r += coordinates3[i,j]**2
            r = math.sqrt ( r )
            nodeData[nodeType].append ( ( r, i ) )
        for node in sorted ( nodeData.keys ( ) ):
            data = nodeData[node]
            data.sort ( )
            ( oldR, i ) = data[0]
            groups = [ [ i ] ]
            if len ( data ) > 1:
                for ( r, i ) in data[1:]:
                    if math.fabs ( r - oldR ) <= self.distanceTolerance:
                        groups[-1].append ( i )
                    else:
                        groups.append ( [ i ] )
                        oldR = r
            nodeGroups[node] = groups
        if LogFileActive ( log ):
            n     = max ( [ len ( nodeGroups[node] ) for node in nodeData.keys ( ) ] )
            table = log.GetTable ( columns = [ 10, 5, 3 ] + n * [ 5 ] )
            table.Start  ( )
            table.Title  ( "Node Groups" )
            for node in sorted ( nodeData.keys ( ) ):
                groups = nodeGroups[node]
                table.Entry ( "{!r}".format ( node           ) )
                table.Entry ( "{:d}".format ( len ( groups ) ) )
                table.Entry ( ":" )
                for group in groups: table.Entry ( "{:d}".format ( len ( group ) ) )
                table.EndRow ( )
            table.Stop ( )
        return nodeGroups

    #===============================================================================================================================
    def ProcessHigherOrderRotationAxes ( self, nodeGroups, coordinates3, properRotations, improperRotations, symmetryOperations ):
        """Process postulated higher-order rotation axes.

        All non-unique operations of lower order are excluded.

        Thus, a S2m is redundant if there is a higher S2n axis such that (2m) * p = (2n) where p is odd.

        A Cm is redundant if there is a higher S2n axis such that m * p = (2n) where p is even.

        A Cm is also redundant if there is a higher Cn such that m * p = n where p is any integer.
        """

        # . Find the orders of the operations.
        improperOrders = list ( improperRotations.keys ( ) )
        improperOrders.sort ( reverse = True )
        properOrders   = list ( properRotations.keys   ( ) )
        properOrders.sort   ( reverse = True )

        # . Maximum orders.
        maximumImproperOrder = 0
        maximumProperOrder   = 0
        if len ( improperOrders ) > 0: maximumImproperOrder = improperOrders[0]
        if len (   properOrders ) > 0: maximumProperOrder   =   properOrders[0]

        # . Improper rotations.
        for order in improperOrders:
            axes       = improperRotations[order]
            operations = []

            # . Loop over axes.
            for axis in axes:
                # . Check whether to test.
                toTest = True
                toTry  = []
                i      = 3
                while ( i * order <= maximumImproperOrder ):
                    toTry.extend ( symmetryOperations.get ( "S{:d}".format ( i * order ), [] ) )
                    i += 2
                for other in toTry:
                    if ( math.fabs ( math.fabs ( other.axis.Dot ( axis ) ) - 1.0 ) <= self.cosineParallelTolerance ):
                        toTest = False
                        break
                if toTest:
                    operation = ImproperRotation.WithOptions ( axis = Clone ( axis ), order = order )
                    if operation.EstablishSymmetryRelatedPairs ( nodeGroups, coordinates3, self.distanceTolerance ):
                        operations.append ( operation )
            if len ( operations ) > 0: symmetryOperations["S{:d}".format ( order )] = operations

        # . Proper rotations.
        for order in properOrders:
            axes       = properRotations[order]
            operations = []

            # . Loop over axes.
            for axis in axes:
                # . Check whether to test.
                toTest = True
                toTry  = []
                i = 2
                while ( i * order <= maximumImproperOrder ):
                    toTry.extend ( symmetryOperations.get ( "S{:d}".format ( i * order ), [] ) )
                    i += 2
                i = 2
                while ( i * order <= maximumProperOrder ):
                    toTry.extend ( symmetryOperations.get ( "C{:d}".format ( i * order ), [] ) )
                    i += 1
                for other in toTry:
                    if ( math.fabs ( math.fabs ( other.axis.Dot ( axis ) ) - 1.0 ) <= self.cosineParallelTolerance ):
                        toTest = False
                        break
                # . Test.
                if toTest:
                    operation = ProperRotation.WithOptions ( axis = Clone ( axis ), order = order )
                    if operation.EstablishSymmetryRelatedPairs ( nodeGroups, coordinates3, self.distanceTolerance ):
                        operations.append ( operation )
            if len ( operations ) > 0: symmetryOperations["C{:d}".format ( order )] = operations

    #===============================================================================================================================
    # . A fudge until something better comes along.
    def ProcessIndeterminateC2Axes ( self, nodeGroups, coordinates3, symmetryOperations ):
        """Process indeterminate C2 axes."""
        # . Only treat if there are more than one mirror plane.
        sigmas = symmetryOperations.get ( "sigma", [] )
        if len ( sigmas ) > 1:

            # . Gather all rotations (C2n and S4n) which can produce a C2.
            existing = []
            for ( key, operations ) in symmetryOperations.items ( ):
                if ( key.startswith ( "C" ) and ( operations[0].order % 2 == 0 ) ) or \
                   ( key.startswith ( "S" ) and ( operations[0].order % 4 == 0 ) ):
                    for operation in operations:
                        existing.append ( operation.axis )

            # . Find possible C2 axes from pairs of perpendicular mirror planes.
            axes = []
            for i in range ( 1, len ( sigmas ) ):
                iNormal = sigmas[i].normal
                for j in range ( i ):
                    jNormal = sigmas[j].normal
                    # . Perpendicular normals.
                    if ( math.fabs ( iNormal.Dot ( jNormal ) ) <= self.cosinePerpendicularTolerance ):
                        axis = Clone ( iNormal )
                        axis.Cross ( jNormal )
                        axis.Normalize ( )
                        # . Check for uniqueness.
                        isUnique = True
                        for other in existing:
                            if ( math.fabs ( math.fabs ( other.Dot ( axis ) ) - 1.0 ) <= self.cosineParallelTolerance ):
                                isUnique = False
                                break
                        if isUnique:
                            axes.append     ( axis )
                            existing.append ( axis )

            # . Process the additional axes.
            if len ( axes ) > 0:
                rotations = symmetryOperations.get ( "C2", [] )
                for axis in axes:
                    operation = ProperRotation.WithOptions ( axis = Clone ( axis ), order = 2 )
                    if operation.EstablishSymmetryRelatedPairs ( nodeGroups, coordinates3, self.distanceTolerance ):
                        rotations.append ( operation )
                if len ( rotations ) > 0: symmetryOperations["C2"] = rotations

    #===============================================================================================================================
    def ProcessMomentsOfInertia ( self, nodeGroups, coordinates3, weights, symmetryOperations ):
        """Check the moments of inertia for basic information about the graph."""

        # . Inertia tolerance.
        inertiaTolerance = 2.0 * sum ( weights ) * self.distanceTolerance**2

        # . Find the number of "zero" adjusted moments of inertia.
        ( moments, axes ) = coordinates3.MomentsOfInertia ( weights = weights )
        mr2    = 0.5 * sum ( moments )
        zeroes = []
        for i in range ( 3 ):
            if ( math.fabs ( mr2 - moments[i] ) <= inertiaTolerance ): zeroes.append ( i )
        nZeroes = len ( zeroes )

        # . Initialization.
        axis = Vector3.Null ( )

        # . Check for planarity.
        if nZeroes > 0:
            for i in zeroes:
                for j in range ( 3 ): axis[j] = axes[j,i]
                operation = Reflection.WithOptions ( normal = Clone ( axis ), order = 1 )
                if operation.EstablishSymmetryRelatedPairs ( nodeGroups, coordinates3, self.distanceTolerance ):
                    symmetryOperations["sigma"] = [ operation ]
                    break

        # . Check for linearity.
        if nZeroes > 1:
            for i in range ( 3 ):
                if ( ( i + 1 ) % 3 in zeroes ) and ( ( i + 2 ) % 3 in zeroes ):
                    for j in range ( 3 ): axis[j] = axes[j,i]
                    operation = ProperRotation.WithOptions ( axis = Clone ( axis ), order = 1 )
                    if operation.EstablishSymmetryRelatedPairs ( nodeGroups, coordinates3, self.distanceTolerance ):
                        symmetryOperations["CInfinity"] = [ operation ]
                        break

    #===============================================================================================================================
    def ProcessOrderTwoOperations ( self, nodeGroups, coordinates3, symmetryOperations ):
        """Find order-2 operations (C2 rotation, inversion, reflection)."""

        # . Initialization.
        hasIndeterminateC2Axes = False
        axis          = Vector3.Null ( )
        vij           = Vector3.Null ( )
        moleculePlane = None
        rotations     = []
        reflections   = []
        if "sigma" in symmetryOperations:
            moleculePlane = symmetryOperations["sigma"][0]
            reflections.append ( moleculePlane )

        # . Gather all rotations (C2n and S4n) which can produce a C2.
        existing = []
        for ( key, operations ) in symmetryOperations.items ( ):
            if ( key.startswith ( "C" ) and ( operations[0].order % 2 == 0 ) ) or \
               ( key.startswith ( "S" ) and ( operations[0].order % 4 == 0 ) ):
                existing.extend ( operations )

        # . Inversion.
        operation = Inversion.WithDefaults ( )
        if operation.EstablishSymmetryRelatedPairs ( nodeGroups, coordinates3, self.distanceTolerance ):
            symmetryOperations["i"] = [ operation ]

        # . C2 rotations and reflections.
        # . Loop over all pairs of possibly symmetry-related nodes.
        for ( node, groups ) in nodeGroups.items ( ):
            for group in groups:
                n = len ( group )
                for i in range ( 1, n ):
                    inode = group[i]
                    for j in range ( i ):
                        jnode = group[j]
                        # . C2.
                        maximumValue       = 0.0
                        tryC2              = True
                        usedMolecularPlane = False
                        for s in range ( 3 ):
                            axis[s]      = coordinates3[inode,s] + coordinates3[jnode,s]
                            maximumValue = max ( maximumValue, math.fabs ( axis[s] ) )
                        if ( maximumValue < self.distanceTolerance ):
                            if ( moleculePlane is not None ):
# . Check.
                                moleculePlane.normal.CopyTo ( axis )
                                usedMolecularPlane = True
                            else:
                                hasIndeterminateC2Axes = True
                                tryC2 = False
                        if tryC2:
                            axis.Normalize ( )
                            isUnique = True
                            for other in ( existing + rotations ):
                                if ( math.fabs ( math.fabs ( other.axis.Dot ( axis ) ) - 1.0 ) <= self.cosineParallelTolerance ):
                                    isUnique = False
                                    break
                            if isUnique:
                                operation = ProperRotation.WithOptions ( axis = Clone ( axis ), order = 2 )
                                if operation.EstablishSymmetryRelatedPairs ( nodeGroups, coordinates3, self.distanceTolerance ):
                                    rotations.append ( operation )
                        # . Additional C2 in cases where there is a molecular plane and the (ij) mid-point passes through the origin.
                        if usedMolecularPlane:
                            for s in range ( 3 ): vij[s] = coordinates3[inode,s] - coordinates3[jnode,s]
                            normal = moleculePlane.normal
                            axis[0]  = normal[1]*vij[2] - vij[1]*normal[2]
                            axis[1]  = normal[2]*vij[0] - vij[2]*normal[0]
                            axis[2]  = normal[0]*vij[1] - vij[0]*normal[1]
                            axis.Normalize ( )
                            isUnique = True
                            for other in ( existing + rotations ):
                                if ( math.fabs ( math.fabs ( other.axis.Dot ( axis ) ) - 1.0 ) <= self.cosineParallelTolerance ):
                                    isUnique = False
                                    break
                            if isUnique:
                                operation = ProperRotation.WithOptions ( axis = Clone ( axis ), order = 2 )
                                if operation.EstablishSymmetryRelatedPairs ( nodeGroups, coordinates3, self.distanceTolerance ):
                                    rotations.append ( operation )
                        # . Reflection.
                        for s in range ( 3 ): axis[s] = coordinates3[inode,s] - coordinates3[jnode,s]
                        axis.Normalize ( )
                        isUnique = True
                        for other in reflections:
                            if ( math.fabs ( math.fabs ( other.normal.Dot ( axis ) ) - 1.0 ) <= self.cosineParallelTolerance ):
                                isUnique = False
                                break
                        if isUnique:
                            operation = Reflection.WithOptions ( normal = Clone ( axis ), order = 2 )
                            if operation.EstablishSymmetryRelatedPairs ( nodeGroups, coordinates3, self.distanceTolerance ):
                                reflections.append ( operation )

        # . Finish up.
        if len ( rotations   ) > 0: symmetryOperations["C2"   ] = rotations
        if len ( reflections ) > 0: symmetryOperations["sigma"] = reflections
        return hasIndeterminateC2Axes

    #===============================================================================================================================
    def RemoveRedundantS4s ( self, symmetryOperations ):
        """Remove S4s that can be formed by combination of a C4 and a reflection.

        This is not done previously as S4s are found before C4s and sigmas.
        """
        C4s    = symmetryOperations.get ( "C4"   , None )
        S4s    = symmetryOperations.get ( "S4"   , None )
        sigmas = symmetryOperations.get ( "sigma", None )
        if ( C4s is not None ) and ( S4s is not None ) and ( sigmas is not None ):
            newS4s = []
            for S4 in S4s:
                axis      = S4.axis
                toInclude = True
                for C4 in C4s:
                    if ( math.fabs ( math.fabs ( C4.axis.Dot ( axis ) ) - 1.0 ) <= self.cosineParallelTolerance ):
                        for sigma in sigmas:
                            if ( math.fabs ( math.fabs ( sigma.normal.Dot ( axis ) ) - 1.0 ) <= self.cosineParallelTolerance ):
                                toInclude = False
                                break
                        if not toInclude: break
                if toInclude: newS4s.append ( S4 )
            if len ( newS4s ) > 0:     symmetryOperations["S4"] = newS4s
            else:                  del symmetryOperations["S4"]

    #===============================================================================================================================
    def SetUpIrreducibleRepresentationCalculation ( self, nodeGroups, coordinates3, symmetryOperations, pointGroup ):
        """Set up the information necessary for determination of irreducible representations."""
        # . Initialization.
        characterSymmetryOperations = {}
        if ( coordinates3 is None ) or ( pointGroup is None ) or ( symmetryOperations is None ):
            raise SymmetryError ( "Insufficient data for determining irreducible representations." )
        # . Build CInfinity operations if required.
        if ( pointGroup.cInfinityRotations is not None ) and ( len ( pointGroup.cInfinityRotations ) > 0 ):
            cInfinity = symmetryOperations.get ( "CInfinity", None )
            if cInfinity is None: raise SymmetryError ( "Unable to find CInfinity rotation for point group: " + pointGroup.label + "." )
            for rotation in pointGroup.cInfinityRotations:
                operation = ProperRotation.WithOptions ( axis = Clone ( cInfinity[0].axis ), order = int ( rotation[1:] ) )
                if operation.EstablishSymmetryRelatedPairs ( nodeGroups, coordinates3, self.distanceTolerance ):
                    characterSymmetryOperations[rotation] = operation
                else: raise SymmetryError ( "CInfinity rotation " + rotation + " for point group " + pointGroup.label + " is not a symmetry operation." )
        # . Process existing operations.
        labels = pointGroup.characterSymmetryOperations
        # . Get sigma-h.
        sigmaH = None
        sigmas = symmetryOperations.get ( "sigma", [] )
        if ( pointGroup.principalAxisOperation is not None ) and ( ( "sigma-h" in labels ) or ( "sigma-m" in labels ) ):
            principalAxisOperation = symmetryOperations[pointGroup.principalAxisOperation][0]
            for sigma in sigmas:
                if ( math.fabs ( math.fabs ( sigma.normal.Dot ( principalAxisOperation.axis ) ) - 1.0 ) <= self.cosineParallelTolerance ):
                    sigmaH = sigma
                    break
            if ( "sigma-h" in labels ):
                if sigmaH is None: raise SymmetryError ( "Unable to find \"sigma-h\" operation for point group: " + pointGroup.label + "." )
                else: characterSymmetryOperations["sigma-h"] = sigmaH
        # . Get sigma-m.
        if ( "sigma-m" in labels ):
            sigmas = list ( sigmas )
            if sigmaH is not None: sigmas.remove ( sigmaH )
            if len ( sigmas ) > 0:
                sigmaM = sigmas[0]
                for sigma in sigmas[1:]:
                    if sigma.selfMappings > sigmaM.selfMappings: sigmaM = sigma
                characterSymmetryOperations["sigma-m"] = sigmaM
            else: raise SymmetryError ( "Unable to find \"sigma-m\" operation for point group: " + pointGroup.label + "." )
        # . Get C2-m.
        if ( "C2-m" in labels ):
            C2s = symmetryOperations.get ( "C2", [] )
            if len ( C2s ) > 0:
                C2M = C2s[0]
                for C2 in C2s[1:]:
                    if C2.selfMappings > C2M.selfMappings: C2M = C2
                characterSymmetryOperations["C2-m"] = C2M
            else: raise SymmetryError ( "Unable to find \"C2-m\" operation for point group: " + pointGroup.label + "." )
        # . Existing operations.
        for operation in pointGroup.characterSymmetryOperations:
            if operation not in characterSymmetryOperations.keys ( ):
                inputOperations = symmetryOperations.get ( operation, None )
                if inputOperations is None: raise SymmetryError ( "Unable to find required operation \"" + operation + "\" for point group: " + pointGroup.label + "." )
                characterSymmetryOperations[operation] = inputOperations[0]
        # . Finish up.
        return [ characterSymmetryOperations[key] for key in pointGroup.characterSymmetryOperations ]

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def Find3DGraphPointGroup ( nodeTypes, coordinates3, **options ):
    """Find the point group of a 3D-graph."""
    doCharacterSymmetryOperations = options.pop ( "doCharacterSymmetryOperations" , True    )
    log                           = options.pop ( "log"                           , logFile )
    weights                       = options.pop ( "weights"                       , None    )
    finder                        = PointGroupFinder.WithOptions ( **options )
    return finder.Find3DGraphPointGroup ( nodeTypes, coordinates3, doCharacterSymmetryOperations = doCharacterSymmetryOperations, log = log, weights = weights )

# . This method requires a set of item values (real) and a function for determining the characters of each item under a given symmetry operation.
def IdentifyIrreducibleRepresentations ( results             ,
                                         itemValues          ,
                                         GetItemCharacters   ,
                                         degeneracyTolerance ,
                                         maximumIRs = 1      ):
    """Identify irreducible representations for a set of quantities pertaining to an already processed 3D-graph."""
    # . Check results.
    characterSymmetryOperations = results.get ( "Character Symmetry Operations", None )
    pointGroup                  = results.get ( "Point Group"                  , None )
    if ( characterSymmetryOperations is None ) or ( pointGroup is None ):
        raise SymmetryError ( "Missing point group and/or character symmetry operations for identification of irreducible representations." )
    # . Initialization.
    characters                 = None
    irreducibleRepresentations = []
    if len ( itemValues ) > 0:
        # . Get characters for each state under each symmetry operation.
        characters = Array.WithExtents ( len ( itemValues ), len ( characterSymmetryOperations ) )
        characters.Set ( 1.0 )
        # . Loop over symmetry operations.
        for ( j, operation ) in enumerate ( characterSymmetryOperations ):
            if operation.label != "E":
                localCharacters = GetItemCharacters ( operation )
                for ( i, character ) in enumerate ( localCharacters ): characters[i,j] = character
        # . Identify symmetries.
        irreducibleRepresentations = pointGroup.IdentifyIrreducibleRepresentations ( characters, itemValues, degeneracyTolerance, maximumIRs = maximumIRs )
    # . Finish up.
    return ( irreducibleRepresentations, characters )

# . Utility function for printing results.
def PrintIrreducibleRepresentations ( results               ,
                                      itemValues            ,
                                      iRs                   ,
                                      characters            ,
                                      itemFormat = "{:.3f}" ,
                                      itemName   = "Item"   , 
                                      log        = logFile  ,
                                      title      = None     ,
                                      valueName  = "Value"  ):
    """Print a summary of the results of finding the irreducible representations."""
    if LogFileActive ( log ):
        characterSymmetryOperations = results.get ( "Character Symmetry Operations", None )
        pointGroup                  = results.get ( "Point Group"                  , None )
        if ( characters is None ) or ( characterSymmetryOperations is None ) or ( pointGroup is None ):
            logFile.Paragraph ( "Point group information missing." )
        else:
            ( n, c ) = characters.shape
            columns  = [ max ( 6, len ( itemName ) + 2 ) ] + ( c + 1 ) * [ 12 ] + [ 8 ]
            if title is None: title = "{:s} {:s} Characters and IRs".format ( itemName, pointGroup.label )
            table = logFile.GetTable ( columns = columns )
            table.Start  ( )
            table.Title   ( title     )
            table.Heading ( itemName  )
            table.Heading ( valueName )
            for operation in characterSymmetryOperations:
                table.Heading ( operation.label )
            table.Heading ( "IR" )
            for i in range ( n ):
                table.Entry ( "{:d}".format ( i ), align = Align.Center )
                table.Entry ( itemFormat.format ( itemValues[i] ) )
                for j in range ( c ):
                    table.Entry ( "{:.3f}".format ( characters[i,j] ) )
                table.Entry ( "{:s}".format ( iRs[i] ) )
            table.Stop ( )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass

