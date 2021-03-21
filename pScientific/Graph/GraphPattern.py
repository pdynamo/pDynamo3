"""A graph pattern is a set of node and edge patterns that can be used to match various parts of a graph."""

import copy

from  pCore       import AttributableObject
from .GraphStatus import GraphError

#
# . Notes:
#
#   - Can generalize matching to allow ranges or sets of values to match, rather than single values.
#   - Can have a method that matches edges principally rather than nodes (including edge results).
#

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . The character indicating that a pattern attribute is not defined.
_Undefined = "."

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class _MatchObject:
    """A local class for holding match data."""

    def __init__ ( self ):
        """Constructor."""
        self.edges = {}
        self.nodes = {}

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class EdgePattern ( AttributableObject ):
    """A class to represent an edge pattern."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "nodeKey1" : -1 ,
                             "nodeKey2" : -1 } )

    @classmethod
    def FromList ( selfClass, attributeList, attributes, undefined ):
        """Constructor from list."""
        options = {}
        for ( key, value ) in zip ( attributes, attributeList ):
            if value != undefined: options[key] = value
        return selfClass.WithOptions ( **options )

    def Match ( self, edge ):
        """Match the pattern against an edge."""
        return True

    def Opposite ( self, nodeKey ):
        """Return the index of the other node."""
        if   self.nodeKey1 == nodeKey: return self.nodeKey2
        elif self.nodeKey2 == nodeKey: return self.nodeKey1
        else: return -1

    def ToList ( self, attributes, undefined ):
        """Get a list representation of the object."""
        return [  getattr ( self, attribute, undefined ) for attribute in attributes ]

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class NodePattern ( AttributableObject ):
    """A class to represent a node pattern."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "attributes" : dict } ) # . Is this really necessary?

    @classmethod
    def FromList ( selfClass, attributeList, attributes, undefined ):
        """Constructor from list."""
        options = {}
        for ( key, value ) in zip ( attributes, attributeList ):
            if value != undefined: options[key] = value
        return selfClass.WithOptions ( **options )

    def Match ( self, node ):
        """Match the pattern against a node."""
        isMatched = True
        for ( key, pAttribute ) in self.attributes.items ( ):
            if pAttribute is not None:
                try:    attribute = getattr ( node, key )
                except: attribute = None
                isMatched = ( attribute == pAttribute )
                if not isMatched: break
        return isMatched

    def SetOptions ( self, **options ):
        """Set options."""
        self.attributes.update ( { key : value for ( key, value ) in options.items ( ) } )

    def ToList ( self, attributes, undefined ):
        """Get a list representation of the object."""
        return [  self.attributes.get ( attribute, undefined ) for attribute in attributes ]

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class GraphPattern ( AttributableObject ):
    """A class to represent a graph pattern."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "edgeConnections" : None ,
                             "edgePatterns"    : None ,
                             "isSupporting"    : None ,
                             "label"           : None ,
                             "nodeConnections" : None ,
                             "nodePatterns"    : None ,
                             "nodeResults"     : None ,
                             "upperBound"      :   -1 } )
    _edgeClass    = EdgePattern
    _edgeTag      = "Edge"
    _nodeClass    = NodePattern
    _nodeTag      = "Node"

    def ClearRepresentations ( self ):
        """Clear all representations."""
        self.nodeConnections = None
        self.edgeConnections = None
        self.upperBound      =   -1

    def FindAllMatches ( self, graph, selection = None ):
        """Find all occurrences of |self| that occur in |graph|.

        This method uses a depth-first method to find all matches. Uniqueness
        is assured only at the end when incorporating matches into the final
        list and not when generating them. The algorithm in this method could
        undoubtedly be made more efficient.
        """
        # . Initialization.
        matches = []
        # . Treat the selection.
        if selection is None: selection = graph.nodes
        # . Simple checks to avoid too much work.
        if ( len ( graph.nodes ) >= len ( self.nodePatterns ) ) and ( len ( selection ) > 0 ):
            # . Find matches for the first node pattern.
            firstNode    = self.nodePatterns[0]
            isSingleNode = ( len ( self.nodePatterns ) == 1 )
            for iNode in selection:
                if firstNode.Match ( iNode ):
                    # . Single node patterns.
                    if isSingleNode: matches.append ( [ iNode ] )
                    # . Multiple node patterns.
                    else:
                        mObject          = _MatchObject ( )
                        mObject.nodes[0] = iNode
                        mObjects         = self.MatchBranches ( iNode, 0, mObject, graph, selection )
                        for mObject in mObjects:
                            if ( len ( mObject.nodes ) == len ( self.nodePatterns ) ) and ( len ( mObject.edges ) == len ( self.edgePatterns ) ):
                                match = []
                                for i in range ( len ( mObject.nodes ) ): match.append ( mObject.nodes[i] )
                                if match not in matches: matches.append ( match )
        return matches

    @classmethod
    def FromMapping ( selfClass                   ,
                      mapping                     ,
                      nodeAttributes       = None ,
                      nodeResultAttributes = None ,
                      edgeAttributes       = None ):
        """Constructor from a mapping."""
        # . Basic construction.
        self         = selfClass ( )
        edgePatterns = mapping.get ( "{:s} Patterns".format ( self.__class__._edgeTag ), [] )
        nodePatterns = mapping.get ( "{:s} Patterns".format ( self.__class__._nodeTag ), [] )
        self.label   = mapping.get ( "Label", None )
        # . Class-specific processing.
        self.ProcessFromMappingAttributes ( nodeAttributes       ,
                                            nodePatterns         ,
                                            nodeResultAttributes ,
                                            edgeAttributes       ,
                                            edgePatterns         ,
                                            _Undefined           )
        # . Node and edge patterns.
        self.nodePatterns = [ self.__class__._nodeClass.FromList ( patternList, nodeAttributes, _Undefined ) for patternList in nodePatterns ]
        self.edgePatterns = [ self.__class__._edgeClass.FromList ( patternList, edgeAttributes, _Undefined ) for patternList in edgePatterns ]
        return self

    def GetConnections ( self, i ):
        """Get the connections for a node pattern."""
        if ( self.edgeConnections is not None ):
            try:    iNode = int ( i )
            except: raise TypeError ( "Node index must be an integer." )
            if ( iNode >= 0 ) and ( iNode < self.upperBound ):
                try:    return self.edgeConnections[iNode]
                except: return []
            else: raise IndexError ( "Connection index out of range." )
        else: raise AttributeError ( "Connection representation does not exist." )

    def MakeConnections ( self, upperBound = -1 ):
        """Make the connections representation of the edge container."""
        # . Make the representation.
        if ( self.nodeConnections is None ) or ( self.edgeConnections is None ):
            self.ClearRepresentations ( )
            # . Initialization.
            maximumIndex = len ( self.nodePatterns )
            self.nodeConnections = {}
            self.edgeConnections = {}
            for ( b, edge ) in enumerate ( self.edgePatterns ):
                i = edge.nodeKey1
                j = edge.nodeKey2
                if i not in self.nodeConnections: self.nodeConnections[i] = []
                if j not in self.nodeConnections: self.nodeConnections[j] = []
                if i not in self.edgeConnections: self.edgeConnections[i] = []
                if j not in self.edgeConnections: self.edgeConnections[j] = []
                self.nodeConnections[i].append ( j )
                self.nodeConnections[j].append ( i )
                self.edgeConnections[i].append ( b )
                self.edgeConnections[j].append ( b )
                maximumIndex = max ( maximumIndex, i, j )
            # . Sorting.
            for value in self.nodeConnections.values ( ): value.sort ( )
            for value in self.edgeConnections.values ( ): value.sort ( )
            # . Set the upper bound of the representation.
            self.upperBound = max ( maximumIndex, upperBound )
        # . Check the upperBound.
        else:
            if upperBound > self.upperBound: self.upperBound = upperBound

    def MatchBranches ( self, cNode, pNode, mObject, graph, selection ):
        """Match nodes connected to a root node."""
        # . Get the edges for the current nodes.
        cEdges = graph.adjacentEdges[cNode]
        pEdges = copy.copy ( self.GetConnections ( pNode ) )
        # . Remove matched edges from the data.
        for pEdge in pEdges:
            if pEdge in mObject.edges: pEdges.remove ( pEdge )
        # . There are no edges to match so return the current data as is.
        if len ( pEdges ) == 0:
            mObjects = [ mObject ]
        # . Match all edge patterns and descending branches in order.
        else:
            # . Initialization.
            mObjects = []
            # . Always do full double loop without exiting so that all pattern/real edge permutations are explored.
            for pEdge in pEdges:
                qNode         = self.edgePatterns[pEdge].Opposite ( pNode )
                isQMatched    = ( qNode in mObject.nodes )
                isQSupporting = self.isSupporting[qNode]
                # . Loop over edges for cNode.
                for cEdge in cEdges:
                    dNode = cEdge.Opposite ( cNode )
                    # . dNode must be in selection or the match node must be supporting.
                    if ( dNode in selection ) or isQSupporting:
                        isDMatched = ( dNode in mObject.nodes.values ( ) )
                        # . Check for an edge match.
                        if self.edgePatterns[pEdge].Match ( cEdge ):
                            # . Check for a node match.
                            if self.nodePatterns[qNode].Match ( dNode ):
                                # . Copy and update the match.
                                nObject = copy.copy ( mObject )
                                nObject.edges[pEdge] = cEdge
                                # . A previously matched pattern node.
                                if isQMatched:
                                    # . The match is OK if the real node corresponds.
                                    if mObject.nodes[qNode] == dNode: mObjects.append ( nObject )
                                # . A previously unmatched real node.
                                elif not isDMatched:
                                    nObject.nodes[qNode] = dNode
                                    # . Branches arising from the new node.
                                    tObjects = self.MatchBranches ( dNode, qNode, nObject, graph, selection )
                                    # . Branches arising from unmatched edges of the current node (if there are any).
                                    if len ( tObjects ) > 0:
                                        if tObjects[0] is nObject:
                                            mObjects.append ( nObject )
                                        else:
                                            for tObject in tObjects:
                                                mObjects.extend ( self.MatchBranches ( cNode, pNode, tObject, graph, selection ) )
        # . Finish up.
        return mObjects

    def ProcessFromMappingAttributes ( self                 ,
                                       nodeAttributes       ,
                                       nodePatterns         ,
                                       nodeResultAttributes ,
                                       edgeAttributes       ,
                                       edgePatterns         ,
                                       undefined            ):
        """Process class-specific from mapping attributes."""
        # . Node keys.
        try:
            k    = nodeAttributes.index ( "key" )
            keys = {}
            for ( i, pattern ) in enumerate ( nodePatterns ):
                key = pattern[k]
                if key != undefined:
                    keys[key]  = i
                    pattern[k] = undefined
            tag = self.__class__._nodeTag.lower ( )
            k1  = edgeAttributes.index ( "{:s}Key1".format ( tag ) )
            k2  = edgeAttributes.index ( "{:s}Key2".format ( tag ) )
            for pattern in edgePatterns:
                key1 = pattern[k1]
                key2 = pattern[k2]
                pattern[k1] = keys[key1]
                pattern[k2] = keys[key2]
        except:
            pass
        # . Is supporting - specific and then default values.
        self.isSupporting = self.ProcessFromMappingNodeArray ( "isSupporting" ,
                                                               nodeAttributes ,
                                                               nodePatterns   ,
                                                               undefined      )
        if self.isSupporting is None:
            self.isSupporting = [ False for i in range ( len ( nodePatterns ) ) ]
        # . Node results.
        self.nodeResults = {}
        if nodeResultAttributes is not None:
            for attribute in nodeResultAttributes:
                self.nodeResults[attribute] = self.ProcessFromMappingNodeArray ( attribute      ,
                                                                                 nodeAttributes ,
                                                                                 nodePatterns   ,
                                                                                 undefined      )

    def ProcessFromMappingNodeArray ( self           ,
                                      attribute      ,
                                      nodeAttributes ,
                                      nodePatterns   ,
                                      undefined      ):
        """Process a from mapping node array."""
        values = None
        try:
            index  = nodeAttributes.index ( attribute )
            values = []
            for pattern in nodePatterns:
                value = pattern[index]
                if value != undefined:
                    values.append ( value )
                    pattern[index] = undefined
        except:
            pass
        return values

    def ProcessToMappingAttributes ( self, nodeAttributes       ,
                                           nodePatterns         ,
                                           nodeResultAttributes ,
                                           edgeAttributes       ,
                                           edgePatterns         ):
        """Process class-specific to mapping attributes."""
        # . Post-process node patterns.
        try   : k = nodeAttributes.index ( "key"          )
        except: k = -1
        try   : s = nodeAttributes.index ( "isSupporting" )
        except: s = -1
        for ( i, patternList ) in enumerate ( nodePatterns ):
            if k >= 0: patternList[k] = i
            if s >= 0: patternList[s] = self.isSupporting[i]
        # . Node results.
        if nodeResultAttributes is not None:
            for attribute in nodeResultAttributes:
                if attribute in nodeAttributes:
                    index   = nodeAttributes.index ( attribute )
                    results = self.nodeResults[attribute]
                    for ( i, patternList ) in enumerate ( nodePatterns ):
                        patternList[index] = results[i]
        # . Post-process edge patterns - none here.

    def ToMapping ( self                        ,
                    nodeAttributes       = None ,
                    nodeResultAttributes = None ,
                    edgeAttributes       = None ):
        """Get a mapping representation of the pattern."""
        mapping = {}
        # . Node patterns.
        nodePatterns = {}
        for ( i, pattern ) in enumerate ( self.nodePatterns ):
            patternList = pattern.ToList ( nodeAttributes, _Undefined )
            nodePatterns.append ( patternList )
        # . Edge patterns.
        edgePatterns = {}
        for pattern in self.edgePatterns:
            patternList = pattern.ToList ( edgeAttributes, _Undefined )
            edgePatterns.append ( patternList )
        # . Post-processing.
        self.ProcessToMappingAttributes ( nodeAttributes       ,
                                          nodePatterns         ,
                                          nodeResultAttributes ,
                                          edgeAttributes       ,
                                          edgePatterns         )
        # . Finish up.
        mapping["{:s} Patterns".format ( self.__class__._edgeTag )] = edgePatterns
        mapping["{:s} Patterns".format ( self.__class__._nodeTag )] = nodePatterns
        if self.label is not None: mapping["Label"] = self.label
        return mapping

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class GraphPatternContainer ( AttributableObject ):
    """A container for graph patterns."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "edgeAttributes"       : None ,
                             "edgeFields"           : None ,
                             "items"                : None ,
                             "label"                : None ,
                             "nodeAttributes"       : None ,
                             "nodeFields"           : None ,
                             "nodeResultAttributes" : None ,
                             "nodeResultFields"     : None ,
                             "properties"           : None ,
                             "rawItems"             : None ,
                             "termLabel"            : "Graph Pattern" } )
    _patternClass = GraphPattern

    #yaml_tag = "!GraphPatternContainer"

    def __getstate__ ( self ):
        """Get state as a mapping."""
        mapping = {}
        # . Keys.
        edgeTag = self.__class__._patternClass._edgeTag
        nodeTag = self.__class__._patternClass._nodeTag
        mapping["{:s} Fields".format        ( edgeTag )] = self.edgeFields
        mapping["{:s} Fields".format        ( nodeTag )] = self.nodeFields
        mapping["{:s} Result Fields".format ( nodeTag )] = self.nodeResultFields
        # . Patterns.
        patterns = []
        for pattern in self.rawItems:
            patterns.append ( pattern.ToMapping ( nodeAttributes       = self.nodeAttributes       , \
                                                  nodeResultAttributes = self.nodeResultAttributes , \
                                                  edgeAttributes       = self.edgeAttributes       ) )
        mapping["Patterns"] = patterns
        # . Other data.
        if self.label is not None: mapping["Label"] = self.label
        return mapping

    def __setstate__ ( self, mapping ):
        #"""Set state from a mapping."""
        # . There are assumed to be no errors!
        #try:
            # . Basic construction.
            edgeTag = self.__class__._patternClass._edgeTag
            nodeTag = self.__class__._patternClass._nodeTag
            self.__init__ ( )
            self.edgeFields       = mapping.pop ( "{:s} Fields".format        ( edgeTag ) , None )
            self.label            = mapping.pop ( "Label"                                 , None )
            self.nodeFields       = mapping.pop ( "{:s} Fields".format        ( nodeTag )        )
            self.nodeResultFields = mapping.pop ( "{:s} Result Fields".format ( nodeTag ) , None )
            patterns              = mapping.pop ( "Patterns"                                     )
            self.properties       = dict ( mapping )
            # . Convert keys to attributes.
            self.edgeAttributes       = self.FieldsToAttributes ( self.edgeFields       )
            self.nodeAttributes       = self.FieldsToAttributes ( self.nodeFields       )
            self.nodeResultAttributes = self.FieldsToAttributes ( self.nodeResultFields )
            if not set ( self.nodeResultAttributes ).issubset ( set ( self.nodeAttributes ) ): raise # . All result attributes must be in attributes.
            # . Create the raw items.
            self.rawItems = []
            for mapping in patterns:
                self.rawItems.append ( self.__class__._patternClass.FromMapping ( mapping                                          , \
                                                                                  nodeAttributes       = self.nodeAttributes       , \
                                                                                  nodeResultAttributes = self.nodeResultAttributes , \
                                                                                  edgeAttributes       = self.edgeAttributes       ) )
            # . Create the items that are to be used.
            self.ProcessRawItems ( )
#        except Exception as e:
#            print e[0]
        #except:
        #    raise GraphError ( "Error reconstructing instance of " + self.__class__.__name__ + "." )

    def FieldsToAttributes ( self, fields ):
        """Convert fields to attributes."""
        attributes = []
        if fields is not None:
            for field in fields:
                attribute = "".join ( field.split ( ) )
                attributes.append ( attribute[0:1].lower ( ) + attribute[1:] )
        return attributes

    def ProcessRawItems ( self ):
        """Process the raw items."""
        self.items = self.rawItems

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
