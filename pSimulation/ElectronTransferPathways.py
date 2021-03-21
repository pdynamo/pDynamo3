"""Classes to calculate electron transfer pathways using the pathways approach."""

import math

from pCore             import AttributableObject , \
                              logFile            , \
                              LogFileActive
from pScientific.Graph import Edge               , \
                              Graph              , \
                              Node               , \
                              Path

# . Completely unfinished - do not use!

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================

#===================================================================================================================================
# . Error.
#===================================================================================================================================
class PathwaysError ( Exception ):
    """A pathways error."""
    pass

#===================================================================================================================================
# . Edges.
#===================================================================================================================================
class PathwaysEdge ( Edge ):
    """A link between states."""

    _attributable = dict ( Edge._attributable )
    _attributable.update ( { "coupling"         :  0.0 ,
                             "couplingFunction" : None } )

    def CalculateCoupling ( self ):
        """Calculate the coupling and the weight."""
        # . The coupling function must return a number whose magnitude is less than or equal to one.
        self.coupling = self.couplingFunction ( self.node1, self.node2 )
        self.weight   = - math.log ( math.fabs ( self.coupling ) )

    @classmethod
    def FromCouplingFunctions ( selfClass, node1, node2, couplingFunctions, missingCouplingFunctions ):
        """Constructor given a coupling function dictionary."""
        label1 = node1.definition.label
        label2 = node2.definition.label
        key    = ( max ( label1, label2 ), min ( label1, label2 ) )
        couplingFunction = couplingFunctions.get ( key, None )
        if couplingFunction is not None:
            self = selfClass.WithNodes ( node1, node2, couplingFunction = couplingFunction )
        else:
            missingCouplingFunctions.add ( key )
            self = None
        return self

class ThroughBondLink ( PathwaysEdge ):
    """A through bond link."""
    pass

class ThroughSpaceLink ( PathwaysEdge ):
    """A through space link."""
    pass

#===================================================================================================================================
# . Nodes.
#===================================================================================================================================
class PathwaysNode ( Node ):
    """A state."""

    _attributable = dict ( Node._attributable )
    _attributable.update ( { "definition" : None ,
                             "x"          :  0.0 ,
                             "y"          :  0.0 ,
                             "z"          :  0.0 } )

    def BuildCoordinates ( self, coordinates3 ): pass

    def Distance ( self, other ):
        """Calculate the distance between two states."""
        dX = self.x - other.x
        dY = self.y - other.y
        dZ = self.z - other.z
        return math.sqrt ( dX*dX + dY*dY + dZ*dZ )

class AtomCenteredState ( PathwaysNode ):
    """An atom-centered state."""

    _attributable = dict ( PathwaysNode._attributable )
    _attributable.update ( { "atom"      : None ,
                             "atomIndex" :   -1 } )

    def BuildCoordinates ( self, coordinates3 ):
        """Build coordinates of the state."""
        i = self.atomIndex
        self.x = coordinates3[i,0]
        self.y = coordinates3[i,1]
        self.z = coordinates3[i,2]

class BondCenteredState ( PathwaysNode ):
    """A bond-centered state."""

    _attributable = dict ( PathwaysNode._attributable )
    _attributable.update ( { "atom1"      : None ,
                             "atom2"      : None ,
                             "atomIndex1" :   -1 ,
                             "atomIndex2" :   -1 ,
                             "bond"       : None ,
                             "bondIndex"  :   -1 ,
                             "bondLength" :  0.0 ,
                             "nX"         :  0.0 ,
                             "nY"         :  0.0 ,
                             "nZ"         :  0.0 } )

    def BuildCoordinates ( self, coordinates3 ):
        """Build coordinates of the state."""
        i  = self.atomIndex1
        j  = self.atomIndex2
        x1 = coordinates3[i,0]
        y1 = coordinates3[i,1]
        z1 = coordinates3[i,2]
        x2 = coordinates3[j,0]
        y2 = coordinates3[j,1]
        z2 = coordinates3[j,2]
        dX = x1 - x2
        dY = y1 - y2
        dZ = z1 - z2
        r  = math.sqrt ( dX*dX + dY*dY + dZ*dZ )
        self.bondLength = r
        self.nX         = dX / r
        self.nY         = dY / r
        self.nZ         = dZ / r
        self.x          = x2 + 0.5 * dX
        self.y          = y2 + 0.5 * dY
        self.z          = z2 + 0.5 * dZ

#===================================================================================================================================
# . Paths.
#===================================================================================================================================
class PathwaysPath ( Path ):
    "An electron transfer path."""

    pass

#===================================================================================================================================
# . State definitions.
#===================================================================================================================================
class PathwaysStateDefinition ( AttributableObject ):
    """A state definition."""
 
    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "energy" :  0.0 ,
                             "label"  : None } )

class AtomCenteredStateDefinition ( PathwaysStateDefinition ):
    """An atom-centered state definition."""
    pass

class BondCenteredStateDefinition ( PathwaysStateDefinition ):
    """A bond-centered state definition."""
    pass

#===================================================================================================================================
# . Graph.
#===================================================================================================================================
class PathwaysGraph ( Graph ):
    """A pathways graph."""
 
    _attributable = dict ( Graph._attributable )
    _attributable.update ( { "atomStateMapping" : None ,
                             "bondStateMapping" : None ,
                             "system"           : None } )

    def _CheckOptions ( self ):
        super ( PathwaysGraph, self )._CheckOptions ( )
        if len ( system.connectivity.bonds ) <= 0: raise PathwaysError ( "System must have a full connectivity." )

# . Through space too?

    def ApplySystemCoordinates ( self ):
        """Build the coordinates of the states and the link couplings."""
        coordinates3 = self.system.coordinates3
        if coordinates3 is None: raise PathwaysError ( "System must have coordinates3." )
        else:
            for state in self.nodes:
                state.BuildCoordinates ( coordinates3 )
            for link in self.edges:
                link.CalculateCoupling ( )

    def BuildStates ( self, atomTypes, atomCenteredStateDefinitions, bondCenteredStateDefinitions ):
        """Build states."""
        # . Initialization.
        atomNodes                   = set ( )
        missingAtomStateDefinitions = set ( )
        missingBondStateDefinitions = set ( )
        self.Clear ( )
        # . Atom-centered.
        if atomCenteredStateDefinitions is not None:
            for ( i, ( atom, atomType ) ) in enumerate ( zip ( self.system.atoms, atomTypes ) ):
                definition = atomCenteredStateDefinitions. get ( atomType, None )
                if definition is not None:
                    state = AtomCenteredState ( atom       = atom       ,
                                                atomIndex  = i          ,
                                                definition = definition )
                    self.AddNode ( state )
                    self.atomStateMapping[atom] = state
                    atomNodes.add ( i )
                else:
                    missingAtomStateDefinitions.add ( atomType )
        # . Bond-centered.
        if bondCenteredStateDefinitions is not None:
            for ( bondIndex, bond ) in enumerate ( self.system.bonds ):
                i = bond.i
                j = bond.j
                # . Skip only if both atoms are already assigned to atom states.
                if ( i not in atomNodes ) or ( j not in atomNodes ):
                    iType      = atomTypes[i]
                    jType      = atomTypes[j]
                    key        = ( max ( iType, jType ), min ( iType, jType ) )
                    definition = bondCenteredStateDefinitions. get ( key, None )
                    if definition is not None:
                        state = BondCenteredState ( atom1      = self.system.atoms[i],
                                                    atom2      = self.system.atoms[j],
                                                    atomIndex1 = i                   ,
                                                    atomIndex2 = j                   ,
                                                    bond       = bond                ,
                                                    bondIndex  = bondIndex           ,
                                                    definition = definition          )
                        self.AddNode ( state )
                        self.bondStateMapping[bond] = state
                    else:
                        missingBondStateDefinitions.add ( key )
            # . Missing bond definitions.
            if len ( missingBondStateDefinitions ) > 0:
                raise PathwaysError ( "There are missing bond state definitions.", missingBondStateDefinitions )
        # . Are there missing atom state definitions?
        elif len ( missingAtomStateDefinitions ) > 0:
            raise PathwaysError ( "There are missing atom state definitions.", missingAtomStateDefinitions )

    def BuildThroughBondLinks ( self, couplingFunctions ):
        """Build through-bond links."""
        if len ( self.nodes ) > 0:
            self.ClearEdges ( )
            missingCouplingFunctions = set ( )
            # . Loop over states.
            for state in self.nodes:
                if isinstance ( state, AtomCenteredState ):
                    bondIndices = self.system.bonds.GetConnections ( state.atomIndex )
                    for bondIndex in bondIndices:
                        bond       = self.system.bonds[bondIndex]
                        otherState = self.bondStateMapping.get ( bond, None )
                        if otherState is None:
                            otherAtomIndex = bond.Other ( state.atomIndex )
                            otherAtom      = self.system.atoms[otherAtomIndex]
                            otherState     = self.atomStateMapping[otherAtom]
                            if state.atomIndex <= otherState.atomIndex: otherState = None
                        if otherState is not None:
                            link = ThroughBondLink.FromCouplingFunctions ( state, otherState, couplingFunctions, missingCouplingFunctions )
                            if link is not None: self.AddEdge ( link )                          
                elif isinstance ( state, BondCenteredState ):
                    bond             = state.bond
                    otherBondIndices = set ( )
                    for ( atomIndex, atom ) in ( ( bond.i, self.system.atoms[bond.i] ) ,
                                                 ( bond.j, self.system.atoms[bond.j] ) ):
                        if atom not in self.atomStateMapping:
                            otherBondIndices.update ( self.system.bonds.GetConnections ( atomIndex ) )
                    otherBondIndices.remove ( state.bondIndex )
                    for otherBondIndex in otherBondIndices:
                        otherBond  = self.system.bonds[otherBondIndex]
                        otherState = self.bondStateMapping[otherBond]
                        if state.bondIndex > otherState.bondIndex:
                            link = ThroughBondLink.FromCouplingFunctions ( state, otherState, couplingFunctions, missingCouplingFunctions )
                            if link is not None: self.AddEdge ( link )
            # . Check for missing parameters.
            if len ( missingCouplingFunctions ) > 0: raise PathwaysError ( "There are missing through-bond link definitions.", missingCouplingFunctions )

#    def BuildThroughSpaceLinks ( self ):
#                """Build through-space links."""
#                Get coordinates of states
#                Get pairlist of excluded interactions (throughbond)
#                Build pairlist of through-space interactions
#                Build links given pairlist

#        For this need index of each state (in addition to key)

#* Define TS link types as function of two state types. Requires coupling function.
#  See method for pair generation.

#   def GetEdgeByKey ( self, key ): pass
#   def GetNodeByKey ( self, key ): pass


# . Need these - or do directly?
# . Do best path first to find weight, then set limits.
#   def FindAllPaths ( self, acceptor, donor, cutOff = _MinimumWeight, maximumWeight = _MaximumWeight, minimumWeight = _MinimumWeight ):
#       """Find all paths between donor and acceptor with weights in a specific range."""
#       allPaths = DijkstraSingleSource ( self, donor, cutOff        = cutOff       ,
#                                                      maximumWeight = maximumWeight,
#                                                      minimumWeight = minimumWeight,
#                                                      target        = acceptor )
#       return allPaths

#   def FindBestPath ( self, acceptor, donor, cutOff = _MinimumWeight ):
#       """Search for the best path between donor and acceptor states."""
#       bestPath = DijkstraSingleSource ( self, donor, cutOff = cutOff, target = acceptor )
#       return bestPath

    def Clear ( self ):
        """Clear all edge and node data."""
        super ( PathwaysGraph, self ).Clear ( )
        self.atomStateMapping = {}
        self.bondStateMapping = {}

    @classmethod
    def FromSystem ( selfClass, system, **options ):
        """Constructor from system."""
        options = dict ( options )
        options["system"] = system
        return selfClass.WithOptions ( **options )

    def Summary ( self, log = logFile ):
        """Summary."""
        if LogFileActive ( log ):
            log.SummaryOfItems ( [ ( "States", "{:d}".format ( len ( self.nodes ) ) )   ,
                                   ( "Links" , "{:d}".format ( len ( self.edges ) ) ) ] ,
                                     title = "Pathways Graph Summary" )

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass
