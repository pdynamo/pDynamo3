"""Utilities for handling a connectivity."""

# . Some are put here temporarily while their proper home is decided ...

from  collections       import defaultdict
from  enum              import Enum
from  itertools         import product
from  pCore             import logFile                           , \
                               LogFileActive
from  pScientific       import IsMainGroup                       , \
                               PeriodicTable
from  pScientific.Graph import EdmondsMaximumMatching            , \
                               VismaraRelevantCycles
from .AromaticPattern   import AromaticPatterns                  , \
                               AromaticPatternsImplicitHydrogens
from .Bond              import Atom                              , \
                               Bond                              , \
                               BondType

#===================================================================================================================================
# . Definitions.
#===================================================================================================================================
# . Atom geometry types.
#
# . Many others possible:
#
#    3 trigonal pyramid, T-shaped                                
#    6 hexagonal planar, trigonal prism                          
#    7 capped octahedron, capped trigonal prism                  
#    8 dodecahedron, cube, square antiprism, hexagonal bipyramid 
#   10 bicapped square antiprism                                 
#   11 all-faced capped trigonal prism                           
#   12 cuboctahedron                                             
#
class AtomGeometryType ( Enum ):
    """Atom geometry types."""
    Cubic                  = ( "Cubic"                  , 8 , "Cub" )
    CubicAntiprism         = ( "CubicAntiprism"         , 8 , "CAp" )
    Isolated               = ( "Isolated"               , 0 , "Iso" )
    Linear                 = ( "Linear"                 , 2 , "Lin" )
    Octahedral             = ( "Octahedral"             , 6 , "Oct" )
    PentagonalBipyramidal  = ( "PentagonalBipyramidal"  , 7 , "PBp" )
    Resonant               = ( "Resonant"               , 3 , "Res" )
    SquarePlanar           = ( "SquarePlanar"           , 4 , "SPl" )
    SquarePyramidal        = ( "SquarePyramidal"        , 5 , "SPy" )
    Terminal               = ( "Terminal"               , 1 , "Ter" )
    Tetrahedral            = ( "Tetrahedral"            , 4 , "Tet" ) 
    TricappedTrigonalPrism = ( "TricappedTrigonalPrism" , 9 , "TTP" )
    Trigonal               = ( "Trigonal"               , 3 , "Tri" )
    TrigonalBipyramidal    = ( "TrigonalBipyramidal"    , 5 , "TBp" )
    Undefined              = ( "Undefined"              , 0 , ""    )

    def __init__ ( self, label, neighbors, tag ):
        """Constructor."""
        self.label     = label
        self.neighbors = neighbors
        self.tag       = tag

# . From tag to type.
AtomGeometryTypeFromTag = { agt.tag : agt for agt in AtomGeometryType }

#===================================================================================================================================
# . Parameters and utility functions.
#===================================================================================================================================
# . Default atom geometries given a coordination number.
_DefaultAtomGeometries = { 0 : "Iso" ,
                           1 : "Ter" ,
                           2 : "Lin" ,
                           3 : "Tri" ,
                           4 : "Tet" ,
                           5 : "TBp" ,
                           6 : "Oct" ,
                           7 : "PBp" ,
                           8 : "Cub" ,
                           9 : "TTP" }

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class ConnectivityError ( Exception ):
    pass

#===================================================================================================================================
# . Convert implicit into explicit hydrogens.
#===================================================================================================================================
def AddExplicitHydrogens ( self ):
    """Add explicit hydrogens."""
    self.ClearDerivableAttributes ( ) # . As the connectivity is being changed.
    nNodes = len ( self.nodes )
    for atom in self.nodes[0:nNodes]: # . As nodes is extended dynamically.
        implicitHs = atom.__dict__.pop ( "implicitHydrogens", 0 )
        for h in range ( implicitHs ):
            hAtom = Atom.WithOptions ( atomicNumber = 1 )
            self.AddNode ( hAtom )
            self.AddEdge ( Bond.WithNodes ( atom, hAtom, type = BondType.Single ) )

#===================================================================================================================================
# . Check for Kekule output - used for MOL, MOL2 and SMILES (which do not have specificiations for aromatic triple bonds).
#===================================================================================================================================
def CheckForKekuleOutput ( self, useKekule ):
    """Determine the Kekule output option."""
    forceKekule = useKekule
    if not forceKekule:
        for bond in self.edges:
            if bond.isAromatic and ( bond.type.bondOrder > 2 ):
                forceKekule = True
                break
    return forceKekule

#===================================================================================================================================
# . Clear various atom attributes.
#===================================================================================================================================
def ClearAtomAttributes ( self, *attributes ):
    """Clear atom attributes."""
    for atom in self.nodes:
        for attribute in attributes: atom.__dict__.pop ( attribute, None )

#===================================================================================================================================
# . Convert an input connectivity.
#===================================================================================================================================
#
# . Many formats define aromaticity and use implicit hydrogens (e.g. MOL, MOL2, Smiles).
# . This method tries to convert one of these input connectivities to pDynamo form.
#
# . Problems that arise include:
#
#   - The specification of atoms or bonds as being aromatic without the equivalent
#     bond or atom definition.
#
#   - Aromatic bonds have uncertain bond order (either 1 or 2 normally in the resulting
#     Kekule structure).
#
#   - The number of implicit hydrogens is not specified but needs to be deduced.
#
# . Aromaticity is checked only for atoms and bonds which have been specified as being aromatic.
# . pDynamo's standard aromaticity method is then applied afterwards on the resulting Kekule structure.
#
def ConvertInputConnectivity ( self, organicValencies, log = logFile ):
    """Convert an input connectivity to pDynamo form."""
    if len ( self.nodes ) > 0:
        # . Ensure atom and bond aromaticity flags are consistent and check for undefined bonds.
        MakeAtomAttributes ( self, "index" )
        for bond in self.edges:
            node1 = bond.node1
            node2 = bond.node2
            if bond.isAromatic:
                node1.isAromatic = True
                node2.isAromatic = True
            elif node1.isAromatic and node2.isAromatic: # . This could be wrong but will be corrected later.
                bond.isAromatic = True
        # . Gather some atom information.
        MakeAtomAttributes ( self, "inRing" )
        certainAromatics      = set ( )
        rejectedAromatics     = set ( )
        undefinableImplicitHs = set ( )
        uncertainAromatics    = set ( )
        for atom in self.nodes:
            bonds                     = self.adjacentEdges[atom]
            others                    = self.adjacentNodes[atom]
            bondOrderSum              = sum ( [ bond.type.bondOrder for bond in bonds ] ) # . This includes undefined bond orders.
            implicitHs                = getattr ( atom, "implicitHydrogens", 0 )
            neighbors                 = len ( others )
            undefinedAromaticBonds    = len ( [ 1 for bond in bonds if (       bond.isAromatic   and ( bond.type is BondType.Undefined ) ) ] )
            undefinedNonAromaticBonds = len ( [ 1 for bond in bonds if ( ( not bond.isAromatic ) and ( bond.type is BondType.Undefined ) ) ] )
            unknownImplicitHs         = getattr ( atom, "isReduced", False )
            # . Aromatic atoms.
            if atom.isAromatic:
                connections = implicitHs + neighbors
                # . Accepted.
                if ( connections in ( 2, 3 ) ) and atom.inRing: # . Basic pre-checks.
                    # . Flag aromatics that might have unknown hydrogen count.
                    # . The assumption is made here that if the all bonds to aromatic atoms
                    # . have been explicitly specified then there is no uncertainty. 
                    if unknownImplicitHs and ( ( connections != 2 ) or ( undefinedAromaticBonds == 0 ) ):
                        atom.__dict__.pop ( "isReduced", None ) # . No longer unknown (= 0).
                        unknownImplicitHs = False
                    if unknownImplicitHs: uncertainAromatics.add ( atom )
                    else:                   certainAromatics.add ( atom )
                # . Rejected.
                else:
                    rejectedAromatics.add ( atom )
                    atom.isAromatic = False
                    for bond in bonds:
                        bond.isAromatic = False
                    undefinedNonAromaticBonds += undefinedAromaticBonds
                    undefinedAromaticBonds     = 0
            # . Determine connections and implicit hydrogens for as many atoms as possible.
            atom.connections = neighbors
            if unknownImplicitHs:
                if not atom.isAromatic:
                    valence = bondOrderSum + undefinedNonAromaticBonds
                    for v in organicValencies[atom.atomicNumber]:
                        hCount = v - valence
                        if hCount >= 0:
                            atom.implicitHydrogens = hCount
                            atom.connections      += hCount
                            break
                    atom.__dict__.pop ( "isReduced", None ) # . No longer uncertain.
                    # . If there are undefined bonds we cannot guess the hydrogen number precisely but we still store the connectivity value.
                    if undefinedNonAromaticBonds > 0: undefinableImplicitHs.add ( atom )
            else:
                atom.connections += implicitHs
        # . Remove all aromaticity flags.
        for atom in self.nodes: atom.isAromatic = False
        for bond in self.edges: bond.isAromatic = False
        # . Verify aromaticity for specified atoms only!
        # . This will automatically exclude partial specifications within the same ring.
        electronsC = {}
        electronsU = {}
        if len ( uncertainAromatics ) > 0:
            electronsU = AromaticPatternsImplicitHydrogens.TypeConnectivity ( self, uncertainAromatics )
            certainAromatics.update ( uncertainAromatics ) # . Update certain aromatics with those atoms that have not been matched.
        if len ( certainAromatics ) > 0:
            electronsC = AromaticPatterns.TypeConnectivity ( self, certainAromatics )
        rejectedAromatics.update ( certainAromatics ) # . Further rejected atoms.
        # . Get the subgraph of possible aromatics and loop over rings.
        kekuleM   = 0
        kekuleP   = 0
        possibleC = set ( electronsC.keys ( ) )
        possibleU = set ( electronsU.keys ( ) )
        possible  = possibleC.union ( possibleU )
        uMultiple = False
        if len ( possible ) > 0:
            aGraph = self.MakeSubgraph ( possible )
            for ringSet in aGraph.ringSets:
                # . Loop over rings and preprocess for certain aromatics.
                piBonds = {}
                rings   = VismaraRelevantCycles ( aGraph, biconnectedComponents = [ ringSet ] )
                # . Add all atoms in the ring set if there are more than one ring.
                if len ( rings ) > 1: rings = [ ringSet ] + rings
                ringsU  = []
                for ring in rings:
                    nElectrons = 0
                    ringAtoms  = set ( ring )
                    ringAtomsC = ringAtoms.intersection ( possibleC )
                    ringAtomsU = ringAtoms.intersection ( possibleU )
                    for atom in ringAtomsC: nElectrons += electronsC[atom][0]
                    # . There are no uncertain aromatics so process immediately.
                    if len ( ringAtomsU ) == 0:
                        if ( nElectrons >= 2 ) and ( ( nElectrons - 2 ) % 4 == 0 ):
                            for atom in ring:
                                atom.isAromatic = True
                                piBonds[atom]   = electronsC[atom][1]
                                for bond in self.adjacentEdges[atom]:
                                    if bond.Opposite ( atom ) in ring:
                                        bond.isAromatic = True
                    # . Save the uncertain aromatics for the ring.
                    else: ringsU.append ( ( ringAtomsU, ring, nElectrons ) )
                # . Find all possible combinations of uncertain aromatics.
                # . Hopefully short except in rare cases.
                uAtoms = sorted ( list ( set ( ringSet ).intersection ( possibleU ) ) )
                if len ( uAtoms ) > 0:
                    uCases  = list ( product ( *[ list ( range ( len ( electronsU[uAtom] ) ) ) for uAtom in uAtoms ] ) )
                    uTotals = defaultdict ( list )
                    # . Loop over the different uncertain cases.
                    for ( u, uCase ) in enumerate ( uCases ):
                        uIndices = { uAtom : index for ( uAtom, index ) in zip ( uAtoms, uCase ) }
                        flagged  = set ( )
                        for ( ringAtomsU, ring, nElectrons ) in ringsU:
                            for uAtom in ringAtomsU: nElectrons += electronsU[uAtom][uIndices[uAtom]][0]
                            if ( nElectrons >= 2 ) and ( ( nElectrons - 2 ) % 4 == 0 ):
                                for atom in ring: flagged.add ( atom )
                        uTotals[len ( flagged )].append ( u )
                    # . Find the case that maximizes the number of aromatic atoms.
                    uMaximum = max ( uTotals.keys ( ) )
                    if uMaximum > 0:
                        uBest = uTotals[uMaximum]
                        if len ( uBest ) > 1: uMultiple = True
                        # . Flag aromatic atoms and bonds for the case.
                        # . Choose the first best case if there are multiple cases.
                        uCase    = uCases[uBest[0]]
                        uIndices = { uAtom : index for ( uAtom, index ) in zip ( uAtoms, uCase ) }
                        for ( ringAtomsU, ring, nElectrons ) in ringsU:
                            for uAtom in ringAtomsU: nElectrons += electronsU[uAtom][uIndices[uAtom]][0]
                            if ( nElectrons >= 2 ) and ( ( nElectrons - 2 ) % 4 == 0 ):
                                for atom in ring:
                                    if not atom.isAromatic:
                                        atom.isAromatic = True
                                        if atom in uIndices:
                                            uData         = electronsU[atom][uIndices[atom]]
                                            hCount        = uData[2]
                                            piBonds[atom] = uData[1]
                                            atom.implicitHydrogens = hCount
                                            atom.connections      += hCount
                                            atom.__dict__.pop ( "isReduced", None ) # . No longer uncertain.
                                        else:
                                            piBonds[atom] = electronsC[atom][1]
                                    for bond in self.adjacentEdges[atom]:
                                        if bond.Opposite ( atom ) in ring:
                                            bond.isAromatic = True
                # . Set up a Kekule representation of the ring set if one does not already exist.
                isOK        = True
                rsAromatics = piBonds.keys ( )
                for atom in rsAromatics:
                    nPiBonds   = 0
                    nUndefined = 0
                    bonds      = [ bond for bond in self.adjacentEdges[atom] if bond.Opposite ( atom ) in rsAromatics ]
                    nPiBonds   = len ( [ 1 for bond in bonds if bond.type.bondOrder > 1         ] )
                    nUndefined = len ( [ 1 for bond in bonds if bond.type is BondType.Undefined ] )
                    if ( nPiBonds != piBonds[atom] ) or ( nUndefined > 0 ):
                        isOK = False
                        for bond in bonds:
                            bond.isAromatic = True
                            bond.type       = BondType.Single
                if not isOK:
                    rsNodes = [ atom for atom in rsAromatics if piBonds[atom] != 0 ]
                    rsGraph = self.MakeSubgraph ( rsNodes )
                    result  = EdmondsMaximumMatching ( rsGraph )
                    if len ( result ) != len ( rsNodes ): kekuleM += 1
                    else:                                 kekuleP += 1
                    for ( iAtom, jAtom ) in result.items ( ):
                        for bond in self.adjacentEdges[iAtom]:
                            if bond.Opposite ( iAtom ) is jAtom:
                                bond.isAromatic = True
                                bond.type       = BondType.Double
        # . Save all possible aromatics that were not flagged as aromatic.
        for atom in possible:
            if not atom.isAromatic: rejectedAromatics.add ( atom )
        # . Treat remaining uncertain atoms (all non-aromatic).
        for atom in self.nodes:
            if getattr ( atom, "isReduced", False ):
                atom.connections = len ( self.adjacentNodes[atom] )
                bondOrderSum     = sum ( [ bond.type.bondOrder for bond in self.adjacentEdges[atom] ] )
                undefinedBonds   = len ( [ 1 for bond in  self.adjacentEdges[atom] if bond.type is BondType.Undefined ] )
                valence          = bondOrderSum + undefinedBonds
                for v in organicValencies[atom.atomicNumber]:
                    hCount = v - valence
                    if hCount >= 0:
                        atom.implicitHydrogens = hCount
                        atom.connections      += hCount
                        break
                atom.__dict__.pop ( "isReduced", None ) # . No longer uncertain.
                if undefinedBonds > 0: undefinableImplicitHs.add ( atom )
        # . Are there any more undefined bonds?
        numberOfUndefinedBonds = 0
        for bond in self.edges:
            if bond.type is BondType.Undefined:
                if ( not bond.isAromatic ) and bond.node1.isAromatic and bond.node2.isAromatic:
                    bond.type = BondType.Single # . Bond between aromatic ring sets. This should be sufficient in almost all cases.
                else:
                    numberOfUndefinedBonds += 1 # . Cannot do much without more data.
        # . Output.
        hasError   = ( numberOfUndefinedBonds > 0 ) or ( len ( undefinableImplicitHs ) > 0 )
        hasWarning = uMultiple or ( kekuleM > 0 ) or ( len ( rejectedAromatics ) > 0 )
        if LogFileActive ( log ) and ( hasError or hasWarning ):
            items = [ ( "Multiple Aromatic Assignment" , "{:s}".format ( repr ( uMultiple            ) ) ) ,
                      ( "Rejected Aromatics"           , "{:d}".format ( len ( rejectedAromatics     ) ) ) ,
                      ( "Undefined Bonds"              , "{:d}".format ( numberOfUndefinedBonds        ) ) ,
                      ( "Undefined Implicit Hydrogens" , "{:d}".format ( len ( undefinableImplicitHs ) ) ) ,
                      ( "Kekule Maximum Matches"       , "{:d}".format ( kekuleM                       ) ) ,
                      ( "Kekule Perfect Matches"       , "{:d}".format ( kekuleP                       ) ) ]
            log.SummaryOfItems ( items, title = "Error/Warning Summary of Input Connectivity Conversion" )
        # . Check for error.
        if hasError: raise ConnectivityError ( "Error converting input connectivity (undefined bonds or implicit hydrogens)." )
        # . Finish up.
        self.CompleteConnectivity ( )

#===================================================================================================================================
# . Aromaticity.
#===================================================================================================================================
def DetermineAromaticity ( self ):
    """Determine which atoms of a connectivity are aromatic."""
    # . Initialization.
    atoms         = self.atoms
    bonds         = self.bonds
    neighborBonds = self.adjacentEdges
    ringSets      = self.ringSets
    for atom in atoms:
        atom.inRing     = False
        atom.isAromatic = False
    for bond in bonds:
        bond.isAromatic = False
    # . Check that there are rings.
    if ( ringSets is not None ) and ( len ( ringSets ) > 0 ):
        # . Loop over each ring set separately.
        for ringSet in ringSets:
            # . Initialization - ringSet duplicated because ringAtoms is consumed by TypeConnectivity.
            ringAtoms = set ( ringSet )
            for atom in ringAtoms: atom.inRing = True
            # . Determine which atoms could be aromatic.
            electrons = AromaticPatterns.TypeConnectivity ( self, ringAtoms )
            possible  = set ( electrons.keys ( ) )
            # . Proceed if there are possible aromatics.
            nPossible = len ( possible )
            if nPossible <= 0: continue
            # . Get the rings.
            rGraph = self.MakeSubgraph ( possible )
            rings  = VismaraRelevantCycles ( rGraph )
            # . If all atoms are possible be optimistic and try all atoms first ...
            if ( len ( ringAtoms ) == 0 ) and ( len ( rings ) > 1 ):
                rings = [ possible ] + rings
            # . For each ring check if the number of electrons obeys Huckel's rule.
            nAromatic = 0
            for ring in rings:
                nElectrons = 0
                for atom in ring:
                    nElectrons += electrons[atom][0]
                # . Flag ring atoms and bonds as aromatic if Huckel's rule is obeyed.
                # . Bonds have to be done here as bonds between aromatic atoms
                # . are not necessarily aromatic (e.g. biphenyl).
                if ( nElectrons >= 2 ) and ( ( nElectrons - 2 ) % 4 == 0 ):
                    for atom in ring:
                        if not atom.isAromatic:
                            atom.isAromatic = True
                            nAromatic += 1
                        for bond in neighborBonds[atom]:
                            if bond.Opposite ( atom ) in ring:
                                bond.isAromatic = True
                    if nAromatic >= nPossible: break

#===================================================================================================================================
# . Atom geometry.
#===================================================================================================================================
# . Needs improvement.
def DetermineAtomGeometry ( self ):
    """Determine the bonding geometry of the atoms in a connectivity."""
    # . Loop over atoms.
    for atom in self.atoms:
        # . Attributes.
        atomicNumber = atom.atomicNumber
        connections  = getattr ( atom, "connections", None  )
        isAromatic   = getattr ( atom, "isAromatic",  False )
        valence      = getattr ( atom, "valence",     None  )
        # . Tests.
        if ( connections is not None ) and ( valence is not None ):
            # . Basic cases.
            # . Isolated atom.
            if connections == 0:
                atom.geometry = AtomGeometryType.Isolated.tag
            # . Resonant.
            elif isAromatic:
                atom.geometry = AtomGeometryType.Resonant.tag
            # . Main-group elements - use VSEPR.
            elif IsMainGroup ( atomicNumber ):
                # . Get the number of available electrons (number of valence electrons - formal charge - one per bond valence.
                electrons = max ( 0, getattr ( PeriodicTable[atomicNumber], "valenceElectrons" ) - atom.formalCharge - valence )
                # . Find the number of lone pairs and single electrons.
                ( lonePairs, unpaired ) = divmod ( electrons, 2 )
                # . Find the number of sites around the atom and get the corresponding geometry.
                sites = connections + lonePairs + unpaired
                try:    atom.geometry = _DefaultAtomGeometries[sites]
                except: pass
            # . Other cases by connections only.
            else:
                try:    atom.geometry = _DefaultAtomGeometries[connections]
                except: pass

#===================================================================================================================================
# . Atom oxidation states.
#===================================================================================================================================
# . Needs improvement.
# . For the definitive definition of OS see: Karen P, McArdle P, Takats J, Pure Appl Chem 86, 1017-1081, 2014.
def DetermineAtomOxidationState ( self ):
    """Determine the oxidation states of the atoms in a connectivity."""
    # . Initialization.
    atoms = self.atoms
    bonds = self.bonds
    # . Initialization the oxidation states of each atom.
    for atom in atoms: atom.oxidationState = atom.formalCharge
    # . Loop over bonds giving the more electronegative atom the electrons.
    for bond in bonds:
        atom1 = bond.node1
        atom2 = bond.node2
        if atom1.atomicNumber != atom2.atomicNumber:
            en1   = getattr ( atom1, "Allen Electronegativity", 2.0 )
            en2   = getattr ( atom2, "Allen Electronegativity", 2.0 )
            order = bond.type.bondOrder
            if en1 > en2:
                atom1.oxidationState -= order
                atom2.oxidationState += order
            else:
                atom1.oxidationState += order
                atom2.oxidationState -= order

#===================================================================================================================================
# . Make various atom attributes.
#===================================================================================================================================
def MakeAtomAttributes ( self, *attributes ):
    """Make miscellaneous atom attributes."""
    for attribute in attributes:
        if attribute == "connections":
            for atom in self.nodes:
                implicitHs       = getattr ( atom, "implicitHydrogens", 0 )
                atom.connections = len ( self.adjacentEdges[atom] ) + implicitHs
        elif attribute == "hydrogens":
            for atom in self.nodes:
                implicitHs     = getattr ( atom, "implicitHydrogens", 0 )
                explicitHs     = len ( [ 1 for neighbor in self.adjacentNodes[atom] if neighbor.atomicNumber == 1 ] )
                atom.hydrogens = explicitHs + implicitHs
        elif attribute == "index":
            for ( i, atom ) in enumerate ( self.nodes ): atom.index = i
        elif attribute == "inRing":
            for atom in self.nodes: atom.inRing = False
            for ringSet in self.ringSets:
                for atom in ringSet: atom.inRing = True
        elif attribute == "valence":
            for atom in self.nodes:
                implicitHs   = getattr ( atom, "implicitHydrogens", 0 )
                atom.valence = sum ( [ bond.type.bondOrder for bond in self.adjacentEdges[atom] ] ) + implicitHs
        else: raise ConnectivityError ( "Unknown atom attribute - {:s}.".format ( attribute ) )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
