"""Script to find CIP labels for a system.

A concise, clear nice summary may be found in:

"Preferred IUPAC Names", Chapter 9, IUPAC (http://www.iupac.org)

The last CIP paper published by one of the eponymous authors:

V. Prelog and G. Helmchen
"Basic Principles of the CIP-System and Proposals for a Revision."
Angew. Chem. Int. Ed. Engl. 1982, 21, 567-583.

A paper indicating some of the complications of applying them:

P. Mata, A. M. Lobo, C. Marshall and A. P. Johnson
"Implementation of the Cahn-Ingold-Prelog System for Stereochemical
Perception in the LHASA Program."
J. Chem. Inf. Comp. Sci. 1994, 34, 491-504.

The rules for assignment of priority are complicated and only those
that make use of connectivity information are implemented here. This
means that priorities are determined using atomic numbers, atomic
masses and bond connectivities converted to hierarchical digraphs.

Omitted are complications involving multiple Kekule structures for
a system (e.g. fractional dummy atoms) and all priority comparisons
involving stereocenters. The latter include things such as the fact
that like descriptor pairs (R/R, S/S) take precedence over (R/S, S/R).

It would be good to try and differentiate between centers that are
clearly achiral (e.g. with > 1 bound Hs) and those that cannot be
determined.

Tetrahedral centers - lowest priority substituent at back with
                      a right-handed coordinate system:

   R (rectus)   - clockwise (highest to lowest).
   S (sinister) - anticlockwise.

Double bond centers - two highest priority substituents:

   E (Entgegen) - trans.
   Z (Zusammen) - cis.

Amino acids in proteins are S at their Calpha atoms, except for
cysteine which is R. Isoleucine/threonine - sidechains?

"""

import copy, math, operator

from pCore     import LogFileActive, logFile
from pMolecule import BondType

#===============================================================================
# . Public functions.
#===============================================================================
def CIPLabelFinder ( system, log = logFile ):
    """Determine CIP labels for a system."""

    # . The system must have full connectivity.
    if len ( system.connectivity.bonds ) > 0:

        # . The system should be kekularized.
        #system.Kekularize ( )

        # . Treat the tetrahedral centers.
        ( tcenters, rtcenters, stcenters, utcenters ) = _TetrahedralCenters ( system )

        # . Treat the double bond centers.
        ( dcenters, edcenters, zdcenters, udcenters ) = _DoubleBondCenters  ( system )

        # . Output results.
        if LogFileActive ( log ):
            items = [ ( "Tetrahedral Centers"  , "{:d}".format ( len ( tcenters  ) ) ) ]
            if len ( tcenters ) > 0:
               items.extend ( [ ( "Undefined Tetrahedral Centers" , "{:d}".format ( len ( utcenters ) ) ) ,
                                ( "R Centers"                     , "{:d}".format ( len ( rtcenters ) ) ) ,
                                ( "S Centers"                     , "{:d}".format ( len ( stcenters ) ) ) ] )
            items.append ( ( "Double Bond Centers"  , "{:d}".format ( len ( dcenters  ) ) ) )
            if len ( dcenters ) > 0:
               items.extend ( [ ( "Undefined Double Bond Centers" , "{:d}".format ( len ( udcenters ) ) ) ,
                                ( "E Centers"                     , "{:d}".format ( len ( edcenters ) ) ) ,
                                ( "Z Centers"                     , "{:d}".format ( len ( zdcenters ) ) ) ] )
            log.SummaryOfItems ( items, order = False, title = "Cahn-Ingold-Prelog Labels" )

        # . Finish up.
        return ( ( tcenters, rtcenters, stcenters, utcenters ), ( dcenters, edcenters, zdcenters, udcenters ) )

    # . Do nothing.
    else:
        return None

#===============================================================================
# . Private functions.
#===============================================================================
def _AssignCIPPriorities ( branches, highPriority, lowPriority ):
    """Order a set of atoms with respect to their CIP priorities."""
    # . Sort the branches in order of increasing priority.
    branches.sort ( key = operator.attrgetter ( "levelData" ) )
    # . Check for branches of lower and higher priority.
    while True:
        if len ( branches ) > 1:
            if branches[0].levelData < branches[1].levelData:
                branch = branches.pop ( 0 )
                lowPriority.insert ( 0, branch.level1 )
            else: break
        else: break
    while True:
        if len ( branches ) > 1:
            if branches[-1].levelData > branches[-2].levelData:
                branch = branches.pop ( -1 )
                highPriority.append ( branch.level1 )
            else: break
        else: break
    # . Check for an odd branch.
    if len ( branches ) == 1:
        branch = branches.pop ( )
        highPriority.append ( branch.level1 )
    # . There are no more branches to assign.
    if len ( branches ) == 0:
        return highPriority + lowPriority
    # . The next level is needed.
    else:
        # . Construct the next level for each branch.
        newbranches = []
        for branch in branches:
            branch.IncrementLevel ( )
            if len ( branch ) > 0: newbranches.append ( branch )
        if len ( newbranches ) > 0: return _AssignCIPPriorities ( newbranches, highPriority, lowPriority )
        else:                       return None

def _CIPPriorityOrderer ( system, rootAtom, branchAtoms ):
    """Order a set of atoms with respect to their CIP priorities."""
    # . Initialization.
    branches     = []
    highPriority = []
    lowPriority  = []
    # . Generate the branches.
    for level1 in branchAtoms:
        branches.append ( _CIPBranch ( system.connectivity, rootAtom, level1 ) )
     # . Assign priorities.
    return _AssignCIPPriorities ( branches, highPriority, lowPriority )

def _DoubleBondCenters ( system ):
    """Treat the double bond centers."""
    # . Find all double bonds between two sp2 carbons.
    dcenters  = []
    edcenters = []
    zdcenters = []
    udcenters = []
    for ( b, bond ) in enumerate ( system.connectivity.bonds ):
        if bond.type == BondType.Double:
            iAtom = bond.node1
            jAtom = bond.node2
            i     = system.connectivity.nodeIndices[iAtom]
            j     = system.connectivity.nodeIndices[jAtom]
            if ( iAtom.atomicNumber == 6 ) and ( iAtom.connections == 3 ) and \
               ( jAtom.atomicNumber == 6 ) and ( jAtom.connections == 3 ):
               # . Save the bond.
               dcenters.append ( b )
               # . Get the connected atoms.
               iAtoms = set ( system.connectivity.adjacentNodes[iAtom] )
               iAtoms.remove ( jAtom )
               jAtoms = set ( system.connectivity.adjacentNodes[jAtom] )
               jAtoms.remove ( iAtom )
               # . Assign priorities.
               oIAtoms = _CIPPriorityOrderer ( system, iAtom, iAtoms )
               oJAtoms = _CIPPriorityOrderer ( system, jAtom, jAtoms )
               if ( oIAtoms is not None ) and ( oJAtoms is not None ):
                   label = _DoubleBondCenterLabel ( system.coordinates3, i, system.connectivity.nodeIndices[oIAtoms[0]] ,
                                                                         j, system.connectivity.nodeIndices[oJAtoms[0]] )
                   if label == "E": edcenters.append ( b )
                   else:            zdcenters.append ( b )
               # . Cannot assign priorities to the atoms.
               else: udcenters.append ( b )
    return ( dcenters, edcenters, zdcenters, udcenters )

def _DoubleBondCenterLabel ( coordinates3, d1, s1, d2, s2 ):
    """Determine the double bond label - E or Z.

    d1 and d2 are the double bond atoms and s1 and s2 are
    the highest priority substitutents on each atom.
    """
    v1 = coordinates3.Displacement ( s1, d1 )
    v2 = coordinates3.Displacement ( s2, d2 )
    if v1.Dot ( v2 ) < 0.0: return "E"
    else:                   return "Z"

def _SignedVolume ( a, b, c ):
    """Get the signed volume of three vectors."""
    sv = a[0] * ( b[1] * c[2] - b[2] * c[1] ) + \
         a[1] * ( b[2] * c[0] - b[0] * c[2] ) + \
         a[2] * ( b[0] * c[1] - b[1] * c[0] )
    return sv

def _TetrahedralCenters ( system ):
    """Treat the tetrahedral centers."""
    # . Find all tetrahedral carbons.
    tcenters  = []
    rtcenters = []
    stcenters = []
    utcenters = []
    for ( i, atom ) in enumerate ( system.atoms ):
        if ( atom.atomicNumber == 6 ) and ( atom.connections == 4 ):
            # . Save the atom.
            tcenters.append ( i )
            # . Order the connected atoms by their priorities.
            iAtoms = system.connectivity.adjacentNodes[atom]
            oAtoms = _CIPPriorityOrderer ( system, atom, iAtoms )
            if oAtoms is not None:
                label = _TetrahedralCenterLabel ( system.coordinates3, i, system.connectivity.nodeIndices[oAtoms[0]] ,
                                                                          system.connectivity.nodeIndices[oAtoms[1]] ,
                                                                          system.connectivity.nodeIndices[oAtoms[2]] )
                if label == "R": rtcenters.append ( i )
                else:            stcenters.append ( i )
            # . Cannot assign priorities to the atoms.
            else: utcenters.append ( i )
    return ( tcenters, rtcenters, stcenters, utcenters )

def _TetrahedralCenterLabel ( coordinates3, t, s1, s2, s3 ):
    """Determine the tetrahedral bond label - R or S.

    t is the center and s1, s2 and s3 are the highest
    priority substituents (1 > 2 > 3).

    Assigning a signed volume of a particular sign to "R" or "S"
    relies on the definition of the coordinate axes. Here the
    standard x (towards), y (right), z (up) system is used.
    """
    a = coordinates3.Displacement ( s1, t )
    b = coordinates3.Displacement ( s2, t )
    c = coordinates3.Displacement ( s3, t )
    sv = _SignedVolume ( a, b, c )
    if sv < 0.0: return "R"
    else:        return "S"

#===============================================================================
# . Class.
#===============================================================================
# . Comparison is done on the levelData.

class _CIPBranch:
    """Class to hold a hierarchical digraph originating from an atom bound to a stereocenter."""

    def __init__ ( self, connectivity, level0, level1 ):
        """Constructor.

        |connectivity| is the connectivity.
        |level0|       is the stereocenter.
        |level1|       is the atom from which the branch originates.
        """
        self.connectivity = connectivity
        self.level0       = level0
        self.level1       = level1
        self.paths        = [ [ level0, level1 ] ]
        self.levelData    = [ ( level1.atomicNumber, level1.mass ) ]

    def __len__ ( self ):
        """Return the number of paths in the branch."""
        return len ( self.paths )

    def IncrementLevel ( self ):
        """Go to the next level."""
        # . Get the new paths.
        oldPaths   = self.paths
        self.paths = []
        for oldPath in oldPaths:
            # . Get the top atom and its source.
            oldSource = oldPath[-2]
            oldTop    = oldPath[-1]
            # . A dummy atom terminates a path.
            if isinstance ( oldTop, _CIPDummyAtom ):
                pass
            # . A regular atom.
            else:
                # . Loop over the bonds for the atom.
                for bond in self.connectivity.adjacentEdges[oldTop]:
                    otherAtom = bond.Opposite ( oldTop )
                    # . Multiple bonds add dummy atom paths no matter whether the bond is old or new.
                    for i in range ( bond.type.bondOrder - 1 ):
                        newpath = copy.copy ( oldPath )
                        newpath.append ( _CIPDummyAtom ( otherAtom.atomicNumber, otherAtom.mass ) )
                        self.paths.append ( newpath )
                    # . This is a new bond.
                    if otherAtom != oldSource:
                        # . A ring adds a dummy atom path.
                        if otherAtom in oldPath:
                            newpath = copy.copy ( oldPath )
                            newpath.append ( _CIPDummyAtom ( otherAtom.atomicNumber, otherAtom.mass ) )
                            self.paths.append ( newpath )
                        # . Add a normal path.
                        else:
                            newpath = copy.copy ( oldPath )
                            newpath.append ( otherAtom )
                            self.paths.append ( newpath )
        # . Get the sorted atomic number and mass data for the level.
        self.levelData = []
        for path in self.paths:
            topAtom = path[-1]
            self.levelData.append ( ( topAtom.atomicNumber, topAtom.mass ) )
        self.levelData.sort    ( )
        self.levelData.reverse ( )

class _CIPDummyAtom:
    """Class to represent a CIP dummy atom."""

    def __init__ ( self, atomicNumber, mass ):
        """Constructor."""
        self.atomicNumber = atomicNumber
        self.mass         = mass

#===============================================================================
# . Testing.
#===============================================================================
if __name__ == "__main__" :
    pass
