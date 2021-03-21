"""Helper functions for building coordinates and other operations.

These functions only use connectivity information. An energy model of a particular type is not required.
"""

import math

from pCore                     import logFile       , \
                                      LogFileActive , \
                                      Selection
from pMolecule                 import System
from pScientific               import PeriodicTable
from pScientific.Geometry3     import Coordinates3  , \
                                      Matrix33      , \
                                      Vector3
from pScientific.Graph         import ConnectedComponentExcludingEdge 
from pScientific.RandomNumbers import RandomNumberGenerator

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Hydrogen building.
# . Default heavy atom - hydrogen atom bond distance (in Angstroms).
_DefaultBondDistance     = 1.0

# . Default angles for different coordinations (in degrees).
_LinearAngle             = 180.0
_TetrahedralAngle        = math.degrees ( 2.0 * math.asin ( math.sqrt ( 2.0 ) / math.sqrt ( 3.0 ) ) )
_TrigonalPlanarAngle     = 120.0

# . The default angles for a given coordination.
_CoordinationAngles      = { 1: 0.0, 2: _LinearAngle, 3: _TrigonalPlanarAngle, 4: _TetrahedralAngle }

# . The default plane angle for a given coordination.
_CoordinationPlaneAngles = { 3: 180.0, 4: 180.0 - 0.5 * _TetrahedralAngle }

# . Possible undefined coordinate value and tolerance.
#_UNDEFINEDCOORDINATE = 9999.0
#_UNDEFINEDTOLERANCE  = 0.1

#===================================================================================================================================
# . Build hydrogen coordinates using connectivity information only.
#===================================================================================================================================
def BuildHydrogenCoordinates3FromConnectivity ( system, log = logFile, randomNumberGenerator = None ):
    """Build hydrogen coordinates.

    The coordinates are built using connectivity information only (bonds
    and angles) which means that no account is taken of non-connectivity
    information (such as hydrogen-bonding). These interactions will have
    to be optimized separately using energy minimization or dynamics.

    Note that bonds to other hydrogens are ignored and unbound hydrogens
    or hydrogens linked to more than one heavy atom will not be built."""

    # . Check for a system object.
    if isinstance ( system, System ) and ( system.connectivity is not None ) and ( len ( system.connectivity.bonds ) > 0 ):

        # . Check whether there are undefined coordinates.
        coordinates3     = system.coordinates3
        numberUndefined0 = coordinates3.numberUndefined
        if numberUndefined0 > 0:

            # . Initialization.
            direction   = Vector3.Null ( )
            neighbors   = system.connectivity.adjacentNodes
            nodeToIndex = system.connectivity.nodeIndices
            if randomNumberGenerator is None:
                randomNumberGenerator = RandomNumberGenerator.WithRandomSeed ( )

            # . Loop over heavy atoms with defined coordinates.
            for ( c, atom ) in enumerate ( system.atoms ):
                if ( atom.atomicNumber != 1 ) and ( c not in coordinates3.undefined ):

                    # . Initialization.
                    built       = []
                    builtH      = []
                    toBuild     = []
                    unbuildable = []

                    # . Loop over the connected atoms.
                    others = neighbors[atom]
                    for other in others:
                        i      = nodeToIndex[other]
                        QBUILT = ( i not in coordinates3.undefined )
                        # . Hydrogens.
                        if other.atomicNumber == 1:
                            if   QBUILT                       : builtH.append      ( i )
                            elif len ( neighbors[other] ) == 1: toBuild.append     ( i )
                            else                              : unbuildable.append ( i )
                        # . Other atoms.
                        else:
                            if QBUILT: built.append       ( i )
                            else:      unbuildable.append ( i )

                    # . Skip this atom if the number of connections is greater than four, there are no hydrogens to build or there are unbuildable atoms.
                    if ( len ( others ) > 4 ) or ( len ( toBuild ) == 0 ) or ( len ( unbuildable ) > 0 ): continue

                    # . Order the lists and put built hydrogens after built heavy atoms as it is reasoned that heavy atom coordinates will be more reliable.
                    built.sort   ( )
                    builtH.sort  ( )
                    toBuild.sort ( )
                    built += builtH

                    # . Get coordination data for the center.
                    nConnections = len ( built ) + len ( toBuild )
                    bondLength   = PeriodicTable[atom.atomicNumber].GetSingleBondDistance ( 1 )
                    angle        = PeriodicTable[atom.atomicNumber].GetCoordinationAngle ( nConnections )
                    if angle      is None: angle      = _CoordinationAngles.get ( nConnections, 0.0 )
                    if bondLength is None: bondLength = _DefaultBondDistance
                    planeAngle = _CoordinationPlaneAngles.get ( nConnections, 0.0 )
                    # . Build the hydrogens.
                    while len ( toBuild ) > 0:

                        # . Get the hydrogen index.
                        h = toBuild.pop  ( 0 )

                        # . Build according to the number of built connected atoms.
                        nbuilt = len ( built )

                        # . Get a random normalized vector.
                        if ( nbuilt == 0 ) or ( nbuilt == 1 ):
                            for i in range ( 3 ): direction[i] = 2.0 * ( randomNumberGenerator.NextReal ( ) - 0.5 )
                            direction.Normalize ( )

                        # . Put the hydrogen in a random direction from the center.
                        # . Works for all cases.
                        if   nbuilt == 0: coordinates3.BuildPointFromDistance ( h, c, bondLength, direction )

                        # . Put the hydrogen at the correct angle from the center and built atom but in a random plane.
                        # . Works for all cases given correct choice of angle.
                        elif nbuilt == 1:
                            coordinates3.BuildPointFromDistanceAngle ( h, c, built[0], bondLength, angle, direction )

                        # . Put the hydrogen away from the other built points at an appropriate angle from their plane.
                        # . The sign of the plane angle is arbitrary.
                        # . Works for cases 3, 4, 5 (square pyramidal), 6 with correct choice of planeAngle.
                        elif nbuilt == 2:
                            coordinates3.BuildPointFromDistancePlaneAngle ( h, c, built[0], built[1], bondLength, planeAngle )

                        # . Put the hydrogen using a tetrahedral tripod.
                        # . Only works for tetrahedral coordination.
                        elif nbuilt == 3:
                            coordinates3.BuildPointFromDistanceTetrahedralTripod ( h, c, built[0], built[1], built[2], bondLength )

                        # . Cannot handle valencies greater than 4 for the moment.
                        else: break

                        # . The hydrogen has been built.
                        built.append ( h )
                        coordinates3.FlagCoordinateAsDefined ( h )

            # . Output a summary.
            if LogFileActive ( log ):
                numberToBuild = coordinates3.numberUndefined
                numberBuilt   = ( numberUndefined0 - numberToBuild )
                if   numberBuilt <= 0: log.Paragraph ( "Coordinates for no hydrogens were built." )
                elif numberBuilt == 1: log.Paragraph ( "Coordinates for one hydrogen were built." )
                else:                  log.Paragraph ( "Coordinates for {:d} hydrogens were built.".format ( numberBuilt ) )

#===================================================================================================================================
# . Identify undefined coordinates.
#===================================================================================================================================
def IdentifyUndefinedCoordinates3 ( system, log = logFile, printHeavies = True, printHydrogens = False, usePDBNotation = True ):
    """Identify any undefined coordinates."""
    # . Initialization.
    numberOfHeavies   = 0
    numberOfHydrogens = 0
    # . See if there are undefined coordinates.
    if system.coordinates3.numberUndefined > 0:
        # . Get the indices of the undefined coordinates.
        undefined = system.coordinates3.undefined
        # . Find numbers.
        heavies   = []
        hydrogens = []
        for i in undefined:
            n = system.atoms[i].atomicNumber
            if n == 1: hydrogens.append ( i )
            else:      heavies.append   ( i )
        numberOfHeavies   = len ( heavies   )
        numberOfHydrogens = len ( hydrogens )
        # . Set some options.
        usePDBNotation = usePDBNotation and ( system.sequence is not None )
        # . Print the results.
        if LogFileActive ( log ):
            # . Basic summary.
            log.SummaryOfItems ( [ ( "Undefined Heavy Atoms", "{:d}".format ( numberOfHeavies   ) )   ,
                                     ( "Undefined Hydrogens",   "{:d}".format ( numberOfHydrogens ) ) ] ,
                                     title = "Undefined Coordinates3 Summary" )
            # . Explicit listings.
            for ( doPrinting, indices, tag ) in ( ( printHeavies, heavies, "Heavy" ), ( printHydrogens, hydrogens, "Hydrogen" ) ):
                n = len ( indices )
                if doPrinting and ( n > 0 ):
                    if usePDBNotation: table = log.GetTable ( columns = min (  5, n ) * [ 20 ] )
                    else:              table = log.GetTable ( columns = min ( 10, n ) * [ 10 ] )
                    table.Start ( )
                    table.Title ( "Undefined " + tag + " Atoms" )
                    for i in indices:
                        table.Entry ( system.atoms[i].path )
                    table.Stop ( )
    # . Return.
    return ( numberOfHeavies, numberOfHydrogens )

#===================================================================================================================================
# . Rotate about a dihedral, either by a certain amount or to a target value.
# . Input dihedral values are in degrees!
#===================================================================================================================================
def RotateDihedral ( system, a2, a3, dihedral ):
    """Rotate a dihedral by a certain amount.."""
    # . a2-a3 are the atoms of the bond to be rotated.
    # . Only a3 and the atoms that are joined to it are rotated.
    if dihedral != 0.0:
        # . Initialization.
        coordinates3 = system.coordinates3
        nodeIndices  = system.connectivity.nodeIndices
        ( j, k )     = ( nodeIndices[a2], nodeIndices[a3] )
        # . Find atoms in the a3 isolate but without the a2-a3 bond.
        isolate  = ConnectedComponentExcludingEdge ( system.connectivity, a3, a2 )
        toRotate = Selection.FromIterable ( [ nodeIndices[a] for a in isolate ] )
        # . Do the rotation.
        axis        = coordinates3.Displacement ( j, k )
        axis.Normalize ( )
        rotation    = Matrix33.MakeRotationAboutAxis ( math.radians ( dihedral ), axis )
        translation = Vector3.Null ( )
        for c in range ( 3 ): translation[c] = - coordinates3[j,c]
        coordinates3.Translate ( translation )
        translation.Scale ( -1.0 )
        coordinates3.Rotate ( rotation, selection = toRotate )
        coordinates3.Translate ( translation )

def RotateDihedralToTarget ( system, a1, a2, a3, a4, target ):
    """Rotate a dihedral so that it has a given value."""
    # . a1-a2-a3-a4 are the atoms of the dihedral.
    # . Only a4 and the atoms that are joined to it are rotated.
    nodeIndices    = system.connectivity.nodeIndices
    ( i, j, k, l ) = ( nodeIndices[a1], nodeIndices[a2], nodeIndices[a3], nodeIndices[a4] )
    dihedral       = target - system.coordinates3.Dihedral ( i, j, k, l )
    RotateDihedralToTarget ( system, a2, a3, dihedral )

#===================================================================================================================================
# . Verify fixed atom coordinates.
#===================================================================================================================================
def VerifyFixedAtomCoordinates ( system, trial, log = logFile, reference = None, tolerance = 1.0e-03 ):
    """Verify that the coordinates of fixed atoms in a trial set are the same as in a reference set."""
    # . Initialization.
    fixedAtoms = None
    if system.freeAtoms is not None: fixedAtoms = system.freeAtoms.Complement ( upperBound = len ( system.atoms ) )
    isOK       = False
    if reference is None: reference  = system.coordinates3
    # . Basic checks.
    if   reference is None: message = "Missing reference coordinate set."
    elif trial     is None: message = "Missing trial coordinate set."
    elif ( trail.rows != len ( system.atoms ) ) or ( len ( reference ) != len ( trial ) ):
        message = "Coordinate sets of incompatible size."
    elif fixedAtoms is None:
        isOK    = True
        message = "No fixed atoms are defined."
    # . Check the coordinates.
    else:
        set1 = Coordinates3 ( len ( fixedAtoms ) ) ; set1.Gather ( reference, fixedAtoms )
        set2 = Coordinates3 ( len ( fixedAtoms ) ) ; set2.Gather ( other    , fixedAtoms )
        set1.Add ( set2, scale = -1.0 )
        biggest = set1.AbsoluteMaximum ( )
        if biggest >= tolerance:
            message = "The maximum difference in the fixed atom coordinates, {:.6g}, is greater than the tolerance, {:.6g}.".format ( biggest, tolerance )
        else:
            isOK    = True
            message = "The fixed atom coordinates are compatible."
    # . Printing.
    if LogFileActive ( log ): log.Paragraph ( message )
    # . Finish up.
    return isOK

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
