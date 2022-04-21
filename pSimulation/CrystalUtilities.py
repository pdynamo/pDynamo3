"""Utilities for dealing with crystals."""

import math

from  pCore                 import Clone                      , \
                                   logFile                    , \
                                   LogFileActive              , \
                                   SelectionContainer
from  pMolecule             import System
from  pScientific.Geometry3 import CrossPairList_FromDoubleCoordinates3
from  pScientific.Symmetry  import CrystalSystemTriclinic     , \
                                   PeriodicBoundaryConditions , \
                                   SymmetryError
from .EditMergePrune        import MergeByAtom

#===================================================================================================================================
# . Analyze the transformations of a crystal.
#===================================================================================================================================
def CrystalAnalyzeTransformations ( system, log = logFile ):
    """Analyze the transformations for a crystal.

    Transformations must be either proper or improper rotations and have inverses.

    An inverse check needs to be added.
    """

    # . Basic checks.
    if LogFileActive ( log ) and isinstance ( system, System ) and ( system.symmetry is not None ):

        # . Get the lattice matrices for the system.
        sp       = system.symmetryParameters
        M        = sp.H
        inverseH = sp.inverseH

        # . Loop over the transformations.
        for ( i, t3 ) in enumerate ( system.symmetry.transformations ):

            # . Output the transformation.
            newT3 = Clone ( t3 )
            newT3.Orthogonalize ( M, inverseH )
            newT3.Print ( log = log, title = "Transformation {:d}".format ( i ) )

            # . Check for a rotation of some sort.
            if   newT3.rotation.isProperRotation  : log.Paragraph ( "Transformation is a proper rotation."    )
            elif newT3.rotation.isImproperRotation: log.Paragraph ( "Transformation is an improper rotation." )
            else: log.Paragraph ( "Transformation is neither a proper nor an improper rotation." )

#===================================================================================================================================
# . Center the coordinates of a system by putting them inside the primary image so that their fractional coordinates are in the
# . range [0,1]. Centering can be done by atom, by isolate or by MM isolate (the default).
#
# . Centering is done with respect to the origin of the coordinate system, but it can also be done with respect to a selected set
# . of indices (in toCenter).
#
# . Centering of each isolate is done with respect to its center of geometry.
#
# . If |selection| is present, only selected atoms or isolates with selected atoms will be centered.
#
# . The return values of the function are the centered coordinates and the translations effected (i.e. the difference between the
# . new and old coordinates).
#===================================================================================================================================
def CrystalCenterCoordinates ( system , log              = logFile       ,
                                        centerFixedAtoms = False         ,
                                        mode             = "byMMIsolate" ,
                                        selection        = None          ,
                                        toCenter         = None          ):
    """Center the coordinates of a system in the primary image."""
    coordinates3 = None
    translations = None
    if isinstance ( system, System ) and \
       ( system.symmetry           is not None ) and \
       ( system.coordinates3       is not None ) and \
       ( system.symmetryParameters is not None ):
        if ( mode not in ( "byAtom", "byIsolate", "byMMIsolate" ) ): raise SymmetryError ( "Unknown centering mode: " + mode + "." )
        # . Get a new set of coordinates to work with.
        coordinates3 = Clone ( system.coordinates3 )
        # . Is there a set of atoms to be centered in the box?
        if toCenter is not None:
            translation = coordinates3.Center ( selection = toCenter )
            translation.Scale ( -1.0 )
            for i in range ( 3 ):
                for j in range ( 3 ): translation[i] += ( 0.5 * system.symmetryParameters.H[i,j] )
            coordinates3.Translate ( translation )
            if selection is None:
                selection = toCenter.Complement ( len ( system.atoms ) )
            else:
                selection = Clone ( selection )
                selection.Exclude ( toCenter  )
        # . Center the coordinates.
        if ( mode == "byAtom" ):
            if ( not centerFixedAtoms ) and ( system.freeAtoms is not None ):
                selection = Selection.FromIterable ( set ( selection ).intersection ( set ( system.freeAtoms ) ) )
            system.symmetryParameters.CenterCoordinatesByAtom ( coordinates3, selection = selection )
        else:
            freeIsolates = None
            isolates     = None
            if mode == "byIsolate":
                isolates = Clone ( system.connectivity.isolateIndices )
            elif system.mmModel is not None :
                n          = len ( system.atoms )
                exclusions = getattr ( system.mmModel, "exclusions", None )
                if exclusions is None: isolates = SelectionContainer.FromCapacity ( n )
                else:                  isolates = exclusions.GetConnectedComponents ( upperBound = n )
            if isolates is None: raise SymmetryError ( "Unable to determine isolates for the system." )
            if system.qcModel is not None: # . QC atoms are kept together.
                toFuse = isolates.MakeMembershipFlags ( system.qcState.qcAtoms )
                isolates.FuseItems ( toFuse )
            if system.freeAtoms is not None:
                if centerFixedAtoms: # . Isolates with fixed atoms are kept together.
                    toFuse = isolates.MakeMembershipFlags ( system.freeAtoms, andTest = True )
                    toFuse.Negate ( )
                    isolates.FuseItems ( toFuse )
                else: # . Isolates with fixed atoms do not move.
                    freeIsolates = isolates.MakeMembershipFlags ( system.freeAtoms, andTest = True )
            system.symmetryParameters.CenterCoordinates3ByFreeIsolate ( isolates, freeIsolates, coordinates3 )
        # . Get the translations.
        translations = Clone ( coordinates3 )
        translations.Add ( system.coordinates3, scale = -1.0 )
    return ( coordinates3, translations )

#===================================================================================================================================
# . Non-P1 to P1 with an arbitrary number of cells in each direction.
#===================================================================================================================================
def CrystalExpandToP1 ( system, arange = range ( 1 ), brange = range ( 1 ), crange = range ( 1 ), imposeTriclinic = False, log = logFile ):
    """Create a system with P1 symmetry from one with arbitrary crystal symmetry.

    |system| is the input system.

    |arange| is the range of unit cells to build in the a direction.
    |brange| is the range of unit cells to build in the a direction.
    |crange| is the range of unit cells to build in the a direction.

    The return value is the new system but note that:
    |None|   is returned if |system| does not have crystal symmetry.
    |system| is returned if it is already P1.

    Although the system is converted to P1, the original crystal class is retained.

    Other options to include are:
    |QTRICLINIC| for selecting a triclinic crystal class for the output system.
    """
    result = None

    # . Get the number of cells.
    ncells = len ( arange ) * len ( brange ) * len ( crange )

    # . Basic checks.
    if isinstance ( system, System ) and ( system.symmetry is not None ):

        # . The system is already P1 and only one cell is required.
        if system.symmetry.transformations.isIdentity and ( ncells == 1 ):
            result = system

        # . Build the new system.
        else:

            # . Get an alias for the system's symmetry parameters.
            sp = system.symmetryParameters

            # . Define the lattice matrix H and its inverse.
            M        = sp.H
            inverseH = sp.inverseH

            # . Create the images of the original system by building over each cell in turn.
            images = []
            for a in arange:
                for b in brange:
                    for c in crange:
                        for t3 in system.symmetry.transformations:
                            # . Create the transformation.
                            newT3 = Clone ( t3 )
                            newT3.translation[0] += a
                            newT3.translation[1] += b
                            newT3.translation[2] += c
                            newT3.Orthogonalize ( M, inverseH )
                            # . Create the image.
                            image = Clone ( system )
                            image.coordinates3.Transform ( newT3 )
                            images.append ( image )

            # . Get the new box lengths.
            a = len ( arange ) * sp.a
            b = len ( brange ) * sp.b
            c = len ( crange ) * sp.c

            # . Create the new system with the correct symmetry.
            if imposeTriclinic: crystalSystem = CrystalSystemTriclinic ( )
            else:               crystalSystem = system.symmetry.crystalSystem
            result                    = MergeByAtom ( images )
            result.symmetry           = PeriodicBoundaryConditions.WithCrystalSystem ( crystalSystem )
            result.symmetryParameters = result.symmetry.MakeSymmetryParameters ( a = a, b = b, c = c, alpha = sp.alpha, beta = sp.beta, gamma = sp.gamma )
            if system.label is None: result.label = "P1 Crystal System"
            else:                    result.label = system.label + " - P1 Crystal Symmetry"
            result.Summary ( log = log )

    return result

#===================================================================================================================================
# . Get lists of possible bonds between the primary and secondary images.
# . It may be advisable to center the coordinates first (or use a large
# . range for the image search).
#===================================================================================================================================
def CrystalGetImageBondPairs ( system, arange = range ( -1, 2 ), brange = range ( -1, 2 ), crange = range ( -1, 2 ), radii = None, safety = 0.45 ):
    """Create pairlists of interactions between the primary and secondary images."""

    # . Initialization.
    results = None

    # . Basic checks.
    if isinstance ( system, System ) and ( system.symmetry is not None ) and hasattr ( system.symmetry, "transformations" ):

        # . Need to have centering option!

        # . Get an alias for the system's symmetry parameters.
        sp = system.symmetryParameters

        # . Define the lattice matrix H and its inverse.
        M        = sp.H
        inverseH = sp.inverseH

        # . Create the images of the original system by building over each cell in turn.
        results = []
        for a in arange:
            for b in brange:
                for c in crange:
                    for ( t, t3 ) in enumerate ( system.symmetry.transformations ):
                        # . Create the transformation.
                        newT3 = Clone ( t3 )
                        newT3.translation[0] += a
                        newT3.translation[1] += b
                        newT3.translation[2] += c
                        newT3.Orthogonalize ( M, inverseH )
                        # . Skip the identity transformation.
                        if not newT3.isIdentity:
                            # . Create the image coordinates.
                            icoordinates3 = Clone ( system.coordinates3 )
                            icoordinates3.Transform ( newT3 )
                            # . Create the pair list.
                            pairlist = CrossPairList_FromDoubleCoordinates3 ( system.coordinates3, icoordinates3, radii1 = radii, radii2 = radii, safety = safety )
                            # . Calculate and save the pairs.
                            if ( pairlist is not None ) and ( len ( pairlist ) > 0 ):
                                pairs = []
                                for ( i, j ) in pairlist:
                                    dx = system.coordinates3[i,0] - icoordinates3[j,0]
                                    dy = system.coordinates3[i,1] - icoordinates3[j,1]
                                    dz = system.coordinates3[i,2] - icoordinates3[j,2]
                                    d  = math.sqrt ( dx * dx + dy * dy + dz * dz )
                                    pairs.append ( ( i, j, d ) )
                                results.append ( ( t, a, b, c, pairs ) )
    return results

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
