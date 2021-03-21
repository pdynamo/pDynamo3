"""Script to calculate the ion mobilities.

Based on the article:

A. A. Shvartsburg and M. F. Jarrold
"An Exact Hard Spheres Scattering Model for the Mobilities of Polyatomic Ions"
Chem. Phys. Lett. 1996, 261, 86.
"""
# . This could be made faster by using the same objects in the loops.

import math

from pCore                     import Clone                 , \
                                      LogFileActive         , \
                                      logFile
from pScientific               import Constants             , \
                                      PeriodicTable         , \
                                      Units
from pScientific.Arrays        import Array
from pScientific.Geometry3     import Matrix33              , \
                                      Vector3
from pScientific.RandomNumbers import RandomNumberGenerator

#===============================================================================
# . Parameters.
#===============================================================================
# . Hard-sphere elemental radii in metres.
_HSRadii = { 1 : 2.200e-10 ,
             6 : 2.700e-10 ,
             7 : 2.700e-10 ,
             8 : 2.700e-10 ,
            11 : 2.853e-10 ,
            14 : 2.950e-10 ,
            15 : 3.500e-10 ,
            26 : 3.500e-10 }

# . Mass constant without the mass term and assuming a charge of one and standard He density.
_MassConstant = ( math.sqrt ( 18.0 * math.pi ) / 16.0 ) * math.sqrt ( Constants.Avogadro_Number * 1000.0 ) * \
                                         ( Constants.Electron_Charge / math.sqrt ( Constants.Boltzmann ) ) * \
                                         ( Constants.Molar_Ideal_Gas_Volume / Constants.Avogadro_Number )

#===============================================================================
# . Public functions.
#===============================================================================
def HardSphereIonMobilities ( molecule                        ,
                              log                   = logFile ,
                              nreflections          = 30      ,
                              ntrajectories         = 600000  ,
                              randomNumberGenerator = None    ,
                              temperature           = 298.0   ):
    """Calculate ion mobilities with a hard-sphere model."""

    # . Get the atom data.
    hsRadii   = _GetHardSphereRadii ( molecule.atoms )
    masses    = Array.FromIterable ( [ atom.mass for atom in molecule.atoms ] )
    totalmass = masses.Sum ( )

    # . Get initial coordinates, move to center of mass and convert to metres.
    xyz0 = Clone ( molecule.coordinates3 )
    xyz0.TranslateToCenter ( weights = masses )
    xyz0.Scale ( 1.0e-10 )

    # . Get the mass constant.
    massHe       = PeriodicTable[2].mass
    massConstant = _MassConstant * math.sqrt ( ( 1.0 / massHe ) + ( 1.0 / totalmass ) )

    # . Get the random number generator.
    if randomNumberGenerator is None: randomNumberGenerator = RandomNumberGenerator.WithRandomSeed ( )
    rotation = Matrix33.Null ( )

    # . Initialize some calculation variables.
    cof          = Array.WithExtent ( nreflections ) ; cof.Set  ( 0.0 )
    crof         = Array.WithExtent ( nreflections ) ; crof.Set ( 0.0 )
    crb          = 0.0
    mreflections = 0

    # . Loop over the trajectories.
    for it in range ( ntrajectories ):

        # . Randomly rotate the coordinate set.
        rotation.RandomRotation ( randomNumberGenerator )
        xyz = Clone ( xyz0 )
        xyz.Rotate ( rotation )

        # . Loop over the collisions.
        hasCollision = False
        for ir in range ( nreflections ):

            # . Initial collision - at a random point in the yz plane along the x-axis.
            if ir == 0:
                ( origin, extents ) = xyz.EnclosingOrthorhombicBox ( radii = hsRadii )
                yzarea = extents[1] * extents[2]
                yc     = origin[1] + extents[1] * randomNumberGenerator.NextReal ( )
                zc     = origin[2] + extents[2] * randomNumberGenerator.NextReal ( )
                xaxis  = Vector3.WithValues ( 1.0, 0.0, 0.0 )
            # . Subsequent collisions - always along the x-axis.
            else:
                yc = 0.0
                zc = 0.0

            # . Initialization.
            ic = -1                     # . The index of the colliding particle.
            xc = origin[0] + extents[0] # . The largest x-coordinate.

            # . Loop over particles.
            for ( i, h ) in enumerate ( hsRadii ):
                # . After the first collision only x-values > 0 are allowed.
                if ( ir == 0 ) or ( xyz[i,0] > 1.0e-16 ):
                    # . yd and zd are the coordinates of the impact points for the ith atom
                    # . with respect to its own coordinates (if such a point exists).
                    # . dev is the impact parameter.
                    h2  = h * h
                    y   = yc - xyz[i,1]
                    z   = zc - xyz[i,2]
                    yz2 = y * y + z * z
                    # . If there is a collision with the ith atom, check to see if it occurs before previous collisions.
                    if yz2 < h2:
                        x = xyz[i,0] - math.sqrt ( h2 - yz2 )
                        if x < xc:
                            xc = x
                            ic = i

            # . Check mreflections.
            if ir >= mreflections: mreflections = ir + 1

            # . There was a collision.
            if ic >= 0:
                hasCollision = True
                # . Translate the coordinates so that the collision point is at the origin.
                xyz.Translate ( Vector3.WithValues ( -xc, -yc, -zc ) )
                # . Rotate the coordinates so that the outgoing vector is along the x-axis.
                h = Vector3.Null ( )
                xyz[ic].CopyTo ( h )                                   # . Normalized vector from the collision point to the ic-th atom.
                h.Normalize ( tolerance = 1.0e-20 )
                axis  = Vector3.WithValues ( 0.0, h[2], -h[1] )       # . Normalized axis of rotation.
                axis.Normalize ( tolerance = 1.0e-20 )
                alpha = math.pi - 2.0 * math.acos ( h[0] )             # . Angle of rotation.
                rotation.RotationAboutAxis ( alpha, axis )
                xyz.Rotate ( rotation )
                rotation.ApplyTo ( xaxis )
                # . Calculate the cosine of the angle between the incoming vector and the normal to a plane,
                # . the reflection from which would be equivalent to the accumulated reflection.
                # . This is equal to h[0] when ir = 0.
                cof[ir] = math.cos ( 0.5 * ( math.pi - math.acos ( xaxis[0] ) ) )
                # . Check outgoing.
                # . Get the outgoing vector (the ingoing vector is always [1,0,0]).
                out = Vector3.WithValues ( 1.0 - 2.0 * h[0] * h[0], - 2.0 * h[0] * h[1], - 2.0 * h[0] * h[2] )
                rotation.ApplyTo ( out )
                out[0] -= 1.0
                if out.Norm2 ( ) > 1.0e-6: print ( "Invalid Rotation: {:10.3f} {:10.3f} {:10.3f}.".format ( out[0], out[1], out[2] ) )
            # . There was no collision.
            else:
                # . Top up the remaining elements of cof with the last valid value of cof.
                if ir == 0: t = 0.0
                else:       t = cof[ir-1]
                for i in range ( ir, nreflections ): cof[i] = t
                # . Exit.
                break

        # . End of collisions.
        # . Projection approximation.
        if hasCollision: crb += yzarea
        # . Hard-sphere approximation.
        for ir in range ( nreflections ):
            crof[ir] += yzarea * cof[ir] * cof[ir]

    # . End of trajectories.
    crof.Scale ( 2.0 / float ( ntrajectories ) )
    pacs  = crb / float ( ntrajectories )
    pamob = massConstant / ( pacs * math.sqrt ( temperature ) )
    hscs  = crof[mreflections-1]
    hsmob = massConstant / ( hscs * math.sqrt ( temperature ) )

    # . Output results.
    if LogFileActive ( log ):
        items = [ ( "MC Trajectories",  "{:d}"  .format ( ntrajectories  ) ) ,
                  ( "Reflection Limit", "{:d}"  .format ( nreflections   ) ) ,
                  ( "PA Mobility",      "{:.4g}".format ( pamob          ) ) ,
                  ( "PA Cross-Section", "{:.4g}".format ( pacs * 1.0e+20 ) ) ,
                  ( "HS Mobility",      "{:.4g}".format ( hsmob          ) ) ,
                  ( "HS Cross-Section", "{:.4g}".format ( hscs * 1.0e+20 ) ) ,
                  ( "Max. Reflections", "{:d}"  .format ( mreflections   ) ) ]
        log.SummaryOfItems ( items, title = "Hard-Sphere Ion Mobilities" )

    # . Finish up.
    results = { "MC Trajectories"     : ntrajectories  ,
                "Reflection Limit"    : nreflections   ,
                "PA Mobility"         : pamob          ,
                "PA Cross-Section"    : pacs * 1.0e+20 ,
                "HS Mobility"         : hsmob          ,
                "HS Cross-Section"    : hscs * 1.0e+20 ,
                "Maximum Reflections" : mreflections   }
    return results

#===============================================================================
# . Private functions.
#===============================================================================
def _GetHardSphereRadii ( atoms ):
    """Get the hard-sphere radii for the atoms."""
    hsRadii       = Array.WithExtent ( len ( atoms ) ) ; hsRadii.Set ( 0.0 )
    atomicNumbers = [ atom.atomicNumber for atom in atoms ]
    for ( i, atomicNumber ) in enumerate ( atomicNumbers ):
        try:    hsRadii[i] = _HSRadii[atomicNumber]
        except: raise KeyError ( "Hard-sphere radius for element " + PeriodicTable.Symbol ( atomicNumber ) + " unknown." )
    return hsRadii

#===============================================================================
# . Testing.
#===============================================================================
if __name__ == "__main__" :
    pass
