"""Utility functions for analyzing and manipulating trajectories."""

import math, os, os.path

from pBabel                import ExportSystem       , \
                                  ExportTrajectory   , \
                                  ImportCoordinates3 , \
                                  ImportTrajectory
from pCore                 import Clone              , \
                                  logFile            , \
                                  LogFileActive      , \
                                  Selection
from pScientific.Arrays    import Array              , \
                                  StorageType
from pScientific.Geometry3 import Coordinates3
from pScientific.Symmetry  import CrystalSystemCubic

# . The trajectories to AveragePositions, CoordinateFluctuations and CovarianceMatrix should have had their rotation/translational motion removed.

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Safety factor for estimating upper bound for RDF calculation.
_DefaultSafety = 1.0

#===================================================================================================================================
# . Helper functions for argument processing.
#===================================================================================================================================
def _GetSystemAndTrajectoriesFromArguments ( path, system, format = None ):
    """Get a system and its trajectories from an argument list."""
    if   isinstance ( path,   str           ): paths = [ path ]
    elif isinstance ( path, ( list, tuple ) ): paths =   path
    else: raise TypeError ( "Invalid |path| argument." )
    return [ ImportTrajectory ( path, system, format = format ) for path in paths ]

#===================================================================================================================================
# . Calculate average positions.
#===================================================================================================================================
def AveragePositions ( path, system, format = None, selection = None ):
    """Calculate the average positions for selected particles."""
    # . Initialization.
    positions    = None
    trajectories = _GetSystemAndTrajectoriesFromArguments ( path, system, format = format )
    # . Get the selection (or all particles otherwise).
    if selection is None: selection = Selection.FromIterable ( range ( len ( system.atoms ) ) )
    # . Get the size of the problem.
    n = len ( selection )
    if n > 0:
        # . Allocate space.
        positions = Coordinates3.WithExtent ( n )
        positions.Set ( 0.0 )
        # . Loop over trajectory frames.
        numberFrames = 0
        for trajectory in trajectories:
            trajectory.ReadHeader ( )
            while trajectory.RestoreOwnerData ( ):
                frame = system.coordinates3
                for ( p, i ) in enumerate ( selection ):
                    positions[p,0] += frame[i,0]
                    positions[p,1] += frame[i,1]
                    positions[p,2] += frame[i,2]
            trajectory.ReadFooter ( )
            trajectory.Close ( )
            numberFrames += len ( trajectory )
        # . Scale.
        if numberFrames > 0: positions.Scale ( 1.0 / float ( numberFrames ) )
    return positions

#===================================================================================================================================
# . Coordinate fluctuations.
#===================================================================================================================================
def CoordinateFluctuations ( path, system, anisotropic = False, asBFactors = False, averagePositions = None, format = None, selection = None ):
    """Calculate the coordinate fluctuations for selected particles."""
    # . Initialization.
    fluctuations = None
    if averagePositions is None:
        averagePositions = AveragePositions ( path, system, selection = selection )
    if selection is None:
        selection = Selection.FromIterable ( range ( len ( system.atoms ) ) )
    # . Continue processing.
    if ( averagePositions is not None ) and ( len ( selection ) > 0 ):
        # . Activate trajectories.
        trajectories = _GetSystemAndTrajectoriesFromArguments ( path, system, format = format )
        # . Allocate space.
        n = len ( selection )
        displacement = Coordinates3.WithExtent ( n )
        if anisotropic: fluctuations = Array.WithExtents ( n, 6 )
        else:           fluctuations = Array.WithExtent  ( n    )
        displacement.Set ( 0.0 )
        fluctuations.Set ( 0.0 )
        # . Loop over trajectory frames.
        numberFrames = 0
        for trajectory in trajectories:
            trajectory.ReadHeader ( )
            while trajectory.RestoreOwnerData ( ):
                frame = system.coordinates3
                if anisotropic:
                    for ( p, i ) in enumerate ( selection ):
                        dx = frame[i,0] - averagePositions[p,0]
                        dy = frame[i,1] - averagePositions[p,1]
                        dz = frame[i,2] - averagePositions[p,2]
                        fluctuations[p,0] += dx * dx
                        fluctuations[p,1] += dy * dx
                        fluctuations[p,2] += dy * dy
                        fluctuations[p,3] += dz * dx
                        fluctuations[p,4] += dz * dy
                        fluctuations[p,5] += dz * dz
                else:
                    for ( p, i ) in enumerate ( selection ):
                        fluctuations[p] += ( ( frame[i,0] - averagePositions[p,0] )**2 + \
                                             ( frame[i,1] - averagePositions[p,1] )**2 + \
                                             ( frame[i,2] - averagePositions[p,2] )**2 )
            trajectory.ReadFooter ( )
            trajectory.Close ( )
            numberFrames += len ( trajectory )
        # . Scale.
        if numberFrames > 0: fluctuations.Scale ( 1.0 / float ( numberFrames ) )
        # . Convert to B-factors if necessary.
        if asBFactors:
            conversionFactor = 8.0 * math.pi**2
            if not anisotropic: conversionFactor /= 3.0
            fluctuations.Scale ( conversionFactor )
    # . Finish up.
    return fluctuations

#===================================================================================================================================
# . Calculate the covariance matrix.
#===================================================================================================================================
def CovarianceMatrix ( path, system, averagePositions = None, format = None, selection = None ):
    """Calculate the covariance matrix for selected particles."""
    # . Initialization.
    covariance   = None
    if averagePositions is None:
        averagePositions = AveragePositions ( path, system, selection = selection )
    if selection is None:
        selection = Selection.FromIterable ( range ( len ( system.atoms ) ) )
    # . Continue processing.
    if ( averagePositions is not None ) and ( len ( selection ) > 0 ):
        # . Activate trajectories.
        trajectories = _GetSystemAndTrajectoriesFromArguments ( path, system, format = format )
        # . Allocate space.
        n = 3 * len ( selection )
        covariance   = Array.WithExtent ( n , storageType = StorageType.Symmetric ) ; covariance.Set   ( 0.0 )
        displacement = Array.WithExtent ( n ) ; displacement.Set ( 0.0 )
        # . Loop over trajectory frames.
        numberFrames = 0
        for trajectory in trajectories:
            trajectory.ReadHeader ( )
            while trajectory.RestoreOwnerData ( ):
                frame = system.coordinates3
                for ( p, i ) in enumerate ( selection ):
                    displacement[3*p  ] = frame[i,0] - averagePositions[p,0]
                    displacement[3*p+1] = frame[i,1] - averagePositions[p,1]
                    displacement[3*p+2] = frame[i,2] - averagePositions[p,2]
                for i in range ( n ):
                    dI = displacement[i]
                    for j in range ( i + 1 ):
                        covariance[i,j] += ( dI * displacement[j] )
            trajectory.ReadFooter ( )
            trajectory.Close ( )
            numberFrames += len ( trajectory )
        # . Scale.
        if numberFrames > 0: covariance.Scale ( 1.0 / float ( numberFrames ) )
    return covariance

#===================================================================================================================================
# . Duplicate a trajectory. Also useful for converting formats.
#===================================================================================================================================
def Duplicate ( inPath, outPath, system, inFormat = None, outFormat = None ):
    """Duplicate a trajectory."""
    inTrajectory  = ImportTrajectory ( inPath , system, format = inFormat  )
    outTrajectory = ExportTrajectory ( outPath, system, format = outFormat )
    inTrajectory.ReadHeader   ( )
    outTrajectory.WriteHeader ( )
    while inTrajectory.RestoreOwnerData ( ): outTrajectory.WriteOwnerData ( )
    inTrajectory.ReadFooter   ( )
    outTrajectory.WriteFooter ( )
    inTrajectory.Close  ( )
    outTrajectory.Close ( )

#===================================================================================================================================
# . Expand a trajectory.
#===================================================================================================================================
# . Rotation and translation should have been removed from the input trajectory if necessary.
def ExpandByLinearInterpolation ( inPath, outPath, system, points, inFormat = None, outFormat = None ):
    """Expand a trajectory by adding linearly interpolated points between the existing points."""
    # . Initialization.
    if points < 1: points = 1
    inTrajectory  = ImportTrajectory ( inPath , system, format = inFormat  )
    outTrajectory = ExportTrajectory ( outPath, system, format = outFormat )
    inTrajectory.ReadHeader   ( )
    outTrajectory.WriteHeader ( )
    # . Save system coordinates.
    coordinates3 = system.coordinates3
    saved3       = Clone ( coordinates3 )
    # . Create the start and step array.
    start = None
    step  = Coordinates3.WithExtent ( len ( system.atoms ) ) ; step.Set ( 0.0 )
    # . Loop over structures.
    while inTrajectory.RestoreOwnerData ( ):
        if start is None:
            start = Clone ( coordinates3 )
        else:
            coordinates3.CopyTo ( step )
            step.Add ( start, scale = -1.0 )
            step.Scale ( 1.0 / float ( points + 1 ) )
            start.CopyTo ( coordinates3 )
            for p in range ( points ):
                coordinates3.Add ( step )
                outTrajectory.WriteOwnerData ( )
            coordinates3.Add ( step )
        outTrajectory.WriteOwnerData ( )
    inTrajectory.ReadFooter   ( )
    outTrajectory.WriteFooter ( )
    inTrajectory.Close  ( )
    outTrajectory.Close ( )
    # . Finish up.
    saved3.CopyTo ( coordinates3 )

#===================================================================================================================================
# . Make a trajectory from a set of coordinate files.
#===================================================================================================================================
def FromCoordinateFiles ( inPaths, outPath, system, inFormat = None, outFormat = None ):
    """Make a trajectory from a list of coordinate files."""
    outTrajectory = ExportTrajectory ( outPath, system, format = outFormat )
    outTrajectory.WriteHeader ( )
    for inPath in inPaths:
        system.coordinates3 = ImportCoordinates3 ( inPath, format = inFormat )
        outTrajectory.WriteOwnerData ( )
    outTrajectory.WriteFooter ( )
    outTrajectory.Close ( )

#===================================================================================================================================
# . Make a trajectory by linearly interpolating between two points.
#===================================================================================================================================
# . The first and last structures should be superimposed before entry if required.
def MakeByLinearInterpolation ( outPath, system, points, first, last, outFormat = None ):
    """Make a trajectory by linearly interpolating between two coordinate files."""
    # . Initialization.
    if points < 2: points = 2
    # . Save system coordinates.
    coordinates3 = system.coordinates3
    saved3       = Clone ( coordinates3 )
    # . Find coordinate step.
    if first is not coordinates3: first.CopyTo ( coordinates3 )
    step = Clone ( last )
    step.Add ( first, scale = -1.0 )
    step.Scale ( 1.0 / float ( points - 1 ) )
    # . Make the trajectory.
    outTrajectory = ExportTrajectory ( outPath, system, format = outFormat )
    outTrajectory.WriteHeader ( )
    for p in range ( points ):
        if p > 0: coordinates3.Add ( step )
        outTrajectory.WriteOwnerData ( )
    outTrajectory.WriteFooter ( )
    outTrajectory.Close ( )
    # . Finish up.
    saved3.CopyTo ( coordinates3 )

#===================================================================================================================================
# . Calculate a radial distribution function.
#===================================================================================================================================
def RadialDistributionFunction ( path                        ,
                                 system                      ,
                                 bins       = 100            ,
                                 format     = None           ,
                                 log        = logFile        ,
                                 safety     = _DefaultSafety ,
                                 selection1 = None           ,
                                 selection2 = None           ,
                                 upper      = None           ,
                                 volume     = None           ):
    """Calculate a radial distribution function."""
    # . Initialization.
    trajectories = _GetSystemAndTrajectoriesFromArguments ( path, system, format = format )
    if selection1 is None: selection1 = Selection.FromIterable ( range ( len ( system.atoms ) ) )
    if selection2 is None: selection2 = selection1
    # . Get the numbers of particles.
    np1 = len ( selection1 )
    np2 = len ( selection2 )
    # . Check for symmetry.
    hasSymmetry = ( volume is None )
    if hasSymmetry:
        try:
            cc = system.symmetry.crystalSystem
            if not isinstance ( cc, CrystalSystemCubic ): raise
        except:
            raise ValueError ( "System does not have cubic symmetry." )
    # . Estimate an upper bound.
    if upper is None:
        try:    upper = ( system.symmetryParameters.a - safety ) / 2.0
        except: raise ValueError ( "Please supply an upper bound for the RDF calculation." )
    # . Initialization.
    lower     = 0.0
    distances = []
    histogram = [ 0 for i in range ( bins ) ]
    rdf       = []
    # . Calculate the width of each bin.
    width = ( upper - lower ) / float ( bins )
    v     = 0.0
    # . Loop over trajectory frames.
    numberFrames = 0
    for trajectory in trajectories:
        trajectory.ReadHeader ( )
        while trajectory.RestoreOwnerData ( ):
            # . Get the frame.
            frame         = system.coordinates3
            numberFrames += 1
            # . Deal with symmetry.
            if hasSymmetry:
                a  = system.symmetryParameters.a
                v += a**3
                if upper > 0.5 * a: raise ValueError ( "Box does not satisfy minimum image convention." )
            # . Bin the distances.
            # . This loop involves duplicate work for self rdfs.
            for i in selection1:
                for j in selection2:
                    if i != j:
                        dr = frame.Displacement ( i, j )
                        if hasSymmetry:
                            for c in range ( 3 ):
                                x      = a * round ( dr[c] / a )
                                dr[c] -= x
                        r = dr.Norm2 ( )
                        if ( r >= lower ) and ( r < upper ):
                            b = int ( ( r - lower ) / width )
                            histogram[b] += 1
        # . Finish up.
        trajectory.ReadFooter ( )
        trajectory.Close ( )
    # . Calculate g(r).
    if hasSymmetry: volume = v / float ( numberFrames )
    fact = 4.0 * math.pi * float ( np1 * np2 * numberFrames ) / ( 3.0 * volume )
    for i in range ( bins ):
        rlower = lower  + float ( i ) * width
        rupper = rlower + width
        distances.append ( rlower + 0.5 * width )
        rdf.append ( float ( histogram[i] ) / ( fact * ( rupper**3 - rlower**3 ) ) )
    # . Output the results.
    if LogFileActive ( log ):
        table = log.GetTable ( columns = [ 20, 20 ] )
        table.Start ( )
        table.Title ( "Radial Distribution Function" )
        table.Heading ( "Distance" )
        table.Heading ( "g(r)"     )
        for ( r, g ) in zip ( distances, rdf ):
            table.Entry ( "{:20.4f}".format ( r ) )
            table.Entry ( "{:20.4f}".format ( g ) )
        table.Stop ( )
    # . Return results.
    return ( distances, rdf )

#===================================================================================================================================
# . Remove rotational and translational motion from a trajectory.
# . This function can also be used for superimposition.
#===================================================================================================================================
# . To remove rotation and translation the weights should be the masses using:
#       Array.FromIterable ( [ atom.mass for atom in system.atoms ] )  
def RemoveRotationTranslation ( inPath, outPath, system, inFormat = None, outFormat = None, reference3 = None, weights = None ):
    """Superimpose the trajectory on a reference structure using the first trajectory frame by default."""
    # . Initialization.
    inTrajectory  = ImportTrajectory ( inPath , system, format = inFormat  )
    outTrajectory = ExportTrajectory ( outPath, system, format = outFormat )
    inTrajectory.ReadHeader   ( )
    outTrajectory.WriteHeader ( )
    # . Save system coordinates.
    coordinates3 = system.coordinates3
    saved3       = Clone ( coordinates3 )
    if reference3 is coordinates3: reference3 = saved3
    # . Loop over trajectory frames.
    results = []
    while inTrajectory.RestoreOwnerData ( ):
        if reference3 is None: reference3 = Clone ( coordinates3 )
        rms0 = coordinates3.RootMeanSquareDeviation ( reference3, weights = weights )
        coordinates3.Superimpose ( reference3, weights = weights )
        rms1 = coordinates3.RootMeanSquareDeviation ( reference3, weights = weights )
        outTrajectory.WriteOwnerData ( )
        results.append ( ( rms0, rms1 ) )
    inTrajectory.ReadFooter   ( )
    outTrajectory.WriteFooter ( )
    inTrajectory.Close  ( )
    outTrajectory.Close ( )
    # . Finish up.
    saved3.CopyTo ( coordinates3 )
    return results

#===================================================================================================================================
# . Calculate the self diffusion function.
#===================================================================================================================================
def SelfDiffusionFunction ( path                  ,
                            system                ,
                            format      = None    ,
                            log         = logFile ,
                            maximumTime = None    ,
                            selection   = None    ,
                            timeStep    = 1.0     ):
    """Calculate the self diffusion function."""
    # . Initialization.
    trajectories = _GetSystemAndTrajectoriesFromArguments ( path, system, format = format )
    if selection is None: selection = Selection.FromIterable ( range ( len ( system.atoms ) ) )
    # . Initialization.
    np     = len ( selection ) # . The number of particles.
    dSelf  = []
    frames = []
    times  = []
    # . Loop over trajectory frames.
    for trajectory in trajectories:
        trajectory.ReadHeader ( )
        while trajectory.RestoreOwnerData ( ):
            frames.append ( system.coordinates3.Prune ( selection ) )
        trajectory.ReadFooter ( )
        trajectory.Close ( )
    # . Get the number of frames and tStop.
    numberFrames = len ( frames )
    if maximumTime is None:
        tStop = numberFrames - 1
    else:
        tStop = int ( maximumTime / timeStep )
        tStop = min ( numberFrames - 1, tStop )
    # . Calculate the function - slow version.
    # . Loop over time differences to calculate.
    for t in range ( tStop + 1 ):
        # . Initialization.
        tmax  = numberFrames - t
        total = 0.0
        # . Loop over allowed increments.
        for t0 in range ( tmax ):
            # . Calculate differences.
            fa = frames[t0+t]
            fb = frames[t0]
            for i in range ( np ):
                total += ( fa[i,0] - fb[i,0] )**2 + ( fa[i,1] - fb[i,1] )**2 + ( fa[i,2] - fb[i,2] )**2
        # . Calculate dSelf.
        dSelf.append ( total / float ( 3 * np * tmax ) )
        times.append ( float ( t ) * timeStep )
    # . Output the results.
    if LogFileActive ( log ):
        table = log.GetTable ( columns = [ 20, 20 ] )
        table.Start ( )
        table.Title ( "Self-Diffusion Function" )
        table.Heading ( "Time"  )
        table.Heading ( "Dself" )
        for ( t, d ) in zip ( times, dSelf ):
            table.Entry ( "{:20.4f}".format ( t ) )
            table.Entry ( "{:20.4f}".format ( d ) )
        table.Stop ( )
    # . Return results.
    return ( times, dSelf )

#===================================================================================================================================
# . Separate a trajectory into a set of coordinate files in a given directory with names "frame<int>.extension".
#===================================================================================================================================
def ToCoordinateFiles ( inPath, outPath, system, inFormat = None, outFormat = None, padding = 0 ):
    """Separate a trajectory into separate coordinate files."""
    if not isinstance ( outFormat, str ): raise TypeError ( "An out format must be specified." )
    if not os.path.exists ( outPath ): os.mkdir ( outPath )
    if padding <= 0: format = "frame{:d}" + outFormat
    else:            format = "frame{:0" + str ( padding ) + "d}"
    format += ( "." + outFormat )
    inTrajectory = ImportTrajectory ( inPath, system, format = inFormat )
    inTrajectory.ReadHeader ( )
    i = 1
    while inTrajectory.RestoreOwnerData ( ):
        ExportSystem ( os.path.join ( outPath, format.format ( i ) ), system )
        i += 1
    inTrajectory.ReadFooter ( )
    inTrajectory.Close ( )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__" :
    pass
