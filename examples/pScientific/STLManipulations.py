"""STL surface checking."""

import glob, os, os.path

from Definitions           import dataPath            , \
                                  outPath
from pBabel                import ExportSystem
from pCore                 import logFile             , \
                                  TestScriptExit_Fail
from pMolecule             import System
from pScientific.Surfaces  import STLFileReader       , \
                                  STLFileWriter       , \
                                  STLSurface

# . Header.
logFile.Header ( )

# . Options.
_InitialXYZ      = False
_Randomize       = True
_ReorientedFiles = True
_Source          = "stl"
_Step            = 0.1
_STL             = "stl"

# . Paths.
inPaths = glob.glob ( os.path.join ( dataPath, _Source, "*.stl" ) )
outPath = os.path.join ( outPath, _STL )
if not os.path.exists ( outPath ): os.mkdir ( outPath )

# . Loop over paths.
nOK    = 0
nTotal = 0
for inPath in sorted ( inPaths ):

    # . Get the basic surface.
    surface = STLFileReader.PathToPolygonalSurface ( inPath )
    surface.MakePolygonNormals ( )
    surface.Summary ( )

    # . STL surface manipulations.
    stl  = STLSurface.FromPolygonalSurface ( surface )
    if _Randomize: stl.RandomizePolygonVertices ( )
    isOK = stl.OrientPolygons ( )
    if isOK: nOK += 1
    nTotal += 1
    stl.Summary ( )
    if _ReorientedFiles:
        ( head, tail ) = os.path.split ( inPath )
        path = os.path.join ( outPath, tail[0:-4] )
        STLFileWriter.PathFromPolygonalSurface ( path + ".stl", surface, name = surface.label )
        system              = System.FromAtoms ( surface.polygons.rows * [ 2 ] )
        system.coordinates3 = stl.MakePolygonCentroids ( normalStep = _Step )
        system.label        = surface.label
        ExportSystem ( path + ".xyz", system )

# . Finish up.
logFile.Paragraph ( "Number of successes = {:d} / {:d}.".format ( nOK, nTotal ) )
logFile.Footer ( )
if nOK < nTotal: TestScriptExit_Fail ( )
