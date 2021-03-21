"""Marching cubes for surface generation."""

import os, os.path

from math import exp, sqrt, pi

from Definitions               import outPath
from pCore                     import logFile                    , \
                                      TestScriptExit_Fail
from pScientific.Arrays        import Array                      , \
                                      ArrayPrint
from pScientific.Geometry3     import Coordinates3               , \
                                      RegularGrid
from pScientific.Surfaces      import MarchingCubes_Isosurface3D , \
                                      STLFileReader              , \
                                      STLFileWriter

"""
Formulae from: https://winter.group.shef.ac.uk/orbitron/AOs/

For 4f z^3 function have:

n     = 4 (principal quantum number)
rho   = 2*Z*r/n
R4f   = (1/96 sqrt ( 35 )) rho^3 Z^(3/2) exp(-rho/2)
Y4fz3 = sqrt(7/4) z(5z2 - 3r2) / r3 sqrt (1/4 pi)
"""

# . Options.
_Bins      = 50
_BinSize   = 0.2
_IsoValue  = 0.05
_Lower     = - 0.5 * float ( _Bins ) * _BinSize
_N         = 4.0
_Surfaces  = "surfaces"
_Z         = 8.0

# . Factors.
_2ZN       = 2.0 * _Z / _N
_Prefactor = sqrt ( 7.0 * _Z**3 / ( 16.0 * 35.0 * pi ) ) * _2ZN**3 / 96.0

def Function ( x, y, z ):
    """4f orbital."""
    z2 = z * z
    r2 = x * x + y * y + z2
    r  = sqrt ( r2 )
    return ( _Prefactor * z * ( 5.0 * z2 - 3.0 * r2 ) * exp ( - 0.5 * _2ZN * r ) )

# . Header.
logFile.Header ( )
numberFailed = 0

# . Paths.
outPath = os.path.join ( outPath, _Surfaces )
if not os.path.exists ( outPath ): os.mkdir ( outPath )

# . Create the grid and grid coordinates.
axes = [ { "bins" : _Bins, "binSize" : _BinSize, "lower" : _Lower } for d in range ( 3 ) ]
grid = RegularGrid.FromDimensionData ( axes )
xyz  = Coordinates3.FromGrid         ( grid )

# . Create the data to interpolate.
data = Array.WithShape ( grid.shape )
n    = 0
for iX in range ( axes[0]["bins"] ):
    for iY in range ( axes[1]["bins"] ):
            for iZ in range ( axes[2]["bins"] ):
                x = xyz[n,0]             
                y = xyz[n,1]             
                z = xyz[n,2]             
                data[iX,iY,iZ] = Function ( x, y, z )
                n += 1
if False: ArrayPrint ( data, itemFormat = "{:.3f}" )

# . Calculate the surfaces and save them.
for ( isoValue, name ) in ( ( -_IsoValue, "fOrbitalM" ), ( +_IsoValue, "fOrbitalP" ) ):

    # . Generation.
    surface       = MarchingCubes_Isosurface3D ( grid, data, isoValue )
    surface.label = name
    path          = os.path.join ( outPath, "{:s}.stl".format ( name ) )
    STLFileWriter.PathFromPolygonalSurface ( path, surface, name = name )

    # . Reloading.
    reloaded       = STLFileReader.PathToPolygonalSurface ( path, name = name )
    reloaded.label = name + " reloaded"
    surface.Summary  ( )
    reloaded.Summary ( )

# . Footer.
logFile.Footer ( )
if numberFailed > 0: TestScriptExit_Fail ( )
