"""Classes and functions for writing Gaussian cube files.

The following is taken from VMD 1.8:

Write a formatted cubefile very similar to those created by the gaussian program or the cubegen utility.
The format is as follows (last checked against gaussian 98):

Line   Format      Contents
===============================================================
 1     a           Title.
 2     a           Description of the quantity stored in the cubefile.
 3     i5,3f12.6   #atoms, x-,y-,z-coordinates of origin, X0, Y0, Z0
 4     i5,3f12.6   #gridPoints, increment vector, N1, X1, Y1, Z1 - slowest running direction.
 5     i5,3f12.6   #gridPoints, increment vector, N2, X2, Y2, Z2
 6     i5,3f12.6   #gridPoints, increment vector, N3, X3, Y3, Z3 - fastest running direction.

 #atoms lines of atom coordinates:
 ...   i5,4f12.6   Atom number, charge, x-,y-,z-coordinate.

 if #atoms is negative:
 1 line 10i5       Number of orbitals and their numbers.

 rest: 6e13.5      Cube data (with z increment moving fastest,
                   then y and then x) unless there is more than
                   one orbital in which case the orbital number
                   increments fastest. Each row is written out
                   separately so there will be partly empty lines
                   unless N3 is a multiple of 6.

 Other formats that appear to be possible are 6f12.6 and 3d23.16.

All coordinates are in atomic units.

The coordinates of the grid point (I1,I2,I3) are (Fortran convention):

X0+(I1-1)*X1+(I2-1)*X2+(I3-1)*X3
Y0+(I1-1)*Y1+(I2-1)*Y2+(I3-1)*Y3
Z0+(I1-1)*Z1+(I2-1)*Z2+(I3-1)*Z3

"""

import math

from pCore                 import logFile            , \
                                  LogFileActive      , \
                                  SummarizableObject , \
                                  TextFileWriter
from pMolecule             import System
from pMolecule.QCModel     import SpinType
from pScientific           import Units
from pScientific.Arrays    import Array
from pScientific.Geometry3 import Coordinates3       , \
                                  Vector3

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
_DefaultGridSpacing  = 0.2 # . Atomic units.
_DefaultRadiusFactor = 3.0
_DefaultZeroDensity  = 1.0e-10

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class RegularCubicGrid3 ( SummarizableObject ):
    """Defines a regular cubic grid."""

    _attributable = dict ( SummarizableObject._attributable )
    _classLabel   = "Regular 3D Cubic Grid"
    _attributable.update ( { "coordinates3"     : None ,
                             "gridCoordinates3" : None ,
                             "gridSpacing"      : 0.0  ,
                             "nPoints"          : list ,
                             "origin"           : None ,
                             "radii"            : None } )

    def _CheckOptions ( self ):
        """Check options."""
        super ( RegularCubicGrid3, self )._CheckOptions ( )
        ( self.origin, extents ) = self.coordinates3.EnclosingOrthorhombicBox ( radii = self.radii )
        # . Determine the number of points and readjust the origin and extents accordingly.
        for i in range ( len ( extents ) ):
            n = int ( math.ceil ( extents[i] / self.gridSpacing ) )
            l = float ( n ) * self.gridSpacing
            self.origin[i] -= 0.5 * ( l - extents[i] )
            extents[i] = l
            self.nPoints.append ( n )

    @classmethod
    def Construct ( selfClass, gridSpacing, coordinates3, radii ):
        """Constructor.

        |gridSpacing| defines the grid spacing. It is fixed.
        |coordinates3| are the coordinates of the points to be enclosed.
        |radii| are the ranges around each point that must be included in the grid.

        The grid takes the same units as those of coordinates3 and radii.
        """
        return selfClass.WithOptions ( coordinates3 = coordinates3 ,
                                       gridSpacing  = gridSpacing  , 
                                       radii        = radii        )

    def GetGridPointCoordinates ( self ):
        """Return the grid point coordinates."""
        if self.gridCoordinates3 is None:
            coordinates3 = Coordinates3.WithExtent ( self.numberOfPoints )
            n = 0
            for iX in range ( self.nPoints[0] ):
                for iY in range ( self.nPoints[1] ) :
                    for iZ in range ( self.nPoints[2] ) :
                        point = self.Point ( iX, iY, iZ )
                        coordinates3[n,0] = point[0]
                        coordinates3[n,1] = point[1]
                        coordinates3[n,2] = point[2]
                        n += 1
            self.gridCoordinates3 = coordinates3
        return self.gridCoordinates3

    def Point ( self, iX, iY, iZ ):
        """Return the coordinates of a point."""
        point   = Vector3.Null ( )
        point[0] = self.origin[0] + float ( iX ) * self.gridSpacing
        point[1] = self.origin[1] + float ( iY ) * self.gridSpacing
        point[2] = self.origin[2] + float ( iZ ) * self.gridSpacing
        return point

    def SummaryItems ( self ):
        """Summary items."""
        return [ ( "Points"       , "{:d}"  .format (         self.numberOfPoints ) ) ,
                 ( "Grid Spacing" , "{:.2f}".format (         self.gridSpacing    ) ) ,
                 ( "Origin x"     , "{:.2f}".format (         self.origin [0]     ) ) ,
                 ( "Extent x"     , "{:.2f}".format ( float ( self.nPoints[0] ) * self.gridSpacing ) ) ,
                 ( "Origin y"     , "{:.2f}".format (         self.origin [1]     ) ) ,
                 ( "Extent y"     , "{:.2f}".format ( float ( self.nPoints[1] ) * self.gridSpacing ) ) ,
                 ( "Origin z"     , "{:.2f}".format (         self.origin [2]     ) ) ,
                 ( "Extent z"     , "{:.2f}".format ( float ( self.nPoints[2] ) * self.gridSpacing ) ) ,
                 ( "x Points"     , "{:d}"  .format (         self.nPoints[0]     ) ) ,
                 ( "y Points"     , "{:d}"  .format (         self.nPoints[1]     ) ) ,
                 ( "z Points"     , "{:d}"  .format (         self.nPoints[2]     ) ) ]

    @property
    def numberOfPoints ( self ):
        """The number of points in the grid."""
        return ( self.nPoints[0] * self.nPoints[1] * self.nPoints[2] )

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class GaussianCubeFileWriter ( TextFileWriter ):
    """Class for writing Gaussian cube files."""

    _attributable = dict ( TextFileWriter._attributable )
    _attributable.update (  { "atomicNumbers"  : None ,
                              "grid"           : None ,
                              "qcCoordinates3" : None ,
                              "system"         : None } )

    def _CheckOptions ( self ):
        """Constructor."""
        super ( GaussianCubeFileWriter, self )._CheckOptions ( )
        if not isinstance ( self.system, System ): raise TypeError ( "Invalid system." )
        self.atomicNumbers  = self.system.qcState.atomicNumbers
        self.qcCoordinates3 = self.system.scratch.qcCoordinates3AU

    def DefineGrid ( self, gridSpacing = _DefaultGridSpacing, radiusFactor = _DefaultRadiusFactor ):
        """Define the grid."""
        allRadii = Array.FromIterable ( [ atom.vdwRadius for atom in self.system.atoms ] )
        if len ( self.atomicNumbers ) == len ( allRadii ):
            radii = allRadii
        else:
            radii = Array.WithExtent ( len ( self.atomicNumbers ) )
            for ( i, s ) in enumerate ( self.qcAtoms ): radii[i] = allRadii[s]
        radii.Scale ( radiusFactor / Units.Length_Bohrs_To_Angstroms )
        self.grid = RegularCubicGrid3.Construct ( gridSpacing, self.qcCoordinates3, radii )

    @classmethod
    def PathFromSystemDensity ( selfClass, path, system, gridSpacing  = _DefaultGridSpacing  ,
                                                         log          = logFile              ,
                                                         radiusFactor = _DefaultRadiusFactor ,
                                                         spinType     = None                 ):
        """Write a system's electronic density to a Gaussian cube file."""
        outFile = selfClass.FromPath ( path, system = system )
        outFile.DefineGrid      ( gridSpacing, radiusFactor )
        outFile.grid.Summary    ( log = log )
        outFile.WriteHeader     ( )
        outFile.WriteSystemData ( "Density", spinType = spinType )

    @classmethod
    def PathFromSystemOrbitals ( selfClass, path, system, orbitalIndices, gridSpacing  = _DefaultGridSpacing  ,
                                                                          log          = logFile              ,
                                                                          radiusFactor = _DefaultRadiusFactor ,
                                                                          spinType     = None                 ):
        """Write a system's orbitals to a Gaussian cube file."""
        outFile = selfClass.FromPath ( path, system = system )
        outFile.DefineGrid         ( gridSpacing, radiusFactor )
        outFile.grid.Summary       ( log = log )
        outFile.WriteHeader        ( doOrbitals = True )
        outFile.WriteOrbitalHeader ( orbitalIndices )
        outFile.WriteSystemData    ( "Orbitals", orbitalIndices = orbitalIndices, spinType = spinType )

    @classmethod
    def PathFromSystemPotential ( selfClass, path, system, gridSpacing  = _DefaultGridSpacing  ,
                                                           log          = logFile              ,
                                                           radiusFactor = _DefaultRadiusFactor ,
                                                           spinType     = None                 ):
        """Write a system's electrostatic potential to a Gaussian cube file."""
        outFile = selfClass.FromPath ( path, system = system )
        outFile.DefineGrid      ( gridSpacing, radiusFactor )
        outFile.grid.Summary    ( log = log )
        outFile.WriteHeader     ( )
        outFile.WriteSystemData ( "Potential", spinType = spinType )

    def WriteSystemData ( self, dataType, orbitalIndices = None, spinType = None ):
        """Get the system data."""
        # . Grid points.
        gridPoints = self.grid.GetGridPointCoordinates ( )
        # . Generate the data.
        qcModel = self.system.qcModel
        if   dataType == "Density"  : data = qcModel.GridPointDensities  ( self.system, gridPoints ,                 spinType = spinType )
        elif dataType == "Orbitals" : data = qcModel.GridPointOrbitals   ( self.system, gridPoints , orbitalIndices, spinType = spinType )
        elif dataType == "Potential": data = qcModel.GridPointPotentials ( self.system, gridPoints ,                 spinType = spinType )
        else: raise ValueError ( "Unknown grid data type: " + dataType + "." )
        # . Get the multiplicity of the data.
        if dataType == "Orbitals": nDatum = len ( orbitalIndices )
        else:                      nDatum = 1
        # . Write the data - explicit loop required here as each row is written out separately.
        iData = iter ( data )
        n     = 0
        for iX in range ( self.grid.nPoints[0] ):
            for iY in range ( self.grid.nPoints[1] ) :
                for iZ in range ( self.grid.nPoints[2] ) :
                    for i in range ( nDatum ):
                        self.file.write ( "{:13.5e}".format ( next ( iData ) ) )
                        n += 1
                        if ( n > 5 ):
                            self.file.write ( "\n" )
                            n = 0
                if n != 0 :
                    self.file.write ( "\n" )
                    n = 0
        # . Finish up.
        self.Close ( )

    def WriteHeader ( self, doOrbitals = False ):
        """Write the header to the file."""
        self.Open ( )
        if self.system.label is None: self.file.write ( "Gaussian Cube File.\n" )
        else:                         self.file.write ( "Gaussian Cube File for " + self.system.label + ".\n" )
        self.file.write ( "X: outer loop, Y: middle loop, Z: inner loop.\n" )
        n = len ( self.atomicNumbers )
        if doOrbitals: n *= -1
        self.file.write ( "{:5d}{:12.6f}{:12.6f}{:12.6f}\n".format ( n, self.grid.origin[0], self.grid.origin[1], self.grid.origin[2] ) )
        self.file.write ( "{:5d}{:12.6f}{:12.6f}{:12.6f}\n".format ( self.grid.nPoints[0], self.grid.gridSpacing, 0.0, 0.0 ) )
        self.file.write ( "{:5d}{:12.6f}{:12.6f}{:12.6f}\n".format ( self.grid.nPoints[1], 0.0, self.grid.gridSpacing, 0.0 ) )
        self.file.write ( "{:5d}{:12.6f}{:12.6f}{:12.6f}\n".format ( self.grid.nPoints[2], 0.0, 0.0, self.grid.gridSpacing ) )
        for ( i, n ) in enumerate ( self.atomicNumbers ):
            self.file.write ( "{:<5d}{:12.6f}{:12.6f}{:12.6f}{:12.6f}\n".format ( n, float ( n ), self.qcCoordinates3[i,0], self.qcCoordinates3[i,1], self.qcCoordinates3[i,2] ) )

    # . No checking - assumed to be correct.
    def WriteOrbitalHeader ( self, orbitalIndices ):
        """Write the orbital header."""
        self.file.write ( "{:5d}".format ( len ( orbitalIndices ) ) )
        for i in orbitalIndices: self.file.write ( "{:5d}".format ( i ) )
        self.file.write ( "\n" )

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass
