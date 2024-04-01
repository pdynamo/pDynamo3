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

from .ExportImport import _Exporter
from  pCore        import  TextFileWriter
from  pMolecule    import  System

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class GaussianCubeFileWriter ( TextFileWriter ):
    """Class for writing Gaussian cube files."""

    _attributable = dict ( TextFileWriter._attributable )
    _attributable.update (  { "atomicNumbers"  : None ,
                              "qcCoordinates3" : None ,
                              "system"         : None ,
                              "title"          : None } )

    def _CheckOptions ( self ):
        """Constructor."""
        super ( GaussianCubeFileWriter, self )._CheckOptions ( )
        try:
            self.atomicNumbers  = self.system.qcState.atomicNumbers
            self.qcCoordinates3 = self.system.scratch.qcCoordinates3AU
        except: raise TypeError ( "Invalid QC system." )
        if self.title is None:
            if self.system.label is None:
                self.title = "Gaussian Cube File."
            else:
                self.title = "Gaussian Cube File for \"{:s}\".".format ( self.system.label )

    @classmethod
    def PathFromSystem ( selfClass, path, system, grid           = None ,
                                                  gridValues     = None ,
                                                  orbitalIndices = None ,
                                                  title          = None ):
        """Write grid data to a Gaussian cube file."""
        if ( grid is None ) or ( gridValues is None ):
            raise ValueError ( "Missing grid definitions and/or grid data values." )
        self = selfClass.FromPath ( path, system = system, title = title )
        if orbitalIndices is None:
            nDatum = 1
            self.WriteHeader ( grid )
        else:
            nDatum = len ( orbitalIndices )
            self.WriteHeader ( grid, doOrbitals = True )
            self.WriteOrbitalHeader ( orbitalIndices )
        self.WriteData ( grid, nDatum, gridValues )

    def WriteData ( self, grid, nDatum, data ):
        """Write the data."""
        ( nX, nY, nZ ) = grid.shape
        iData = iter ( data )
        n     = 0
        for iX in range ( nX ):
            for iY in range ( nY ) :
                for iZ in range ( nZ ) :
                    for i in range ( nDatum ):
                        self.file.write ( "{:13.5e}".format ( next ( iData ) ) )
                        n += 1
                        if ( n > 5 ):
                            self.file.write ( "\n" )
                            n = 0
                if n != 0 :
                    self.file.write ( "\n" )
                    n = 0
        self.Close ( )

    def WriteHeader ( self, grid, doOrbitals = False ):
        """Write the header to the file."""
        self.Open ( )
        self.file.write ( "{:s}\n".format ( self.title ) )
        self.file.write ( "X: outer loop, Y: middle loop, Z: inner loop.\n" )
        n = len ( self.atomicNumbers )
        if doOrbitals: n *= -1
        ( oX, oY, oZ ) = grid.origin
        ( nX, nY, nZ ) = grid.shape
        ( hX, hY, hZ ) = grid.spacing
        self.file.write ( "{:5d}{:12.6f}{:12.6f}{:12.6f}\n".format ( n, oX, oY, oZ ) )
        self.file.write ( "{:5d}{:12.6f}{:12.6f}{:12.6f}\n".format ( nX, hX, 0.0, 0.0 ) )
        self.file.write ( "{:5d}{:12.6f}{:12.6f}{:12.6f}\n".format ( nY, 0.0, hY, 0.0 ) )
        self.file.write ( "{:5d}{:12.6f}{:12.6f}{:12.6f}\n".format ( nZ, 0.0, 0.0, hZ ) )
        for ( i, n ) in enumerate ( self.atomicNumbers ):
            self.file.write ( "{:<5d}{:12.6f}{:12.6f}{:12.6f}{:12.6f}\n".format ( n, float ( n ), self.qcCoordinates3[i,0] ,
                                                                                                  self.qcCoordinates3[i,1] , 
                                                                                                  self.qcCoordinates3[i,2] ) )

    # . No checking - assumed to be correct.
    def WriteOrbitalHeader ( self, orbitalIndices ):
        """Write the orbital header."""
        self.file.write ( "{:5d}".format ( len ( orbitalIndices ) ) )
        for i in orbitalIndices: self.file.write ( "{:5d}".format ( i ) )
        self.file.write ( "\n" )

#===================================================================================================================================
# . Exporter definitions.
#===================================================================================================================================
_Exporter.AddHandler ( { System : GaussianCubeFileWriter.PathFromSystem } , [ "cub", "CUB", "cube", "CUBE" ], "Gaussian Cube File" )

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass
