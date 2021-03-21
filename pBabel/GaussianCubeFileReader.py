"""Classes and functions for reading Gaussian cube files.

See the preamble to GaussianCubeFileWriter for more information about the format.
"""

import math

from  pCore                 import logFile        , \
                                   LogFileActive  , \
                                   TextFileReader
from  pMolecule             import System
from  pScientific           import Units
from  pScientific.Arrays    import Array
from  pScientific.Geometry3 import Coordinates3   , \
                                   RegularGrid
from .ExportImport          import _Importer

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class GaussianCubeFileReader ( TextFileReader ):
    """Class for reading Gaussian cube files."""

    _classLabel = "Gaussian Cube File Reader"

    def Parse ( self, log = logFile ):
        """Parsing."""
        if not self.isParsed:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
            # . Open the file.
            self.Open ( )
            # . Parse the data.
            try:
                # . Header.
                self.label = self.GetLine ( )
                # . More description - ignored.
                self.GetLine ( )
                # . First definition line.
                ( nAtoms, oX, oY, oZ ) = self.GetFixedFormatTokens ( ( 0, 5, int, 0 ), ( 5, 17, float, 0.0 ), ( 17, 29, float, 0.0 ), ( 29, 41, float, 0.0 ) )
                hasOrbitals = ( nAtoms < 0 )
                if hasOrbitals: nAtoms = abs ( nAtoms )
                # . Grid axis definitions - number of points and increment vector.
                gridData = []
                ngrid    = 1
                for ( i, o ) in enumerate ( ( oX, oY, oZ ) ):
                    items = self.GetFixedFormatTokens ( ( 0, 5, int, 0 ), ( 5, 17, float, 0.0 ), ( 17, 29, float, 0.0 ), ( 29, 41, float, 0.0 ) )
                    n     = items.pop ( 0 )
                    h     = items.pop ( i ) * Units.Length_Bohrs_To_Angstroms
                    o    *= Units.Length_Bohrs_To_Angstroms
                    if ( items[0] != 0.0 ) or ( items[1] != 0.0 ): self.Warning ( "Unable to handle grids whose increment vectors are not aligned along the Cartesian axes.", False )
                    gridData.append ( { "bins" : n, "binSize" : h, "lower" : o - 0.5 * h } )
                    ngrid *= n
                # . Get the grid.
                self.grid = RegularGrid.FromDimensionData ( gridData )
                # . Atom data.
                self.atomicNumbers = []
                self.coordinates3  = Coordinates3.WithExtent ( nAtoms )
                for i in range ( nAtoms ):
                    ( n, q, x, y, z ) = self.GetFixedFormatTokens ( ( 0, 5, int, 0 ), ( 5, 17, float, 0.0 ), ( 17, 29, float, 0.0 ), ( 29, 41, float, 0.0 ), ( 41, 53, float, 0.0 ) )
                    self.atomicNumbers.append ( n )
                    self.coordinates3[i,0] = x
                    self.coordinates3[i,1] = y
                    self.coordinates3[i,2] = z
                self.coordinates3.Scale ( Units.Length_Bohrs_To_Angstroms )
                # . Orbital data.
                if hasOrbitals:
                    indices = self.GetTokens ( )
                    for ( i, index ) in enumerate ( indices ): indices[i] = int ( index )
                    norbitals = indices.pop ( 0 )
                    if norbitals != len ( indices ): self.Warning ( "The number of orbitals does not match the number of orbital indices.", True )
                    self.orbitalIndices = indices
                else:
                    self.orbitalIndices = None
                # . Field data.
                # . Unfortunately each row is written out separately so reading in has to be done by row as well.
                keys = []
                if hasOrbitals:
                    for index in self.orbitalIndices:
                        keys.append ( "Orbital {:d}".format ( index ) )
                else:
                    keys.append ( "Density" )
                self.fieldData = {}
                for key in keys:
                    self.fieldData[key] = Array.WithExtents ( *self.grid.shape )
                nFields = len ( keys )
                if nFields == 1:
                    fieldData = self.fieldData[keys[0]] # []
                else:
                    fieldData = []
                    for key in keys:
                        fieldData.append ( self.fieldData[key] )
                nitems  = gridData[2]["bins"] * nFields
                for ix in range ( gridData[0]["bins"] ):
                    for iy in range ( gridData[1]["bins"] ):
                        newData = self.GetFixedFormatArray ( nitems, 6, 13, converter = float, default = 0.0 )
                        if nFields == 1:
#                            fieldData.extend ( newData )
                            for iz in range ( nitems ):
                                fieldData[ix,iy,iz] = newData[iz]
                        else:
                            for ifield in range ( nFields ):
                                for ( iz, i ) in enumerate ( range ( ifield, nitems, nFields ) ):
                                    fieldData[ifield][ix,iy,iz] = newData[i]
#                            for ifield in range ( nFields ):
#                                fieldData[ifield].extend ( newData[ifield::nFields] )
            except EOFError:
                pass
            # . Close the file.
            self.WarningStop ( )
            self.Close ( )
            # . Set the parsed flag and some other options.
            self.log     = None
            self.isParsed = True

    @classmethod
    def PathToCoordinates3 ( selfClass, path, log = logFile ):
        """Return the coordinates from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse ( log = log )
        return inFile.ToCoordinates3 ( )

    @classmethod
    def PathToSystem ( selfClass, path, log = logFile, volumetricData = None ):
        """Return the system from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse ( log = log )
        if volumetricData is not None: volumetricData.update ( inFile.ToVolumetricData ( ) )
        return inFile.ToSystem ( )

    def ToCoordinates3 ( self ):
        """Return a coordinates3 object."""
        if self.isParsed: return self.coordinates3
        else:             return None

    def ToSystem ( self ):
        """Return a system."""
        if self.isParsed:
            system              = System.FromAtoms ( self.atomicNumbers )
            system.label        = self.label
            system.coordinates3 = self.ToCoordinates3 ( )
            return system
        else:
            return None

    def ToVolumetricData ( self ):
        """Return the volumetric data."""
        # . Temporarily the grid will go here too.
        data = {}
        if self.isParsed:
            data = self.fieldData
            data["Grid"] = self.grid
        return data

#===================================================================================================================================
# . Importer definitions.
#===================================================================================================================================
_Importer.AddHandler ( { Coordinates3 : GaussianCubeFileReader.PathToCoordinates3 ,
                         System       : GaussianCubeFileReader.PathToSystem       } ,
                       [ "cub", "CUB", "cube", "CUBE" ], "Gaussian Cube File", defaultFunction = GaussianCubeFileReader.PathToSystem )

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass
