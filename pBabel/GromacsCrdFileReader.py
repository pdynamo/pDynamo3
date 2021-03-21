#===================================================================================================================================
#
# . File      : GromacsCrdFileReader.py
#
# . Author    : Guilherme M. Arantes (University of Sao Paulo, Brazil, 2011)
#
# . Based on and to be used with the pDynamo library, copyright CEA, CNRS, Martin J. Field 
#
# . Web       : http://www.pdynamo.org
#
#===================================================================================================================================
"""Read data from a Gromacs .gro coordinate file."""

import sys

from  pCore                 import logFile                  , \
                                   LogFileActive            , \
                                   TextFileReader
from  pScientific.Geometry3 import Coordinates3
from  pScientific.Symmetry  import SymmetryParameters
from .ExportImport          import _Importer

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class GromacsCrdFileReader ( TextFileReader ):
    """Class for reading Gromacs .gro coordinate files."""

    _classLabel = "Gromacs Crd File Reader"

    def GetAtomLineFormat ( self ):
        """Get the format of the atom lines.

        This is:
        RESNO RES   TYPE ATOMNO   X     Y     Z  
        I5    A4 2X A4   I5       F8.3  F8.3  F8.3 
        """
        a = 4 ; f = 8 ; i = 5 ; x = 2
        p = 0
        format =      [ ( p, p+i, int  , 0   ) ] ; p +=  i    # . Res. number.
        format.append ( ( p, p+a, None , ""  ) ) ; p += (a+x) # . Res. name.
        format.append ( ( p, p+a, None , ""  ) ) ; p += (a+i) # . Atom name.
        format.append ( ( p, p+f, float, 0.0 ) ) ; p +=  f    # . X.
        format.append ( ( p, p+f, float, 0.0 ) ) ; p +=  f    # . Y.
        format.append ( ( p, p+f, float, 0.0 ) )              # . Z.
        return format

    def Parse ( self, log = logFile ):
        """Parsing."""
        if not self.isParsed:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
            # . Get the atom line format.
            atomlineformat = self.GetAtomLineFormat ( )
            # . Open the file.
            self.Open ( )
            # . Parse all the lines.
            try:
                # . Keyword line.
                self.title = self.GetLine ( )
                # . Number of atoms.
                items    = self.GetTokens ( converters = ( int, ) )
                natoms   = items[0]
                # . The coordinate data.
                self.xyz = Coordinates3.WithExtent ( natoms )
                for n in range ( natoms ):
                    tokens = self.GetFixedFormatTokens ( *atomlineformat )
                    for i in range ( 3 ): self.xyz[n,i] = float ( tokens[i+3]*10.0 )
                # . Symmetry data.
                self.symmetryItems = self.GetTokens ( )
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
        inFile.Parse   ( log = log )
        inFile.Summary ( log = log )
        return inFile.ToCoordinates3 ( )

    @classmethod
    def PathToSymmetryParameters ( selfClass, path, log = logFile ):
        """Return the symmetry paramters from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse   ( log = log )
        inFile.Summary ( log = log )
        return inFile.ToSymmetryParameters ( )

    def SummaryItems ( self ):
        """Summary items."""
        items = []
        if self.isParsed:
            items.append ( ( "Atom Lines", "{:d}".format ( self.xyz.rows ) ) )
        return items

    def ToCoordinates3 ( self ):
        """Return the coordinates."""
        if self.isParsed: return self.xyz
        else:            return None

    def ToSymmetryParameters ( self ):
        """Return the symmetry parameters."""
        # . Will assign only cubic, orthorhombic, dodecahedron or octahedron boxes. 
        # . Will fail for general triclinic, but this could be hard-coded once sizes/angles are known.
        if self.isParsed: 
            # . triclinic box
            items = self.symmetryItems
            if len ( items ) == 9: 
                specialTriclinic = items[8] == items[7] or "-" + items[7] == items[5]
                if not specialTriclinic: self.Warning ( "Invalid general triclinic box symmetry.", True )
                # . Dodecahedron 
                if   items[0] == items[1]:
                    alpha = 60.0
                    beta  = 60.0
                    gamma = 90.0
                # . Octahedron
                else: 
                    alpha = 70.53
                    beta  = 109.47
                    gamma = 70.53
                a = float ( items[0] ) * 10.0
                b = a
                c = a
            # . Cubic or orthorhombic box
            elif len ( items ) == 3  :
                alpha = 90.0
                beta  = 90.0
                gamma = 90.0
                items = [ float ( items[i] ) * 10.0 for i in range(3) ]
                a = items[0]
                b = items[1]
                c = items[2]
            else: self.Warning ( "Invalid or unrecognized box symmetry.", True )
            self.symmetryParameters = SymmetryParameters ( )
            self.symmetryParameters.SetCrystalParameters (  a = a, b = b, c = c, alpha = alpha, beta = beta, gamma = gamma )
            return self.symmetryParameters
        else:            return None

#===================================================================================================================================
# . Importer definitions.
#===================================================================================================================================
_Importer.AddHandler ( { Coordinates3 : GromacsCrdFileReader.PathToCoordinates3 } , [ "gro", "GRO" ], "Gromacs Coordinates" )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass

