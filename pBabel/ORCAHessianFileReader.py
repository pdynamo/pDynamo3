"""Classes and functions for reading ORCA hessian files."""

import math, os.path

from  pCore                 import logFile         , \
                                   LogFileActive   , \
                                   TextFileReader
from  pMolecule             import System
from  pScientific           import PeriodicTable
from  pScientific.Arrays    import Array
from  pScientific.Geometry3 import Coordinates3
from .ExportImport          import _Importer

#
# . Known headers (most not treated):
#   $act_atom
#   $act_coord
#   $act_energy
#   $actual_temperature
#   $atoms
#   $dipole_derivatives
#   $end
#   $hessian
#   $ir_spectrum
#   $normal_modes
#   $orca_hessian_file
#   $vibrational_frequencies
#

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
_DefaultSystemLabel = "ORCA Hessian File"

#===================================================================================================================================
# . ORCA hessian file reader class.
#===================================================================================================================================
class ORCAHessianFileReader ( TextFileReader ):
    """ORCAHessianFileReader is the class for ORCA hessian files that are to be read."""

    _classLabel = "ORCA Hessian File Reader"

    def GetLine ( self, signalWarnings = False ):
        """Get a non-empty line removed of comments."""
        try:
            found      = ""
            isNotFound = True
            while isNotFound:
                found = next ( self.file ).strip ( )
                self.linesParsed += 1
                index = found.find ( "#" )
                if index >= 0: found = found[:index].strip ( )
                isNotFound = ( len ( found ) <= 0 )
            return found
        except:
            if signalWarnings: self.Warning ( "Unexpected end-of-file.", True )
            if ( len ( found ) > 0 ): return found
            else:                     raise EOFError

    def Parse ( self, log = logFile ):
        """Parse the data on the file."""
        if not self.isParsed:
            self.data = {}
            if LogFileActive ( log ): self.log = log
            self.Open ( )
            try:
                while True:
                    line = self.GetLine ( signalWarnings = False )
                    if   line == "$atoms"                   : self.ParseAtoms                  ( )
                    elif line == "$dipole_derivatives"      : self.ParseDipoleDerivatives      ( )
                    elif line == "$end"                     : break
                    elif line == "$hessian"                 : self.ParseHessian                ( )
                    elif line == "$ir_spectrum"             : self.ParseIRSpectrum             ( )
                    elif line == "$normal_modes"            : self.ParseNormalModes            ( )
                    elif line == "$vibrational_frequencies" : self.ParseVibrationalFrequencies ( )
            except EOFError:
                pass
            self.WarningStop ( )
            self.Close       ( )
            self.log     = None
            self.isParsed = True

    def ParseAtoms ( self ):
        """Parse an atoms block."""
        n = int ( self.GetLine ( ) )
        if n > 0:
            atomicNumbers = []
            masses        = Array.WithExtent        ( n )
            coordinates3  = Coordinates3.WithExtent ( n )
            isOK = True
            for i in range ( n ):
                tokens = self.GetTokens ( converters = [ PeriodicTable.AtomicNumber, float, float, float, float ] )
                if len ( tokens ) == 5:
                    atomicNumbers.append ( tokens[0] )
                    masses[i]         = tokens[1]
                    coordinates3[i,0] = tokens[2]
                    coordinates3[i,1] = tokens[3]
                    coordinates3[i,2] = tokens[4]
                else:
                    isOK = False
                    break
            if isOK:
                self.data["Atomic Numbers"] = atomicNumbers
                self.data["Coordinates3"  ] = coordinates3
                self.data["Masses"        ] = masses

    def ParseDipoleDerivatives ( self ):
        """Parse a dipole derivatives block."""
        n = int ( self.GetLine ( ) )
        if n > 0:
            data = Array.WithExtents ( n, 3 )
            isOK = True
            for i in range ( n ):
                tokens = self.GetTokens ( converters = [ float, float, float ] )
                if len ( tokens ) == 3:
                    data[i,0] = tokens[0]
                    data[i,1] = tokens[1]
                    data[i,2] = tokens[2]
                else:
                    isOK = False
                    break
            if isOK:
                self.data["Dipole Derivatives"] = data

    def ParseHessian ( self ):
        """Parse a hessian block."""
        n = int ( self.GetLine ( ) )
        if n > 0:
            data = Array.WithExtents ( n, n )
            isOK = True
            p    = 0
            for b in range ( math.ceil ( n / 6 ) ):
                m = min ( ( n - p ), 6 )
                converters = [ None ] + m * [ float ]
                self.GetLine ( )
                for i in range ( n ):
                    tokens = self.GetTokens ( converters = converters )
                    if len ( tokens ) == ( m + 1 ):
                        for q in range ( m ): data[i,p+q] = tokens[q+1]
                    else:
                        isOK = False
                        break
                if not isOK: break
                p += m
            if isOK:
                self.data["Hessian"] = data

    def ParseIRSpectrum ( self ):
        """Parse an IR spectrum block."""
        n = int ( self.GetLine ( ) )
        if n > 0:
            data = Array.WithExtents ( n, 5 )
            isOK = True
            for i in range ( n ):
                tokens = self.GetTokens ( converters = [ float, float, float, float, float ] )
                if len ( tokens ) == 5:
                    data[i,0] = tokens[0]
                    data[i,1] = tokens[1]
                    data[i,2] = tokens[2]
                    data[i,3] = tokens[3]
                    data[i,4] = tokens[4]
                else:
                    isOK = False
                    break
            if isOK:
                self.data["IR Spectrum"] = data

    def ParseNormalModes ( self ):
        """Parse a normal modes block."""
        tokens = self.GetTokens ( converters = [ int, int ] )
        if ( len ( tokens ) == 2 ) and ( tokens[0] > 0 ) and ( tokens[1] > 0 ):
            n1   = tokens[0]
            n2   = tokens[1]
            data = Array.WithExtents ( n1, n2 )
            isOK = True
            p    = 0
            for b in range ( math.ceil ( n2 / 6 ) ):
                m = min ( ( n2 - p ), 6 )
                converters = [ None ] + m * [ float ]
                self.GetLine ( )
                for i in range ( n1 ):
                    tokens = self.GetTokens ( converters = converters )
                    if len ( tokens ) == ( m + 1 ):
                        for q in range ( m ): data[i,p+q] = tokens[q+1]
                    else:
                        isOK = False
                        break
                if not isOK: break
                p += m
            if isOK:
                self.data["Normal Modes"] = data

    def ParseVibrationalFrequencies ( self ):
        """Parse a vibrational frequencies block."""
        n = int ( self.GetLine ( ) )
        if n > 0:
            data = Array.WithExtent ( n )
            isOK = True
            for i in range ( n ):
                tokens = self.GetTokens ( converters = [ None, float ] )
                if len ( tokens ) == 2:
                    data[i] = tokens[1]
                else:
                    isOK = False
                    break
            if isOK:
                self.data["Vibrational Frequencies"] = data

    @classmethod
    def PathToCoordinates3 ( selfClass, path, log = logFile ):
        """Return the coordinates from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse ( log = log )
        return inFile.ToCoordinates3 ( )

    @classmethod
    def PathToSystem ( selfClass, path, log = logFile ):
        """Return the system from a file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse ( log = log )
        return inFile.ToSystem ( )

    def ToCoordinates3 ( self ):
        """Return a coordinates3 object."""
        coordinates3 = None
        if self.isParsed: 
            coordinates3 = self.data.get ( "Coordinates3", None )
        return coordinates3

    def ToSystem ( self ):
        """Return a system."""
        system = None
        if self.isParsed:
            atomicNumbers = self.data.get ( "Atomic Numbers", None )
            if atomicNumbers is not None:
                ( head, tail )      = os.path.split ( self.path )
                system              = System.FromAtoms ( atomicNumbers )
                system.label        = _DefaultSystemLabel + " " + tail
                system.coordinates3 = self.ToCoordinates3    ( )
        return system

#===================================================================================================================================
# . Importer definitions.
#===================================================================================================================================
_Importer.AddHandler ( { Coordinates3 : ORCAHessianFileReader.PathToCoordinates3 ,
                         System       : ORCAHessianFileReader.PathToSystem       } ,
                       [ "hess", "ohess", "HESS", "OHESS" ], "ORCA Hessian", defaultFunction = ORCAHessianFileReader.PathToSystem )

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
