"""Handle basis set exchange Gaussian basis set files.

This module reads BSE basis set files in Gaussian format.
These files may be downloaded from the basis set exchange at:

           https://www.basissetexchange.org
"""

from pCore                           import logFile          , \
                                            LogFileActive    , \
                                            TextFileReader
from pMolecule.QCModel.GaussianBases import GaussianBasis    , \
                                            ShellLabelDecode
from pScientific                     import PeriodicTable

#===================================================================================================================================
# . Class.
#===================================================================================================================================
class BSEGaussianFileReader ( TextFileReader ):
    """The class for reading BSE basis set files in Gaussian format."""

    _classLabel = "BSE Gaussian File Reader"

    def _ProcessRawData ( self, bases ):
        """Process the input basis lines."""
        # . Process the text data.
        self.bases = []
        for basis in bases:
            # . Title line.
            ( iLine, line ) = basis.pop ( 0 )
            words           = line.split ( )
            atomicNumber    = PeriodicTable.AtomicNumber ( words[0] )
            if atomicNumber <= 0: self.Warning ( "Invalid element, {:s}, for basis on line {:d}.".format ( words[0], iLine ), True )
            # . Shell lines.
            shells = []
            while len ( basis ) > 0:
                # . Header.
                ( iLine, line ) = basis.pop ( 0 )
                words           = line.split ( )
                shellID         = words[0].upper ( )
                try:
                    ( lLow, lHigh ) = ShellLabelDecode ( shellID )
                except:
                    self.Warning ( "Unrecognized shell type, {:s}, for basis on line {:d}.".format ( shellID, iLine ), True )
                nprimitives = int   ( words[1] )
                scale       = float ( words[2].replace ( "D", "e" ) )
                # . Primitive lines.
                primitives = []
                if len ( basis ) >= nprimitives:
                    for p in range ( nprimitives ):
                        primitive       = []
                        ( iLine, line ) = basis.pop ( 0 )
                        words           = line.split ( )
                        if ( len ( words ) != ( len ( shellID ) + 1 ) ): self.Warning ( "Invalid number of primitive data on line {:d}.".format ( iLine ), True )
                        for word in words: primitive.append ( float ( word.replace ( "D", "e" ) ) )
                        primitives.append ( primitive )
                else: self.Warning ( "Invalid number of primitives for basis on line {:d}.".format ( iLine ), True )
                # . Save the data.
                shells.append ( ( lLow, lHigh, nprimitives, scale, primitives ) )
            # . Save the basis.
            self.bases.append ( ( atomicNumber, shells ) )

    def Parse ( self, log = logFile ):
        """Parse the data on the file."""
        if not self.isParsed:
            # . Initialization.
            if LogFileActive ( log ): self.log = log
            # . Open the file.
            self.Open ( )
            # . Read the raw data.
            try:
                basis = []
                bases = []
                for ( iLine, line ) in enumerate ( self.file ):
                    line = line.strip ( )
                    if ( len ( line ) <= 0 ) or line.startswith ( "!" ): continue
                    # . Header line.
                    elif line.startswith ( "BASIS" ):
                        self.basisname = line[7:-1]
                    # . Terminator line.
                    elif line == "****":
                        if len ( basis ) > 0: bases.append ( basis )
                        basis = []
                    # . General line.
                    else:
                        basis.append ( ( iLine, line ) )
            except EOFError:
                pass
            # . Process the data.
            self._ProcessRawData ( bases )
            # . Close the file.
            self.WarningStop ( )
            self.Close ( )
            # . Set the parsed flag and some other options.
            self.log     = None
            self.isParsed = True

    @classmethod
    def PathToGaussianBases ( selfClass, path, isSpherical = True, log = logFile ):
        """Read Gaussian basis sets in Gaussian format from a BSE file."""
        inFile = selfClass.FromPath ( path )
        inFile.Parse ( log = log )
        return inFile.ToGaussianBases ( isSpherical = isSpherical )

    def ToGaussianBases ( self, isSpherical = True ):
        """Return the gaussian bases from the file."""
        if self.isParsed:
            if not hasattr ( self, "gaussianBases" ):
                self.gaussianBases = []
                if self.bases is not None:
                    for ( atomicNumber, shells ) in self.bases:
                        basis  = GaussianBasis.Uninitialized ( atomicNumber, len ( shells ), isSpherical = isSpherical )
                        iShell = 0
                        for ( lLow, lHigh, nprimitives, scale, primitives ) in shells:
                            basis.AllocateShell ( iShell, lHigh, lLow, nprimitives )
                            coefficients = []
                            exponents    = []
                            for i in range ( lHigh - lLow + 1 ): coefficients.append ( [] )
                            for primitive in primitives:
                                exponents.append ( primitive[0] )
                                for i in range ( lHigh - lLow + 1 ): coefficients[i].append ( primitive[i+1] )
                            for ( i, l ) in enumerate ( range ( lLow, lHigh + 1 ) ):
                                basis.SetShellCoefficients ( iShell, l, coefficients[i] )
                            basis.SetShellExponents ( iShell, exponents, scale = scale )
                            iShell += 1
                        self.gaussianBases.append ( basis )
            return self.gaussianBases
        else:
            return None

#===================================================================================================================================
# . Testing.
#===================================================================================================================================
if __name__ == "__main__":
    pass
