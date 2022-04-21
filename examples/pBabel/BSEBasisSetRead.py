"""Test for reading BSE basis set files in Gaussian format."""

import glob, os, os.path

from Definitions import dataPath, outPath
from pBabel      import BSEGaussianFileReader
from pCore       import logFile, TestScriptExit_Fail, YAMLPickle, YAMLPickleFileExtension
from pScientific import PeriodicTable

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . The destination for results.
_Destination = "gaussianBasisSets"

# . The data path.
_Source = "bseGaussian"

# . The file extension.
_Extension = "gbs"

# . Cartesian basis sets.
_IsCartesian = set ( [ "3-21g", "6-31g_st" ] )

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Initialization.
isOK = True

# . Output setup.
dataPath = os.path.join ( dataPath, _Source )
outPath0 = outPath

# . Loop over the files.
for path in sorted ( glob.glob ( os.path.join ( dataPath, "*.{:s}".format ( _Extension ) ) ) ):
    ( head, tail ) = os.path.split ( path )
    label          = tail.rsplit ( ".", 2 )[0]
    isSpherical    = ( label not in _IsCartesian )

    # . Get the bases.
    bases = BSEGaussianFileReader.PathToGaussianBases ( path, isSpherical = isSpherical )
    logFile.Paragraph ( "Processed bases for basis {:s} = {:d}.".format ( label, len ( bases ) ) )

    # . Check for an appropriate outPath.
    if outPath0 is not None:
        outPath = os.path.join ( outPath0, _Destination )
        if not os.path.exists ( outPath ): os.mkdir ( outPath )
        outPath = os.path.join ( outPath, label.lower ( ) )
        if not os.path.exists ( outPath ): os.mkdir ( outPath )

        # . Save the bases.
        for basis in bases:
            symbol = PeriodicTable.Symbol ( basis.atomicNumber )
            YAMLPickle ( os.path.join ( outPath, symbol + YAMLPickleFileExtension ), basis )

# . Footer.
logFile.Footer ( )
if not isOK: TestScriptExit_Fail ( )
