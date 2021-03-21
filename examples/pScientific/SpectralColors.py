"""Tests of various color manipulations."""

import glob, math, os, os.path

from pCore               import Align                    , \
                                Clone                    , \
                                logFile                  , \
                                LogFileActive            , \
                                YAMLMappingFile_ToObject
from pScientific.Spectra import BlackBodySpectrum        , \
                                CIEXYZColorSpace         , \
                                ContinuousSpectrum       , \
                                RayleighSpectrum         , \
                                sRGBRGBColorSpace        , \
                                XYZColor

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================
# . Luminance.
_Luminance = None

# . Reference colors with no luminance.
_ReferenceColors = { "Ansi Status A Blue"  : "#1400E0" ,
                     "Ansi Status A Green" : "#00FF00" ,
                     "Ansi Status A Red"   : "#FF0000" ,
                     "Ansi Status E Blue"  : "#0D00DF" ,
                     "Ansi Status E Green" : "#00FF00" ,
                     "Ansi Status E Red"   : "#FF0800" ,
                     "Ansi Status M Blue"  : "#0800DF" ,
                     "Ansi Status M Green" : "#00FF00" ,
                     "Ansi Status M Red"   : "#FF0000" ,
                     "Ansi Status T Blue"  : "#0007D5" ,
                     "Ansi Status T Green" : "#00FF00" ,
                     "Ansi Status T Red"   : "#FF0800" ,
                     "Black"               : "#000000" ,
                     "Black Body (5778K)"  : "#5F534E" ,
                     "Illuminant A"        : "#C05618" ,
                     "Illuminant D65"      : "#545454" ,
                     "Illuminant E"        : "#66504D" ,
                     "Rayleigh (D65)"      : "#1C3B8A" ,
                     "White"               : "#66504D" }

# . Tolerance for differences.
_Tolerance = 1.0e-03

#===================================================================================================================================
# . Script.
#===================================================================================================================================
# . Header.
logFile.Header ( )

# . Initialization.
colors = []
isOK   = True

# . Define the color space.
cieXYZ = CIEXYZColorSpace ( )

# . Check the white points of some predefined illuminants.
for tag in ( "A", "D65", "E" ):
    s = cieXYZ.GetIlluminant           ( tag )
    w = cieXYZ.GetIlluminantWhitePoint ( tag )
    c = cieXYZ.SpectrumToXYZ ( s )
    c.NormalizeYTo1 ( )
    c.Add ( w, scale = -1.0 )
    isOK = isOK and ( c.AbsoluteMaximum ( ) <= _Tolerance )

# . Get the RGB colors of various spectra.
# . Basic colors.
for ( tag, v ) in ( ( "Black", 0.0 ), ( "White", 1.0 ) ):
    c = XYZColor ( )
    c.Set ( v )
    c.Normalize                    (   )
    r = sRGBRGBColorSpace.XYZToRGB ( c )
    colors.append ( ( tag, r.hexString ) )

# . Illuminants.
for tag in ( "A", "D65", "E" ):
    s = cieXYZ.GetIlluminant ( tag )
    c = cieXYZ.SpectrumToXYZ       ( s )
    c.Normalize                    (   )
    r = sRGBRGBColorSpace.XYZToRGB ( c )
    if _Luminance is not None: r.Scale ( _Luminance / r.luminance )
    colors.append ( ( "Illuminant " + tag, r.hexString ) )

# . Special spectra.
d65 = cieXYZ.GetIlluminant ( "D65" )
s1  = BlackBodySpectrum.WithAbscissae ( Clone ( d65.abscissae ) )
s2  = RayleighSpectrum.FromIlluminant ( d65 )
for ( tag, s ) in ( ( "Black Body ({:.0f}K)".format ( s1.temperature ), s1 ), ( "Rayleigh (D65)", s2 ) ):
    c = cieXYZ.SpectrumToXYZ       ( s )
    c.Normalize                    (   )
    r = sRGBRGBColorSpace.XYZToRGB ( c )
    if _Luminance is not None: r.Scale ( _Luminance / r.luminance )
    colors.append ( ( tag, r.hexString ) )

# . ANSI colors.
paths = glob.glob ( os.path.join ( os.getenv ( "PDYNAMO3_PARAMETERS" ), "cieSpectra", "ansiStatus*.yaml" ) )
paths.sort ( )
for path in paths:
    ( head, tail ) = os.path.split ( path )
    tag = "Ansi Status " + tail[10:11] + " " + tail[11:-5]
    s = YAMLMappingFile_ToObject ( path, ContinuousSpectrum )
    c = cieXYZ.SpectrumToXYZ       ( s )
    c.Normalize                    (   )
    r = sRGBRGBColorSpace.XYZToRGB ( c )
    if _Luminance is not None: r.Scale ( _Luminance / r.luminance )
    colors.append ( ( tag, r.hexString ) )

# . Checks.
if _Luminance is None:
    for ( tag, color ) in colors:
        reference = _ReferenceColors.get ( tag, None )
        if reference is not None:
            isOK = isOK and ( color == reference )

# . Print out.
colors.sort ( )
table = logFile.GetTable ( columns = [ 20, 20 ] )
table.Start   ( )
table.Title   ( "Spectral RGB Colors" )
table.Heading ( "Spectrum" )
table.Heading ( "Color"    )
for ( tag, color ) in colors:
    table.Entry ( tag, align = Align.Left )
    table.Entry ( color )
table.Stop ( )

# . Footer.
logFile.Footer ( )
if ( not isOK ): TestScriptExit_Fail ( )
