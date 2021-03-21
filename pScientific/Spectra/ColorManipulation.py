"""Color manipulation."""

import math, os, os.path

from   pCore           import AttributableObject       , \
                              Clone                    , \
                              logFile                  , \
                              LogFileActive            , \
                              YAMLMappingFile_ToObject
from  .SpectraHandling import ContinuousSpectrum
from ..Geometry3       import Matrix33                 , \
                              Vector3

#
# . Notes:
#
#   A good site to verify HTML RGB color codes is: "html-color-codes.info/".
#
#   A good site with general information is "http://www.brucelindbloom.com".
#   Useful existing programs are "specrend.c" (John Walker) and "colorpy" (Mark Kness).
#
#   It could be useful to have general names for groups of RGB color codes (rather than specific names).
#

#===================================================================================================================================
# . Colors.
#===================================================================================================================================
class BaseColor ( Vector3 ):
    """Base class for colors represented by three real values."""

    def __init__ ( self ):
        """Constructor."""
        super ( BaseColor, self ).__init__ ( 3 )
        self.Set ( 0.0 )

class RGBColor ( BaseColor ):
    """RGB color."""

    def AddWhite ( self ):
        """Add enough white to make all values positive or zero and scale the resulting values if necessary.."""
        w = - min ( self )
        if w > 0.0: self.Increment ( w )
        if max ( self ) > 1.0: self.Normalize ( )

    def ClampToRange ( self ):
        """Clamp all values to the range [0,1]."""
        for i in range ( 3 ):
            self[i] = max ( 0.0, min ( 1.0, self[i] ) )

    def Normalize ( self ):
        """Normalize so that the largest component has a value of one."""
        large = max ( self )
        if large != 0.0: self.Scale ( 1.0 / large )

    def ToIntegers ( self ):
        """Convert to an integer representation."""
        integers = []
        for i in range ( 3 ):
            integers.append ( min ( 255, max ( 0, int ( round ( 256.0 * self[i] - 0.5 ) ) ) ) )
        return tuple ( integers )

    @property
    def hexString ( self ): 
        """Convert to a hex string representation."""
        return "#{:02X}{:02X}{:02X}".format ( *self.ToIntegers ( ) )

# . Use XYZ + then convert to RGB.
#    @property
#    def luminance ( self ):
#        """Relative luminance."""
#        return ( 0.2126 * self[0] + 0.7152 * self[1] + 0.0722 * self[2] )

class XYZColor ( BaseColor ):
    """XYZ color."""

    @classmethod
    def FromXY ( selfClass, x, y ):
        """Constructor given x and y."""
        self = selfClass ( )
        self[0] = x
        self[1] = y
        self[2] = 1.0 - ( x + y )
        return self

    def Normalize ( self ):
        """Normalization."""
        xyz = sum ( self )
        if xyz != 0.0: self.Scale ( 1.0 / xyz )

    def NormalizeYTo1 ( self ):
        """Normalization so that y = 1."""
        y = self[1]
        if y != 0.0: self.Scale ( 1.0 / y )

    @property
    def luminance ( self ):
        """Relative luminance."""
        return self[1]

#===================================================================================================================================
# . CIE XYZ color space.
#===================================================================================================================================
class CIEXYZColorSpace:
    """CIE XYZ color space."""

    def __init__ ( self, observer = None ):
        """Constructor."""
        self.adaptationMatrices    = {}
        self.illuminants           = {}
        self.illuminantWhitePoints = {}
        self.observer              = observer
        self._SetUp ( )

    def _SetUp ( self ):
        """Load the CIE XYZ data."""
        # . Check the standard observer.
        if ( self.observer == 10 ) or ( self.observer == 1964 ):
            self.observer = 10
            tag = "cie64"
        else:
            self.observer = 2
            tag = "cie31"
        # . Load the files.
        path = os.path.join ( os.getenv ( "PDYNAMO3_PARAMETERS" ), "cieSpectra", "cie31" )
        for label in ( "X", "Y", "Z" ):
            setattr ( self, label.lower ( ), YAMLMappingFile_ToObject ( path + label + ".yaml", ContinuousSpectrum ) )

    def ChromaticAdaptation ( self, xyz, xyzSourceIlluminant, xyzDestinationIlluminant, option = None ):
        """Apply a chromatic adaptation to an xyz color."""
        if xyzDestinationIlluminant is not xyzSourceIlluminant:
            # . Get the transformations to and from the cone response domain.
            ( m, mI ) = self.GetChromaticAdaptationMatrices ( option )
            # . Generate the central part of the transformation.
            d = Vector3.Null ( ) ; xyzDestinationIlluminant.CopyTo ( d ) ; yD = d[1] ; m.ApplyTo ( d )
            s = Vector3.Null ( ) ; xyzSourceIlluminant.CopyTo      ( s ) ; yS = s[1] ; m.ApplyTo ( s )
            d.Divide ( s )  ; d.Scale ( yS / yD )
            # . Apply the transformation.
            m.ApplyTo ( xyz ) ; xyz.Multiply ( d ) ; mI.ApplyTo ( xyz )

    def GetChromaticAdaptationMatrices ( self, option ):
        """Get the chromatic adaptation matrices for a given option."""
        matrices = self.adaptationMatrices.get ( option, None )
        if matrices is None:
            path = os.path.join ( os.getenv ( "PDYNAMO3_PARAMETERS" ), "cieSpectra", option + "ChromaticAdaptation.yaml" )
            if os.path.exists ( path ):
                m  = YAMLMappingFile_ToObject ( path, Matrix33 )
                mI = Matrix33.Null ( )
                m.Invert ( mI )
                m.Transpose ( )
                mI.Transpose ( )
                matrices = ( m, mI )
                self.adaptationMatrices[option] = matrices
            else:
                raise ValueError ( "Invalid chromatic adaptation option: " + option + "." )
        return matrices

    def GetIlluminant ( self, label ):
        """Get a predefined illuminant."""
        illuminant = self.illuminants.get ( label, None )
        if illuminant is None:
            path = os.path.join ( os.getenv ( "PDYNAMO3_PARAMETERS" ), "cieSpectra", "illuminant" + label + ".yaml" )
            if os.path.exists ( path ):
                illuminant = YAMLMappingFile_ToObject ( path, ContinuousSpectrum )
                self.illuminants[label] = illuminant
            else:
                raise ValueError ( "Invalid illuminant: " + label + "." )
        return illuminant

    # . Temporary method until have full illuminant spectra from which can derive their white points.
    def GetIlluminantWhitePoint ( self, label ):
        """Get a predefined illuminant white point."""
        whitePoint = self.illuminantWhitePoints.get ( label, None )
        if whitePoint is None:
            path = os.path.join ( os.getenv ( "PDYNAMO3_PARAMETERS" ), "cieSpectra", "illuminant" + label + "WhitePoint.yaml" )
            if os.path.exists ( path ):
                whitePoint = YAMLMappingFile_ToObject ( path, XYZColor )
                self.illuminantWhitePoints[label] = whitePoint
            else:
                raise ValueError ( "Invalid illuminant white point: " + label + "." )
        return whitePoint

    def SpectrumToXYZ ( self, spectrum, illuminant = None ):
        """Get the XYZ color of a spectrum."""
        others = [ spectrum ]
        if illuminant is not None: others.append ( illuminant )
        xyz    = XYZColor ( )
        for ( i, w ) in enumerate ( ( self.x, self.y, self.z ) ):
            overlap = w.OverlapSpectrum ( others )
            xyz[i]  = overlap.integral
        return xyz

    def WavelengthToXYZ ( self, wavelength, illuminant = None ):
        """Get the XYZ color of a wavelength (delta-function spectrum)."""
        xyz = XYZColor ( )
        for ( i, w ) in enumerate ( ( self.x, self.y, self.z ) ):
            xyz[i] = w.spline.Evaluate ( wavelength )[0]
        if illuminant is not None:
            xyz.Scale ( illuminant.spline.Evaluate ( wavelength )[0] )
        return xyz

#===================================================================================================================================
# . RGB color spaces.
#===================================================================================================================================
class RGBColorSpace ( AttributableObject ):
    """Color systems."""

    _attributable = dict ( AttributableObject._attributable )
    _attributable.update ( { "gamma"         : 2.2  ,
                             "Gamma"         : None ,
                             "gammaOption"   : ""   ,
                             "InverseGamma"  : None ,
                             "label"         : None ,
                             "phosphorBlue"  : None ,
                             "phosphorGreen" : None ,
                             "phosphorRed"   : None ,
                             "phosphorWhite" : None ,
                             "rgbToXYZ"      : None ,
                             "xyxToRGB"      : None } )

    def _CheckOptions ( self ):
        """Check options."""
        self._SetUp      ( )
        self._SetUpGamma ( )

    def _SetUp ( self ):
        """Set up the color system."""
        # . Ensure the defining colors are normalized correctly.
        self.phosphorBlue.Normalize      ( )
        self.phosphorGreen.Normalize     ( )
        self.phosphorRed.Normalize       ( )
        self.phosphorWhite.NormalizeYTo1 ( )
        # . Initialization.
        a = Matrix33.Null ( )
        b = Matrix33.Null ( )
        w = Vector3.Null  ( )
        self.phosphorWhite.CopyTo ( w )
        # . Set up the RGB -> XYZ matrix in a.
        # . Form matrix (r|g|b), invert it and apply to the white point to get the intensities.
        for ( i, v ) in enumerate ( ( self.phosphorRed, self.phosphorGreen, self.phosphorBlue ) ):
            for j in range ( 3 ): a[j,i] = v[j]
        a.Invert  ( b )
        b.ApplyTo ( w )
        # . Scale columns of a by w.
        for i in range ( 3 ):
            v = w[i]
            for j in range ( 3 ): a[j,i] *= v
        # . Invert a to get the XYZ -> RGB matrix.
        a.Invert ( b )
        # . Finish up.
        self.rgbToXYZ = a
        self.xyzToRGB = b

    def _SetUpGamma ( self ):
        """Set up the gamma conversions."""
        option = self.gammaOption.upper ( )
        if   option == "GAMMA":
            self.Gamma        = self.GammaGamma
            self.InverseGamma = self.InverseGammaGamma
        elif option == "L*":
            self.Gamma        = self.GammaLStar
            self.InverseGamma = self.InverseGammaLStar
        elif option == "REC709":
            self.Gamma        = self.GammaRec709
            self.InverseGamma = self.InverseGammaRec709
        elif option == "SRGB":
            self.Gamma        = self.GammaSRGB
            self.InverseGamma = self.InverseGammaSRGB
        else:
            self.Gamma        = None
            self.InverseGamma = None

    # . Gamma companding: linear -> non-linear RGB.
    def GammaCompand ( self, rgb ):
        """Gamma companding."""
        if self.Gamma is not None:
            for i in range ( 3 ): rgb[i] = self.Gamma ( rgb[i] )

    def GammaGamma ( self, x ): return math.pow ( x, 1.0 / self.gamma )

    def GammaLStar ( self, x ):
        """L*."""
        if x <= 216.0 / 24389.0: return ( 24389.0 * x / 2700.0 )
        else:                    return ( 1.16 * math.pow ( x, ( 1.0 / 3.0 ) ) - 0.16 )

    def GammaRec709 ( self, x ):
        """Rec709."""
        if x <= 0.018: return ( 4.5 * x )
        else:          return ( 1.099 * math.pow ( x, 0.45 ) - 0.099 )

    def GammaSRGB ( self, x ):
        """sRGB."""
        if x <= 0.0031308: return ( 12.92 * x )
        else:              return ( 1.055 * math.pow ( x, ( 1.0 / 2.4 ) ) - 0.055 )

    # . Inverse gamma companding: non-linear -> linear RGB.
    def InverseGammaCompand ( self, rgb ):
        """Inverse gamma companding."""
        if self.InverseGamma is not None:
            for i in range ( 3 ): rgb[i] = self.InverseGamma ( rgb[i] )

    def InverseGammaGamma ( self, x ): return math.pow ( x, self.gamma )

    def InverseGammaLStar ( self, x ):
        """L*."""
        if x <= 0.08: return ( 2700.0 * x / 24389.0 )
        else:         return math.pow ( ( x + 0.16 ) / 1.16, 3.0 )

    def InverseGammaRec709 ( self, x ):
        """Rec709."""
        if x <= 0.081: return ( x / 4.5 )
        else:          return math.pow ( ( x + 0.099 ) / 1.099, 1.0 / 0.45 )

    def InverseGammaSRGB ( self, x ):
        """sRGB."""
        if x <= 0.04045: return ( x / 12.92 )
        else:            return math.pow ( ( x + 0.055 ) / 1.055, 2.4 )

    # . Conversion functions.
    def RGBToXYZ ( self, rgb ):
        """Convert an RGB color to an XYZ color."""
        xyz = XYZColor ( )
        rgb.CopyTo ( xyz )
        self.InverseGammaCompand ( xyz )
        self.rgbToXYZ.ApplyTo ( xyz )
        return xyz

    def XYZToRGB ( self, xyz ):
        """Convert an XYZ color to an RGB color."""
        rgb = RGBColor ( )
        xyz.CopyTo ( rgb )
        self.xyzToRGB.ApplyTo ( rgb )
        self.GammaCompand ( rgb )
        return rgb

# . Default systems.
IlluminantC   = XYZColor.FromXY ( 0.3101    , 0.3162    )
IlluminantD65 = XYZColor.FromXY ( 0.3127    , 0.3291    )
IlluminantE   = XYZColor.FromXY ( 1.0 / 3.0 , 1.0 / 3.0 )
CIERGBColorSpace    = RGBColorSpace.WithOptions ( gamma         = "Rec709"                            ,
                                                  label         = "CIE"                               ,
                                                  phosphorRed   = XYZColor.FromXY ( 0.7355 , 0.2645 ) ,
                                                  phosphorGreen = XYZColor.FromXY ( 0.2658 , 0.7243 ) ,
                                                  phosphorBlue  = XYZColor.FromXY ( 0.1669 , 0.0085 ) ,
                                                  phosphorWhite = IlluminantE                         )
EBURGBColorSpace    = RGBColorSpace.WithOptions ( gamma         = "Rec709"                            ,
                                                  label         = "EBU (PAL/SECAM)"                   ,
                                                  phosphorRed   = XYZColor.FromXY ( 0.64   , 0.33   ) ,
                                                  phosphorGreen = XYZColor.FromXY ( 0.29   , 0.60   ) ,
                                                  phosphorBlue  = XYZColor.FromXY ( 0.15   , 0.06   ) ,
                                                  phosphorWhite = IlluminantD65                       )
HDTVRGBColorSpace   = RGBColorSpace.WithOptions ( gamma         = "Rec709"                            ,
                                                  label         = "HDTV"                              ,
                                                  phosphorRed   = XYZColor.FromXY ( 0.670  , 0.330  ) ,
                                                  phosphorGreen = XYZColor.FromXY ( 0.210  , 0.710  ) ,
                                                  phosphorBlue  = XYZColor.FromXY ( 0.150  , 0.060  ) ,
                                                  phosphorWhite = IlluminantD65                       )
NTSCRGBColorSpace   = RGBColorSpace.WithOptions ( gamma         = "Rec709"                            ,
                                                  label         = "NTSC"                              ,
                                                  phosphorRed   = XYZColor.FromXY ( 0.67   , 0.33   ) ,
                                                  phosphorGreen = XYZColor.FromXY ( 0.21   , 0.71   ) ,
                                                  phosphorBlue  = XYZColor.FromXY ( 0.14   , 0.08   ) ,
                                                  phosphorWhite = IlluminantC                         )
Rec709RGBColorSpace = RGBColorSpace.WithOptions ( gamma         = "Rec709"                            ,
                                                  label         = "CIE REC 709"                       ,
                                                  phosphorRed   = XYZColor.FromXY ( 0.64   , 0.33   ) ,
                                                  phosphorGreen = XYZColor.FromXY ( 0.30   , 0.60   ) ,
                                                  phosphorBlue  = XYZColor.FromXY ( 0.15   , 0.06   ) ,
                                                  phosphorWhite = IlluminantD65                       )
SMPTERGBColorSpace  = RGBColorSpace.WithOptions ( gamma         = "Rec709"                            ,
                                                  label         = "SMPTE"                             ,
                                                  phosphorRed   = XYZColor.FromXY ( 0.630  , 0.340  ) ,
                                                  phosphorGreen = XYZColor.FromXY ( 0.310  , 0.595  ) ,
                                                  phosphorBlue  = XYZColor.FromXY ( 0.155  , 0.070  ) ,
                                                  phosphorWhite = IlluminantD65                       )
sRGBRGBColorSpace   = RGBColorSpace.WithOptions ( gamma         = "sRGB"                              ,
                                                  label         = "sRGB"                              ,
                                                  phosphorRed   = XYZColor.FromXY ( 0.64   , 0.33   ) ,
                                                  phosphorGreen = XYZColor.FromXY ( 0.30   , 0.60   ) ,
                                                  phosphorBlue  = XYZColor.FromXY ( 0.15   , 0.06   ) ,
                                                  phosphorWhite = IlluminantD65                       )

#===================================================================================================================================
# . Helper functions.
#===================================================================================================================================
def EstimateAbsorbances ( absorptionFactors, illuminant, spectrum, wavelength ):
    """Estimate the absorbances of a spectrum at a particular wavelength."""
    absorbances = []
    iValue      = max ( illuminant.spline.Evaluate ( wavelength )[0], 0.0 ) / max ( illuminant.ordinates.data ) 
    sValue      = max ( spectrum.spline.Evaluate   ( wavelength )[0], 0.0 ) / max ( spectrum.ordinates.data   )
    for factor in absorptionFactors:
        transmittance = max ( iValue - factor * sValue, 0.0 ) / iValue
        if math.fabs ( transmittance ) > 1.0e-10: absorbances.append ( - math.log10 ( transmittance ) )
        else:                                     absorbances.append ( None )
    return absorbances

def EstimateRGBColorOfSpectrum ( spectrum, absorptionFactors = None    ,
                                           illuminant        = None    ,
                                           log               = logFile ,
                                           rgbColorSpace     = None    ,
                                           xyzColorSpace     = None    ):
    """Estimate the RGB color of a spectrum.

    Transmission spectra are calculated if absorption factors are given.
    """

    # . Suitable defaults.
    doPrinting     = LogFileActive ( log )
    doTransmission = ( absorptionFactors is not None )
    if rgbColorSpace is None: rgbColorSpace = sRGBRGBColorSpace
    if xyzColorSpace is None: xyzColorSpace = CIEXYZColorSpace ( )
    if illuminant    is None: illuminant    = xyzColorSpace.GetIlluminant ( "E" ) # . Equal energy illuminant.

    # . Determine the units of the input spectrum.
    # . This needs to be improved.
    doConvert = ( spectrum.abscissae.units == "cm^-1" )

    # . Create a new spectrum covering the appropriate region.
    localSpectrum = ContinuousSpectrum.WithAbscissae ( Clone ( xyzColorSpace.x.abscissae ) )
    for ( i, x ) in enumerate ( localSpectrum.abscissae.data ):
        if doConvert: x = 10000000.0 / x
        localSpectrum.ordinates.data[i] = spectrum.spline.Evaluate ( x )[0]

    # . Transmission spectra.
    if doTransmission:

        # . Calculation.
        results = []
        for factor in absorptionFactors:
            tSpectrum     = localSpectrum.TransmissionSpectrum ( factor = factor, illuminant = illuminant )
            transmittance = max ( tSpectrum.integral / illuminant.integral, 0.0 )
            if math.fabs ( transmittance ) > 0.0: absorbance = - math.log10 ( transmittance )
            else:                                 absorbance = None
            c = xyzColorSpace.SpectrumToXYZ ( tSpectrum, illuminant = illuminant )
            c.Normalize ( )
            r = rgbColorSpace.XYZToRGB ( c )
            results.append ( ( factor, transmittance, absorbance, r.hexString ) )

        # . Printing.
        if doPrinting and ( len ( absorptionFactors ) > 0 ):
            table = log.GetTable ( columns = [ 20, 20, 20, 20 ] )
            table.Start   ( )
            table.Title   ( "Transmission Spectra RGB Colors" )
            table.Heading ( "Absorption Factor" )
            table.Heading ( "Transmittance (%)" )
            table.Heading ( "Total Absorbance"  )
            table.Heading ( "Color"             )
            for ( factor, transmittance, absorbance, color ) in results:
                table.Entry ( "{:.3f}".format ( factor                 ) )
                table.Entry ( "{:.1f}".format ( transmittance * 100.0  ) )
                if absorbance is None: table.Entry ( "{:.1f}".format ( transmittance * 100.0  ) )
                else:                  table.Entry ( "{:.3f}".format ( absorbance             ) )
                table.Entry ( color )
            table.Stop ( )

    # . Normal spectrum.
    else:
        c       = xyzColorSpace.SpectrumToXYZ ( localSpectrum, illuminant = illuminant )
        c.Normalize ( )
        r       = rgbColorSpace.XYZToRGB ( c )
        results = r.hexString
        if doPrinting: log.Paragraph ( "RGB color code of spectrum = {:s}.".format ( results ) ) 

    # . Finish up.
    return results

#===================================================================================================================================
# . Test.
#===================================================================================================================================
if __name__ == "__main__":
    pass


