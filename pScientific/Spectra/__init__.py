"""A package for miscellaneous operations on spectra."""

from .ColorManipulation import CIERGBColorSpace                       , \
                               CIEXYZColorSpace                       , \
                               EBURGBColorSpace                       , \
                               EstimateAbsorbances                    , \
                               EstimateRGBColorOfSpectrum             , \
                               HDTVRGBColorSpace                      , \
                               NTSCRGBColorSpace                      , \
                               Rec709RGBColorSpace                    , \
                               RGBColor                               , \
                               RGBColorSpace                          , \
                               SMPTERGBColorSpace                     , \
                               sRGBRGBColorSpace                      , \
                               XYZColor
from .SpectraHandling   import BlackBodySpectrum                      , \
                               ContinuousSpectrum                     , \
                               DataRange                              , \
                               DataSet                                , \
                               GaussianLineShape                      , \
                               LorentzianLineShape                    , \
                               RayleighSpectrum                       , \
                               StickSpectrum
