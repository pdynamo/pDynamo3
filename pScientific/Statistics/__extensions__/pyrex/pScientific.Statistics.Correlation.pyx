"""Functions for handling correlations."""

from pCore import logFile, LogFileActive

#
# . To do:
#
#   Add FFT option (standalone FFT). Require complex type. Compare Direct vs. FFT.
#   Block algorithms for large data sets.
#   CFs of dot products of vectors, etc.
#   Einstein relations (direct and FFT).
#   Integration of CF.
#   Mean and normalization options.
#   Peak location.
#   Smoothing.
#   Spectra.
#

#===================================================================================================================================
# . Parameters.
#===================================================================================================================================

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def Correlation_AutoCorrelation ( RealArray1D x, CInteger correlationLength = -1, normalize = True, removeMean = True, useFFT = False ):
    """Calculate an auto-correlation function."""
    cdef RealArray1D correlation
    cdef CBoolean    cNormalize
    cdef CBoolean    cRemoveMean
    cdef CBoolean    cUseFFT
    # . Get all options.
    if correlationLength <= 0: correlationLength = len ( x ) - 1
    if normalize:  cNormalize  = CTrue
    else:          cNormalize  = CFalse
    if removeMean: cRemoveMean = CTrue
    else:          cRemoveMean = CFalse
    if useFFT:     cUseFFT     = CTrue
    else:          cUseFFT     = CFalse
    # . Calculate the function.
    correlation = RealArray1D.WithExtent ( correlationLength )
    Correlation_MakeSimple ( x.cObject, NULL, cUseFFT, cNormalize, cRemoveMean, correlationLength, NULL, correlation.cObject, NULL )
    # . Finish up.
    return correlation

def Correlation_CrossCorrelation ( RealArray1D x, RealArray1D y, CInteger correlationLength = -1, normalize = True, removeMean = True, useFFT = False ):
    """Calculate a cross-correlation function."""
    cdef RealArray1D correlation
    cdef CBoolean    cNormalize
    cdef CBoolean    cRemoveMean
    cdef CBoolean    cUseFFT
    # . Get all options.
    if correlationLength <= 0: correlationLength = min ( len ( x ) - 1, len ( y ) - 1 )
    if normalize:  cNormalize  = CTrue
    else:          cNormalize  = CFalse
    if removeMean: cRemoveMean = CTrue
    else:          cRemoveMean = CFalse
    if useFFT:     cUseFFT     = CTrue
    else:          cUseFFT     = CFalse
    # . Calculate the function.
    correlation = RealArray1D.WithExtent ( correlationLength )
    Correlation_MakeSimple ( x.cObject, y.cObject, cUseFFT, cNormalize, cRemoveMean, correlationLength, NULL, correlation.cObject, NULL )
    # . Finish up.
    return correlation

def Correlation_DotProductAutoCorrelation ( RealArray2D x, CInteger correlationLength = -1, normalize = True, removeMean = True, useFFT = False ):
    """Calculate a dot-product auto-correlation function."""
    cdef RealArray1D correlation
    cdef CBoolean    cNormalize
    cdef CBoolean    cRemoveMean
    cdef CBoolean    cUseFFT
    # . Get all options.
    if correlationLength <= 0: correlationLength = x.rows - 1
    if normalize:  cNormalize  = CTrue
    else:          cNormalize  = CFalse
    if removeMean: cRemoveMean = CTrue
    else:          cRemoveMean = CFalse
    if useFFT:     cUseFFT     = CTrue
    else:          cUseFFT     = CFalse
    # . Calculate the function.
    correlation = RealArray1D.WithExtent ( correlationLength )
    Correlation_MakeDotProduct ( x.cObject, NULL, cUseFFT, cNormalize, cRemoveMean, correlationLength, NULL, correlation.cObject, NULL )
    # . Finish up.
    return correlation

def Correlation_DotProductCrossCorrelation ( RealArray2D x, RealArray2D y, CInteger correlationLength = -1, normalize = True, removeMean = True, useFFT = False ):
    """Calculate a dot-product cross-correlation function."""
    cdef RealArray1D correlation
    cdef CBoolean    cNormalize
    cdef CBoolean    cRemoveMean
    cdef CBoolean    cUseFFT
    # . Get all options.
    if correlationLength <= 0: correlationLength = min ( x.rows - 1, y.rows - 1 )
    if normalize:  cNormalize  = CTrue
    else:          cNormalize  = CFalse
    if removeMean: cRemoveMean = CTrue
    else:          cRemoveMean = CFalse
    if useFFT:     cUseFFT     = CTrue
    else:          cUseFFT     = CFalse
    # . Calculate the function.
    correlation = RealArray1D.WithExtent ( correlationLength )
    Correlation_MakeDotProduct ( x.cObject, y.cObject, cUseFFT, cNormalize, cRemoveMean, correlationLength, NULL, correlation.cObject, NULL )
    # . Finish up.
    return correlation
