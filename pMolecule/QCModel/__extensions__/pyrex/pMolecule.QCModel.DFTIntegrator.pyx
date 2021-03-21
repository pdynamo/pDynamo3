"""Fock construction functions."""

from  pCore              import Clone
from  pScientific.Arrays import Array        , \
                                StorageType
from .QCModelError       import QCModelError

# . This function requires alpha/beta representation of densities and Fock matrices.
# . CDFTIntegrator_Integrate does not initialize Fock matrices on entry.
#
# . Formulae:
#
#   Pa = ( Pt + Ps ) / 2
#   Pb = ( Pt - Ps ) / 2
#
#   Ft = ( Fa + Fb ) / 2
#   Fs = ( Fa - Fb ) / 2
#

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def DFTIntegrator_Integrate ( DFTFunctionalModel     functionalModel not None ,
                              DFTGrid                grid            not None ,
                              GaussianBasisContainer gaussianBases   not None ,
                              Coordinates3           qcCoordinates   not None ,
                              SymmetricMatrix        dTotal          not None ,
                              SymmetricMatrix        dSpin                    ,
                                                     inCore                   ,
                                                     isSpinUnrestricted       ,
                              SymmetricMatrix        fTotal                   ,
                              SymmetricMatrix        fSpin                    ,
                              Coordinates3           gradients3               ):
    """DFT integrator."""
    cdef CBoolean          cInCore
    cdef CBoolean          cIsSpinUnrestricted
    cdef CReal             eQuad
    cdef CReal             rhoQuad
    cdef CRealArray2D     *cGradients3 = NULL
    cdef CStatus           cStatus     = CStatus_OK
    cdef CSymmetricMatrix *cDBeta      = NULL
    cdef CSymmetricMatrix *cFAlpha     = NULL
    cdef CSymmetricMatrix *cFBeta      = NULL
    cdef SymmetricMatrix   dAlpha
    cdef SymmetricMatrix   dBeta
    cdef SymmetricMatrix   fAlpha
    cdef SymmetricMatrix   fBeta
    ( dAlpha, dBeta ) = _MakeDensities ( dTotal, dSpin )
    ( fAlpha, fBeta ) = _MakeFock0     ( fTotal, fSpin )
    if dBeta      is not None: cDBeta      = dBeta.cObject
    if fAlpha     is not None: cFAlpha     = fAlpha.cObject
    if fBeta      is not None: cFBeta      = fBeta.cObject
    if gradients3 is not None: cGradients3 = gradients3.cObject
    if inCore:             cInCore             = CTrue
    else:                  cInCore             = CFalse
    if isSpinUnrestricted: cIsSpinUnrestricted = CTrue
    else:                  cIsSpinUnrestricted = CFalse
    CDFTIntegrator_Integrate ( functionalModel.cObject ,
                               grid.cObject            ,
                               gaussianBases.cObject   ,
                               qcCoordinates.cObject   ,
                               dAlpha.cObject          ,
                               cDBeta                  ,
                               cInCore                 ,
                               cIsSpinUnrestricted     ,
                               &eQuad                  ,
                               &rhoQuad                ,
                               cFAlpha                 ,
                               cFBeta                  ,
                               cGradients3             ,
                               &cStatus                )
    if cStatus != CStatus_OK: raise QCModelError ( "Error in DFT integrator." )
    _MakeFock1 ( fTotal, fSpin, fAlpha, fBeta )
    return ( eQuad, rhoQuad )

#===================================================================================================================================
# . Utility functions.
#===================================================================================================================================
# . It is best to create new alpha and beta densities and Fock matrices rather
# . than make them in place as future multiprocessing may work on the same
# . densities at the same time.
def _MakeDensities ( dTotal, dSpin ):
    """Make alpha and beta densities."""
    if dSpin is None:
        return ( dTotal, None )
    else:
        dAlpha = Clone ( dTotal ) ; dAlpha.Add ( dSpin                ) ; dAlpha.Scale (  0.5 )
        dBeta  = Clone ( dSpin  ) ; dBeta.Add  ( dTotal, scale = -1.0 ) ; dBeta.Scale  ( -0.5 )
        return ( dAlpha, dBeta )

def _MakeFock0 ( fTotal, fSpin ):
    """Create alpha and beta Fock matrices."""
    if fSpin is None:
        return ( fTotal, None )
    else:
        n = fTotal.extent
        fAlpha = Array.WithExtent  ( n, storageType = StorageType.Symmetric ) ; fAlpha.Set ( 0.0 )
        fBeta  = Array.WithExtent  ( n, storageType = StorageType.Symmetric ) ; fBeta.Set  ( 0.0 )
        return ( fAlpha, fBeta )

def _MakeFock1 ( fTotal, fSpin, fAlpha, fBeta ):
    """Increment total and spin Fock matrices."""
    if fSpin is not None:
        fTotal.Add ( fAlpha, scale =  0.5 )
        fTotal.Add ( fBeta , scale =  0.5 )
        fSpin.Add  ( fAlpha, scale =  0.5 )
        fSpin.Add  ( fBeta , scale = -0.5 )
