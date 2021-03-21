"""DFT-D2 QC dispersion model energy."""

#===================================================================================================================================
# . Function.
#===================================================================================================================================
def QCDispersionDFTD2_Energy ( s6, sR, dR, RealArray1D  sqrtC6       not None ,
                                           RealArray1D  r0           not None ,
                                           Coordinates3 coordinates3 not None ,
                                           Coordinates3 gradients3            ):
    """DFT-D2 energy."""
    cdef CReal         energy
    cdef CRealArray2D *cGradients3 = NULL
    if gradients3 is not None: cGradients3 = gradients3.cObject
    energy = CQCDispersionDFTD2_Energy ( s6                   ,
                                         sR                   ,
                                         dR                   ,
                                         sqrtC6.cObject       ,
                                         r0.cObject           ,
                                         coordinates3.cObject ,
                                         cGradients3          )
    return energy
