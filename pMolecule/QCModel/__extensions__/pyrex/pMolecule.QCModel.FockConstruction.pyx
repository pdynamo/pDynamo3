"""Fock construction functions."""

from .QCModelError import QCModelError

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
# . dTotal is temporarily modified. This needs to be changed.
# . fTotal is not initialized.
def FockConstruction_MakeCoefficientsFromFitIntegrals ( SymmetricMatrix dTotal          not None ,
                                                        BlockStorage    fitIntegrals    not None ,
                                                        SymmetricMatrix fitMatrix       not None ,
                                                        CReal           totalCharge              ,
                                                        RealArray1D     fitCoefficients not None ,
                                                        RealArray1D     fitVector       not None ):
    """Make fit coefficients from the fit integrals."""
    cdef CStatus  cStatus = CStatus_OK
    Fock_MakeCoefficientsFromFitIntegrals ( fitMatrix.cObject       ,
                                            fitIntegrals.cObject    ,
                                            dTotal.cObject          ,
                                            totalCharge             ,
                                            fitCoefficients.cObject ,
                                            fitVector.cObject       ,
                                            &cStatus                )
    if cStatus != CStatus_OK: raise QCModelError ( "Error making coefficients from fit integrals." )

def FockConstruction_MakeFockFromFitIntegrals ( BlockStorage    fitIntegrals not None ,
                                                RealArray1D     fitVector    not None ,
                                                SymmetricMatrix fTotal       not None ):
    """Make Fock matrix from the fit integrals and a fit vector."""
    cdef CStatus  cStatus = CStatus_OK
    Fock_MakeFockFromFitIntegrals ( fitIntegrals.cObject ,
                                    fitVector.cObject    ,
                                    fTotal.cObject       ,
                                    &cStatus             )
    if cStatus != CStatus_OK: raise QCModelError ( "Error constructing Fock matrix from fit integrals." )

def FockConstruction_MakeFromFitIntegralsCoulomb ( SymmetricMatrix dTotal          not None ,
                                                   BlockStorage    fitIntegrals    not None ,
                                                   SymmetricMatrix fitMatrix       not None ,
                                                   CReal           totalCharge              ,
                                                   RealArray1D     fitCoefficients not None ,
                                                   SymmetricMatrix fTotal          not None ):
    """Coulomb Fock matrix from Coulomb fit integrals."""
    cdef CReal   eFit
    cdef CStatus cStatus = CStatus_OK
    eFit = Fock_MakeFromFitIntegralsCoulomb ( fitIntegrals.cObject    ,
                                              fitMatrix.cObject       ,
                                              totalCharge             ,
                                              fitCoefficients.cObject ,
                                              dTotal.cObject          ,
                                              fTotal.cObject          ,
                                              &cStatus                )
    if cStatus != CStatus_OK: raise QCModelError ( "Error constructing Fock matrix from Coulomb fit integrals." )
    return eFit

def FockConstruction_MakeFromFitIntegralsNonCoulomb ( SymmetricMatrix dTotal           not None ,
                                                      BlockStorage    fitIntegrals     not None ,
                                                      SymmetricMatrix fitMatrix        not None ,
                                                      SymmetricMatrix fitCoulombMatrix not None ,
                                                      CReal           totalCharge               ,
                                                      RealArray1D     fitCoefficients  not None ,
                                                      RealArray1D     fitVectorD       not None ,
                                                      SymmetricMatrix fTotal           not None ):
    """Coulomb Fock matrix from non-Coulomb fit integrals."""
    cdef CReal   eFit
    cdef CStatus cStatus = CStatus_OK
    eFit = Fock_MakeFromFitIntegralsNonCoulomb ( fitIntegrals.cObject     ,
                                                 fitMatrix.cObject        ,
                                                 fitCoulombMatrix.cObject ,
                                                 totalCharge              ,
                                                 fitCoefficients.cObject  ,
                                                 fitVectorD.cObject       ,
                                                 dTotal.cObject           ,
                                                 fTotal.cObject           ,
                                                 &cStatus                 )
    if cStatus != CStatus_OK: raise QCModelError ( "Error constructing Fock matrix from non-Coulomb fit integrals." )
    return eFit

# . All the following initialize fTotal and fSpin if appropriate.
# . Unfortunately they need to do this as they calculate the energy from the Fock matrices and also modify them.
# . fTotal/fSpin.
def FockConstruction_MakeFromTEIs ( SymmetricMatrix dTotal not None               ,
                                    SymmetricMatrix dSpin                         ,
                                    BlockStorage    twoElectronIntegrals not None ,
                                    CReal           exchangeScaling               ,
                                    SymmetricMatrix fTotal not None               ,
                                    SymmetricMatrix fSpin                         ):
    """Fock matrices from TEIs."""
    cdef CReal             eTEI
    cdef CSymmetricMatrix *cDSpin = NULL
    cdef CSymmetricMatrix *cFSpin = NULL
    if dSpin is not None: cDSpin = dSpin.cObject
    if fSpin is not None: cFSpin = fSpin.cObject
    eTEI = Fock_MakeFromTEIs ( twoElectronIntegrals.cObject ,
                               dTotal.cObject               ,
                               cDSpin                       ,
                               exchangeScaling              ,
                               fTotal.cObject               ,
                               cFSpin                       )
    return eTEI

# . fTotal.
def FockConstruction_MakeFromTEIsCoulomb ( SymmetricMatrix dTotal not None               ,
                                           BlockStorage    twoElectronIntegrals not None ,
                                           SymmetricMatrix fTotal not None               ):
    """Coulomb Fock matrix from TEIs."""
    cdef CReal             eTEI
    eTEI = Fock_MakeFromTEIsCoulomb ( twoElectronIntegrals.cObject ,
                                      dTotal.cObject               ,
                                      fTotal.cObject               )
    return eTEI

# . fSpin.
def FockConstruction_MakeFromTEIsExchange ( SymmetricMatrix dTotal not None               ,
                                            SymmetricMatrix dSpin                         ,
                                            BlockStorage    twoElectronIntegrals not None ,
                                            CReal           exchangeScaling               ,
                                            SymmetricMatrix fTotal not None               ,
                                            SymmetricMatrix fSpin                         ):
    """Exchange Fock matrices from TEIs."""
    cdef CReal             eTEI
    cdef CSymmetricMatrix *cDSpin = NULL
    cdef CSymmetricMatrix *cFSpin = NULL
    if dSpin is not None: cDSpin = dSpin.cObject
    if fSpin is not None: cFSpin = fSpin.cObject
    eTEI = Fock_MakeFromTEIsExchange ( twoElectronIntegrals.cObject ,
                                       dTotal.cObject               ,
                                       cDSpin                       ,
                                       exchangeScaling              ,
                                       fTotal.cObject               ,
                                       cFSpin                       )
    return eTEI
