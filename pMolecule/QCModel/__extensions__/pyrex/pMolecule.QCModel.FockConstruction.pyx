"""Fock construction functions."""

from .QCModelError import QCModelError

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
# . dTotal is temporarily modified. This needs to be changed.
# . fTotal is not initialized.
def FockConstruction_MakeFromFitIntegrals ( SymmetricMatrix dTotal       not None ,
                                            BlockStorage    fitIntegrals not None ,
                                            SymmetricMatrix fitMatrix    not None ,
                                            CReal           totalCharge           ,
                                            RealArray1D     fitPotential not None ,
                                            SymmetricMatrix fTotal       not None ):
    """Coulomb Fock matrix from fit integrals."""
    cdef CReal   eTEI
    cdef CStatus cStatus = CStatus_OK
    eTEI = Fock_MakeFromFitIntegrals ( fitIntegrals.cObject ,
                                       fitMatrix.cObject    ,
                                       totalCharge          ,
                                       fitPotential.cObject ,
                                       dTotal.cObject       ,
                                       fTotal.cObject       ,
                                       &cStatus             )
    if cStatus != CStatus_OK: raise QCModelError ( "Error constructing Coulomb Fock matrix from fit integrals." )
    return eTEI

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
