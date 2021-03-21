"""CI CPHF functions."""

from .QCModelError import QCModelError

#===================================================================================================================================
# . Functions.
#===================================================================================================================================
def CICPHF_ApplyCPHFMatrix (                 n1                            ,
                             IntegerArray2D  in1                  not None ,
                                             n2                            ,
                             IntegerArray2D  in2                  not None ,
                             RealArray1D     aDiagonal                     ,
                             RealArray1D     b                    not None ,
                             RealArray2D     orbitals             not None ,
                             BlockStorage    twoElectronIntegrals not None ,
                             SymmetricMatrix work1                not None ,
                             SymmetricMatrix work2                not None ,
                             RealArray1D     x                    not None ):
    """Apply the CPHF matrix to a vector."""
    cdef CRealArray1D *cADiagonal = NULL
    if aDiagonal is not None: cADiagonal = aDiagonal.cObject
    CCICPHF_ApplyCPHFMatrix ( n1                           ,
                              in1.cObject                  ,
                              n2                           ,
                              in2.cObject                  ,
                              cADiagonal                   ,
                              b.cObject                    ,
                              orbitals.cObject             ,
                              twoElectronIntegrals.cObject ,
                              work1.cObject                ,
                              work2.cObject                ,
                              x.cObject                    )

def CICPHF_CalculateCPHFVectors ( nActive   ,
                                  nCore     ,
                                  nOrbitals ,
                                  BlockStorage          twoElectronIntegrals      not None ,
                                  DoubleSymmetricMatrix twoPDM                    not None ,
                                  RealArray1D           energies                  not None ,
                                  RealArray1D           occupancies               not None ,
                                  RealArray2D           orbitals                  not None ,
                                  RealArrayND           moTEI234                  not None ,
                                  SymmetricMatrix       fCore                     not None ,
                                  SymmetricMatrix       onePDM                    not None ,
                                  SymmetricMatrix       onePDMMO                  not None ,
                                  SymmetricMatrix       work1                     not None ,
                                  SymmetricMatrix       work2                     not None ,
                                                        numberDegenerateRedundant          ,
                                                        numberNonRedundant                 ,
                                                        numberRedundant                    ,
                                  IntegerArray2D        indicesNR                 not None ,
                                  IntegerArray2D        indicesR                  not None ,
                                  RealArray1D           aDiagonal                 not None ,
                                  RealArray1D           qNR                       not None ,
                                  RealArray1D           qR                        not None ,
                                  RealArray1D           preconditioner            not None ):


    """Calculate the CPHF vectors."""
    cdef CInteger  cDR = 0
    cdef CInteger  cNR = 0
    cdef CInteger  cR  = 0
    cdef CStatus   cStatus = CStatus_OK
    cDR = numberDegenerateRedundant
    cNR = numberNonRedundant
    cR  = numberRedundant
    CCICPHF_CalculateCPHFVectors ( nActive                       ,
                                   nCore                         ,
                                   nOrbitals                     ,
                                   twoElectronIntegrals.cObject  ,
                                   twoPDM.cObject                ,
                                   energies.cObject              ,
                                   occupancies.cObject           ,
                                   orbitals.cObject              ,
                                   moTEI234.cObject              ,
                                   fCore.cObject                 ,
                                   onePDM.cObject                ,
                                   onePDMMO.cObject              ,
                                   work1.cObject                 ,
                                   work2.cObject                 ,
                                   &cDR                          ,
                                   &cNR                          ,
                                   &cR                           ,
                                   indicesNR.cObject             ,
                                   indicesR.cObject              ,
                                   aDiagonal.cObject             ,
                                   qNR.cObject                   ,
                                   qR.cObject                    ,
                                   preconditioner.cObject        ,
                                   &cStatus                      )
    if cStatus != CStatus_OK: raise QCModelError ( "Error calculating CPHF vectors." )
    return ( cDR, cR )

def CICPHF_Transform (                 n1                ,
                       IntegerArray2D  in1      not None ,
                       RealArray1D     x1       not None ,
                                       n2                ,
                       IntegerArray2D  in2      not None ,
                       RealArray1D     x2       not None ,
                       RealArray2D     orbitals not None ,
                                       doScale           ,
                       SymmetricMatrix work     not None ,
                       SymmetricMatrix z        not None ):
    """Transform the CPHF vectors."""
    cdef CBoolean cDoScale
    if doScale: cDoScale = CTrue
    else:       cDoScale = CFalse
    CCICPHF_Transform ( n1               ,
                        in1.cObject      ,
                        x1.cObject       ,
                        n2               ,
                        in2.cObject      ,
                        x2.cObject       ,
                        orbitals.cObject ,
                        cDoScale         ,
                        work.cObject     ,
                        z.cObject        )
