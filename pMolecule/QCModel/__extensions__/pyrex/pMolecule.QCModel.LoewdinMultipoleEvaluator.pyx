"""Loewdin multipole evaluator."""

from .QCModelError import QCModelError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class LoewdinMultipoleEvaluator:

    def BondOrders ( self                                ,
                     target                              ,
                     SymmetricMatrix density    not None ,
                     SymmetricMatrix bondOrders not None ):
        """Bond orders."""
        cdef IntegerArray1D  basisIndices
        cdef SymmetricMatrix loewdinT
        cdef CStatus         cStatus = CStatus_OK
        basisIndices = target.qcState.orbitalBases.centerFunctionPointers
        loewdinT     = target.scratch.loewdinT
        Loewdin_BondOrders ( basisIndices.cObject ,
                             loewdinT.cObject     ,
                             density.cObject      ,
                             bondOrders.cObject   ,
                             &cStatus             )
        if cStatus != CStatus_OK: raise QCModelError ( "Error calculating bond orders." )

    def FockMultipoleDerivatives ( self                          ,
                                   target                        , 
                                   multipoleOrder                ,
                                   RealArray1D     dXdQ not None ,
                                   SymmetricMatrix dXdP not None ):
        """Multipole density derivatives."""
        cdef IntegerArray1D  basisIndices
        cdef SymmetricMatrix loewdinT
        if multipoleOrder != 0: raise ValueError ( "Invalid multipole order." )
        basisIndices = target.qcState.orbitalBases.centerFunctionPointers
        loewdinT     = target.scratch.loewdinT
        Loewdin_ChargeDensityDerivatives ( basisIndices.cObject ,
                                           dXdQ.cObject         ,
                                           loewdinT.cObject     ,
                                           dXdP.cObject         )

    def FockMultipoles ( self                                ,
                         target                              ,
                         multipoleOrder                      ,
                         RealArray1D     multipoles not None ,
                         SymmetricMatrix density    = None   ,
                                         withNuclei = True   ):
        """Multipoles for interactions with external potentials."""
        cdef IntegerArray1D  basisIndices
        cdef SymmetricMatrix loewdinT
        cdef SymmetricMatrix onePDM
        cdef CStatus         cStatus = CStatus_OK
        if multipoleOrder != 0: raise ValueError ( "Invalid multipole order." )
        basisIndices = target.qcState.orbitalBases.centerFunctionPointers
        loewdinT     = target.scratch.loewdinT
        if density is None: onePDM = target.scratch.onePDMP.density
        else:               onePDM = density
        if withNuclei: target.qcState.nuclearCharges.CopyTo ( multipoles )
        else:          multipoles.Set ( 0.0 )
        Loewdin_AtomicCharges ( basisIndices.cObject ,
                                loewdinT.cObject     ,
                                onePDM.cObject       ,
                                multipoles.cObject   ,
                                &cStatus             )
        if cStatus != CStatus_OK: raise QCModelError ( "Error calculating Fock multipoles." )

    def WeightedDensity ( self                           ,
                          target                         , 
                          multipoleOrder                 ,
                          RealArray1D     dXdQ not None  ,
                          SymmetricMatrix wDM  not None  ,
                          SymmetricMatrix density = None ):
        """Weighted density contributions."""
        cdef IntegerArray1D  basisIndices
        cdef RealArray1D     eigenValues
        cdef RealArray2D     eigenVectors
        cdef SymmetricMatrix loewdinT
        cdef SymmetricMatrix onePDM
        cdef CStatus         cStatus = CStatus_OK
        if multipoleOrder != 0: raise ValueError ( "Invalid multipole order." )
        basisIndices = target.qcState.orbitalBases.centerFunctionPointers
        eigenValues  = target.scratch.overlapEigenValues
        eigenVectors = target.scratch.overlapEigenVectors
        loewdinT     = target.scratch.loewdinT
        if density is None: onePDM = target.scratch.onePDMP.density
        else:               onePDM = density
        Loewdin_WeightedDensity ( basisIndices.cObject ,
                                  dXdQ.cObject         ,
                                  eigenValues.cObject  ,
                                  eigenVectors.cObject ,
                                  loewdinT.cObject     ,
                                  onePDM.cObject       ,
                                  NULL                 ,
                                  wDM.cObject          ,
                                  &cStatus             )
        if cStatus != CStatus_OK: raise QCModelError ( "Error calculating weighted density matrix." )
