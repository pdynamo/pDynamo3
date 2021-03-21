"""Loewdin multipole evaluator."""

from .QCDefinitions import BasisRepresentation
from .QCModelError  import QCModelError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class LoewdinMultipoleEvaluator:

    def BondOrders ( self                                ,
                     target                              ,
                     SymmetricMatrix density    not None ,
                     SymmetricMatrix bondOrders not None ):
        """Bond orders."""
        cdef IntegerArray1D basisIndices
        cdef RealArray2D    loewdinT
        cdef CStatus        cStatus = CStatus_OK
        basisIndices = target.qcState.orbitalBases._centerFunctionPointers[self.basisRepresentation]
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
        cdef IntegerArray1D basisIndices
        cdef RealArray2D    loewdinT
        if multipoleOrder != 0: raise ValueError ( "Invalid multipole order." )
        basisIndices = target.qcState.orbitalBases._centerFunctionPointers[self.basisRepresentation]
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
        cdef RealArray2D     loewdinT
        cdef SymmetricMatrix onePDM
        if multipoleOrder != 0: raise ValueError ( "Invalid multipole order." )
        basisIndices = target.qcState.orbitalBases._centerFunctionPointers[self.basisRepresentation]
        loewdinT     = target.scratch.loewdinT
        if density is None: onePDM = target.scratch.onePDMP.density
        else:               onePDM = density
        if withNuclei: target.qcState.nuclearCharges.CopyTo ( multipoles )
        else:          multipoles.Set ( 0.0 )
        Loewdin_AtomicCharges ( basisIndices.cObject ,
                                loewdinT.cObject     ,
                                onePDM.cObject       ,
                                multipoles.cObject   )

    def WeightedDensity ( self                           ,
                          target                         , 
                          multipoleOrder                 ,
                          RealArray1D     dXdQ not None  ,
                          SymmetricMatrix wDM  not None  ,
                          SymmetricMatrix density = None ):
        """Weighted density contributions."""
        cdef IntegerArray1D  basisIndices
        cdef RealArray1D     eigenValues
        cdef RealArray2D     a2w
        cdef RealArray2D     eigenVectors
        cdef RealArray2D     loewdinT
        cdef RealArray2D     w2a
        cdef SymmetricMatrix onePDM
        cdef CStatus         cStatus = CStatus_OK
        if multipoleOrder != 0: raise ValueError ( "Invalid multipole order." )
        basisIndices = target.qcState.orbitalBases._centerFunctionPointers[self.basisRepresentation]
        w2a          = target.qcState.orbitalBases.w2a
        a2w          = target.qcState.orbitalBases.a2w
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
                                  a2w.cObject          ,
                                  w2a.cObject          ,
                                  onePDM.cObject       ,
                                  NULL                 ,
                                  wDM.cObject          ,
                                  &cStatus             )
        if cStatus != CStatus_OK: raise QCModelError ( "Error calculating weighted density matrix." )

    @property
    def basisRepresentation ( self ): return BasisRepresentation.Actual
