"""Mulliken multipole evaluator."""

from .QCDefinitions import BasisRepresentation
from .QCModelError  import QCModelError

#===================================================================================================================================
# . Class.
#===================================================================================================================================
cdef class MullikenMultipoleEvaluator:

    def BondOrders ( self                                ,
                     target                              ,
                     SymmetricMatrix density    not None ,
                     SymmetricMatrix bondOrders not None ):
        """Bond orders."""
        cdef IntegerArray1D  basisIndices
        cdef SymmetricMatrix overlap
        cdef CStatus         cStatus = CStatus_OK
        basisIndices = target.qcState.orbitalBases._centerFunctionPointers[self.basisRepresentation]
        overlap      = target.scratch.overlapMatrix
        Mulliken_BondOrders ( basisIndices.cObject ,
                              density.cObject      ,
                              overlap.cObject      ,
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
        cdef SymmetricMatrix overlap
        if multipoleOrder != 0: raise ValueError ( "Invalid multipole order." )
        basisIndices = target.qcState.orbitalBases._centerFunctionPointers[self.basisRepresentation]
        overlap      = target.scratch.overlapMatrix
        Mulliken_ChargeDensityDerivatives ( basisIndices.cObject ,
                                            dXdQ.cObject         ,
                                            overlap.cObject      ,
                                            dXdP.cObject         )

    def FockMultipoles ( self                                ,
                         target                              ,
                         multipoleOrder                      ,
                         RealArray1D     multipoles not None ,
                         SymmetricMatrix density    = None   ,
                                         withNuclei = True   ):
        """Multipoles for interactions with external potentials."""
        # , Multipoles set here.
        cdef IntegerArray1D  basisIndices
        cdef SymmetricMatrix onePDM
        cdef SymmetricMatrix overlap
        if multipoleOrder != 0: raise ValueError ( "Invalid multipole order." )
        basisIndices   = target.qcState.orbitalBases._centerFunctionPointers[self.basisRepresentation]
        overlap        = target.scratch.overlapMatrix
        if density is None: onePDM = target.scratch.onePDMP.density
        else:               onePDM = density
        if withNuclei: target.qcState.nuclearCharges.CopyTo ( multipoles )
        else:          multipoles.Set ( 0.0 )
        Mulliken_AtomicCharges ( basisIndices.cObject ,
                                 onePDM.cObject       ,
                                 overlap.cObject      ,
                                 multipoles.cObject   )

    def WeightedDensity ( self                           ,
                          target                         , 
                          multipoleOrder                 ,
                          RealArray1D     dXdQ not None  ,
                          SymmetricMatrix wDM  not None  ,
                          SymmetricMatrix density = None ):
        """Weighted density contributions."""
        cdef IntegerArray1D  basisIndices
        cdef SymmetricMatrix onePDM
        if multipoleOrder != 0: raise ValueError ( "Invalid multipole order." )
        basisIndices = target.qcState.orbitalBases._centerFunctionPointers[self.basisRepresentation]
        if density is None: onePDM = target.scratch.onePDMP.density
        else:               onePDM = density
        Mulliken_WeightedDensity ( basisIndices.cObject ,
                                   dXdQ.cObject         ,
                                   onePDM.cObject       ,
                                   wDM.cObject          )

    @property
    def basisRepresentation ( self ): return BasisRepresentation.Work
