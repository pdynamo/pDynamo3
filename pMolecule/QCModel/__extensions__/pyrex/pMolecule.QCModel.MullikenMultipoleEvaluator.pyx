"""Mulliken multipole evaluator."""

from  pScientific.Arrays import Array, StorageType
from .QCModelError       import QCModelError

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
        basisIndices = target.qcState.orbitalBases.centerFunctionPointers
        overlap      = target.scratch.overlapMatrix
        Mulliken_BondOrders ( basisIndices.cObject ,
                              density.cObject      ,
                              overlap.cObject      ,
                              bondOrders.cObject   ,
                              &cStatus             )
        if cStatus != CStatus_OK: raise QCModelError ( "Error calculating bond orders." )

    def ChargeRestraintMatrix ( self, target, restraint ):
        """Get the W matrix and the core term appropriate for the charge restraint."""
        cdef IntegerArray1D  basisIndices
        cdef IntegerArray1D  crIndices
        cdef RealArray1D     crWeights
        cdef RealArray1D     nuclearCharges
        cdef SymmetricMatrix overlap
        cdef SymmetricMatrix W
        cdef CBoolean        cIsSpin
        cdef CReal           core
        basisIndices   = target.qcState.orbitalBases.centerFunctionPointers
        nuclearCharges = target.qcState.nuclearCharges
        crIndices      = restraint.indices
        crWeights      = restraint.weights
        overlap        = target.scratch.overlapMatrix
        W              = Array.WithExtent ( len ( target.qcState.orbitalBases ), storageType = StorageType.Symmetric )
        if restraint.isSpin: cIsSpin = CTrue
        else:                cIsSpin = CFalse
        core = Mulliken_ChargeRestraintMatrix ( basisIndices.cObject   ,
                                                nuclearCharges.cObject ,
                                                crIndices.cObject      ,
                                                crWeights.cObject      ,
                                                cIsSpin                ,
                                                overlap.cObject        ,
                                                W.cObject              )
        return ( W, core )

    def ChargeRestraintWeightedDensity ( self                             ,
                                         target                           ,
                                         restraints                       ,
                                         dict            Lambdas not None ,
                                         SymmetricMatrix wDM     not None ):
        """Get the weighted density matrix for the charge restraints."""
        cdef IntegerArray1D    basisIndices
        cdef IntegerArray1D    crIndices
        cdef RealArray1D       crWeights
        cdef SymmetricMatrix   density
        cdef CBoolean          cIsSpin
        cdef CReal             cdRdL
        cdef CSymmetricMatrix *cDensity
        cdef CSymmetricMatrix *cOnePDMP
        cdef CSymmetricMatrix *cOnePDMQ
        basisIndices = target.qcState.orbitalBases.centerFunctionPointers
        onePDM       = target.scratch.Get ( "onePDMP" )
        if onePDM is None:
            cOnePDMP = NULL
        else:
            density  = onePDM.density
            cOnePDMP = density.cObject
        onePDM = target.scratch.Get ( "onePDMQ" )
        if onePDM is None:
            cOnePDMQ = NULL
        else:
            density  = onePDM.density
            cOnePDMQ = density.cObject
        for key in sorted ( restraints.keys ( ) ):
            restraint    = restraints[key]
            crIndices    = restraint.indices
            crWeights    = restraint.weights
            ( _, cdRdL ) = Lambdas[key]
            if restraint.isSpin:
                cDensity = cOnePDMQ
                cIsSpin  = CTrue
            else:
                cDensity = cOnePDMP
                cIsSpin  = CFalse
            Mulliken_ChargeRestraintWeightedDensity ( basisIndices.cObject ,
                                                      crIndices.cObject    ,
                                                      crWeights.cObject    ,
                                                      cIsSpin              ,
                                                      cdRdL                ,
                                                      cDensity             ,
                                                      wDM.cObject          )

    def FockMultipoleDerivatives ( self                          ,
                                   target                        , 
                                   multipoleOrder                ,
                                   RealArray1D     dXdQ not None ,
                                   SymmetricMatrix dXdP not None ):
        """Multipole density derivatives."""
        cdef IntegerArray1D  basisIndices
        cdef SymmetricMatrix overlap
        if multipoleOrder != 0: raise ValueError ( "Invalid multipole order." )
        basisIndices = target.qcState.orbitalBases.centerFunctionPointers
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
        basisIndices   = target.qcState.orbitalBases.centerFunctionPointers
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
        basisIndices = target.qcState.orbitalBases.centerFunctionPointers
        if density is None: onePDM = target.scratch.onePDMP.density
        else:               onePDM = density
        Mulliken_WeightedDensity ( basisIndices.cObject ,
                                   dXdQ.cObject         ,
                                   onePDM.cObject       ,
                                   wDM.cObject          )
