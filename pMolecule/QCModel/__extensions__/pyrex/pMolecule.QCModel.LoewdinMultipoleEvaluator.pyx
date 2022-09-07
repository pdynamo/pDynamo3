"""Loewdin multipole evaluator."""

from  pScientific.Arrays import Array, StorageType
from .QCModelError       import QCModelError

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

    def ChargeRestraintMatrix ( self, target, restraint ):
        """Get the W matrix and the core term appropriate for the charge restraint."""
        cdef IntegerArray1D  basisIndices
        cdef IntegerArray1D  crIndices
        cdef RealArray1D     crWeights
        cdef RealArray1D     nuclearCharges
        cdef SymmetricMatrix loewdinT
        cdef SymmetricMatrix W
        cdef CBoolean        cIsSpin
        cdef CReal           core
        basisIndices   = target.qcState.orbitalBases.centerFunctionPointers
        nuclearCharges = target.qcState.nuclearCharges
        crIndices      = restraint.indices
        crWeights      = restraint.weights
        loewdinT       = target.scratch.loewdinT
        W              = Array.WithExtent ( len ( target.qcState.orbitalBases ), storageType = StorageType.Symmetric )
        if restraint.isSpin: cIsSpin = CTrue
        else:                cIsSpin = CFalse
        core = Loewdin_ChargeRestraintMatrix ( basisIndices.cObject   ,
                                               nuclearCharges.cObject ,
                                               crIndices.cObject      ,
                                               crWeights.cObject      ,
                                               cIsSpin                ,
                                               loewdinT.cObject       ,
                                               W.cObject              )
        return ( W, core )

    def ChargeRestraintWeightedDensity ( self                             ,
                                         target                           ,
                                         restraints                       ,
                                         dict            Lambdas not None ,
                                         SymmetricMatrix wDM     not None ):
        """Get the weighted density matrix for the charge restraints."""
        cdef IntegerArray1D  basisIndices
        cdef IntegerArray1D  crIndices
        cdef RealArray1D     crWeights
        cdef RealArray2D     eigenVectors
        cdef RealArray2D     Z
        cdef SymmetricMatrix A
        cdef CBoolean        cIsSpin
        cdef CReal           cdRdL
        cdef CRealArray2D   *cZ
        # . Gather data.
        basisIndices    = target.qcState.orbitalBases.centerFunctionPointers
        eigenValueRoots = target.scratch.overlapEigenValueRoots
        eigenVectors    = target.scratch.overlapEigenVectors
        loewdinT        = target.scratch.loewdinT
        # . Find out whether there are charge and spin restraints.
        spins     = [ restraint.isSpin for restraint in restraints.values ( ) ]
        hasCharge = not all ( spins )
        hasSpin   =     any ( spins )
        # . Get Z matrices - wasteful as only need a fraction of these matrix elements!
        PX        = Array.WithShape ( loewdinT.shape )
        zMatrices = {}
        for ( flag, attribute, tag ) in ( ( hasCharge, "onePDMP", "ZP" ) ,
                                          ( hasSpin  , "onePDMQ", "ZQ" ) ):
            if flag:
                Z = Array.WithShape ( loewdinT.shape )
                P = target.scratch.Get ( attribute ).density
                P.MatrixMultiply  ( loewdinT, PX )
                Z.MatrixMultiply ( PX, eigenVectors, xTranspose = True )
                zMatrices[tag] = Z
        # . Get the partial WDM by looping over the restraints.
        A = Array.WithExtent ( loewdinT.rows, storageType = StorageType.Symmetric )
        A.Set ( 0.0 )
        for key in sorted ( restraints.keys ( ) ):
            restraint    = restraints[key]
            crIndices    = restraint.indices
            crWeights    = restraint.weights
            ( _, cdRdL ) = Lambdas[key]
            if restraint.isSpin:
                Z       = zMatrices["ZQ"]
                cZ      = Z.cObject
                cIsSpin = CTrue
            else:
                Z       = zMatrices["ZP"]
                cZ      = Z.cObject
                cIsSpin = CFalse
            Loewdin_ChargeRestraintWeightedDensity ( basisIndices.cObject ,
                                                     crIndices.cObject    ,
                                                     crWeights.cObject    ,
                                                     cIsSpin              ,
                                                     cdRdL                ,
                                                     eigenVectors.cObject ,
                                                     cZ                   ,
                                                     A.cObject            )
        # . Scale the partial WDM by the eigenvalue factors.
        for i in range ( eigenValueRoots.extent ):
            eI = eigenValueRoots[i]
            for j in range ( i+1 ):
                eJ = eigenValueRoots[j]
                A[i,j] /= ( eI + eJ )
        # . Add into the full matrix the partial WDM transformed by V * A * V^T.
        B = Array.WithExtent ( loewdinT.rows, storageType = StorageType.Symmetric )
        A.Transform ( eigenVectors, B, useTranspose = True )
        wDM.Add ( B )
        # . When rationalize Arrays can use the following instead which doesn't require B.
        #wDM.Transform ( eigenVectors, A, alpha = 1.0, beta = 1.0, useTranspose = True )

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
